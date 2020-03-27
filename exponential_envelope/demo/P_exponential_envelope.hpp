#pragma once

#include <algorithm>
#include <cmath>

constexpr double twopi = 6.283185307179586;

template<typename Sample> class PController {
public:
  // float 型での cutoffHz の下限は 3~4 Hz 程度。
  static Sample cutoffToP(Sample sampleRate, Sample cutoffHz)
  {
    auto omega_c = Sample(twopi) * cutoffHz / sampleRate;
    auto y = Sample(1) - cos(omega_c);
    return -y + sqrt((y + Sample(2)) * y);
  }

  void setP(Sample p) { kp = std::clamp<Sample>(p, Sample(0), Sample(1)); }
  void reset(Sample value = 0) { this->value = value; }
  Sample process(Sample input) { return value += kp * (input - value); }

  Sample kp = 1; // Range in [0, 1].
  Sample value = 0;
};

template<typename Sample> class ExpADSREnvelopeP {
public:
  void setup(Sample sampleRate)
  {
    this->sampleRate = sampleRate;
    tailLength = uint32_t(0.01 * sampleRate);
  }

  void reset(Sample attackTime, Sample decayTime, Sample sustainLevel, Sample releaseTime)
  {
    state = State::attack;
    sustain = std::clamp<Sample>(sustainLevel, Sample(0), Sample(1));
    atk = int32_t(sampleRate * attackTime);
    decTime = decayTime;
    relTime = releaseTime;
    pController.setP(PController<Sample>::cutoffToP(sampleRate, Sample(1) / attackTime));
  }

  void set(Sample attackTime, Sample decayTime, Sample sustainLevel, Sample releaseTime)
  {
    switch (state) {
      case State::attack:
        atk = int32_t(sampleRate * attackTime);
        // Fall through.

      case State::decay:
        decTime = decayTime;
        sustain = std::clamp<Sample>(sustainLevel, Sample(0), Sample(1));
        // Fall through.

      case State::release:
        relTime = releaseTime;

      default:
        break;
    }

    if (state == State::attack)
      pController.setP(
        PController<Sample>::cutoffToP(sampleRate, Sample(1) / attackTime));
    else if (state == State::decay)
      pController.setP(PController<Sample>::cutoffToP(sampleRate, Sample(1) / decayTime));
    else if (state == State::release)
      pController.setP(
        PController<Sample>::cutoffToP(sampleRate, Sample(1) / releaseTime));
  }

  void release()
  {
    state = State::release;
    pController.setP(PController<Sample>::cutoffToP(sampleRate, Sample(1) / relTime));
  }

  bool isAttacking() { return state == State::attack; }
  bool isReleasing() { return state == State::release; }
  bool isTerminated() { return state == State::terminated; }

  Sample process()
  {
    switch (state) {
      case State::attack: {
        value = pController.process(Sample(1));
        --atk;
        if (atk == 0) {
          state = State::decay;
          pController.setP(
            PController<Sample>::cutoffToP(sampleRate, Sample(1) / decTime));
        }
      } break;

      case State::decay:
        value = pController.process(sustain);
        break;

      case State::release:
        value = pController.process(0);
        if (value < threshold) {
          value = threshold;
          state = State::tail;
          tailCounter = tailLength;
        }
        break;

      case State::tail:
        --tailCounter;
        value = threshold * tailCounter / float(tailLength);
        if (tailCounter == 0) {
          state = State::terminated;
          pController.reset(0);
        } else {
          pController.reset(value);
        }
        break;

      default:
        return 0;
    }
    return value;
  }

private:
  enum class State : int32_t { attack, decay, release, tail, terminated };
  const Sample threshold = 1e-5;

  uint32_t tailLength = 32;
  uint32_t tailCounter = tailLength;

  PController<Sample> pController;
  State state = State::terminated;
  uint32_t atk = 0;
  Sample decTime = 0;
  Sample relTime = 0;
  Sample sampleRate = 44100;
  Sample sustain = 1;
  Sample value = 0;
};
