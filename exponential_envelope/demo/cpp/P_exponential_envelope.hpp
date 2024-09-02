#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

// Exponential moving average filger.
template<typename Sample> class EmaLowpass {
public:
  // 原則として `double` で呼び出すこと。 `float` では精度が足りないので `timeInSamples`
  // が 5000 を超えるあたりで正しい値が出ない。
  static Sample samplesToP(Sample timeInSamples)
  {
    auto omega_c = std::numbers::pi_v<Sample> / timeInSamples;
    auto y = Sample(1) - cos(omega_c);
    return -y + sqrt((y + Sample(2)) * y);
  }

  static Sample thresholdToP(Sample timeInSample, Sample threshold = Sample(0.9))
  {
    return Sample(1)
      - std::pow(Sample(1) - threshold, Sample(1) / (timeInSample + Sample(1)));
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
    tailLength = int(0.01 * sampleRate);
  }

  void reset(Sample attackTime, Sample decayTime, Sample sustainLevel, Sample releaseTime)
  {
    state = State::attack;
    atk = int(sampleRate * attackTime);
    decTime = decayTime;
    sustain = std::clamp<Sample>(sustainLevel, Sample(0), Sample(1));
    relTime = releaseTime;
    lowpass.setP(EmaLowpass<double>::samplesToP(sampleRate * attackTime));
  }

  void set(Sample attackTime, Sample decayTime, Sample sustainLevel, Sample releaseTime)
  {
    switch (state) {
      case State::attack:
        atk = int(sampleRate * attackTime);
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
      lowpass.setP(EmaLowpass<double>::samplesToP(sampleRate * attackTime));
    else if (state == State::decay)
      lowpass.setP(EmaLowpass<double>::samplesToP(sampleRate * decayTime));
    else if (state == State::release)
      lowpass.setP(EmaLowpass<double>::samplesToP(sampleRate * releaseTime));
  }

  void release()
  {
    state = State::release;
    lowpass.setP(EmaLowpass<double>::samplesToP(sampleRate * relTime));
  }

  bool isAttacking() { return state == State::attack; }
  bool isReleasing() { return state == State::release; }
  bool isTerminated() { return state == State::terminated; }

  Sample process()
  {
    switch (state) {
      case State::attack: {
        value = lowpass.process(Sample(1));
        --atk;
        if (atk <= 0) {
          state = State::decay;
          lowpass.setP(EmaLowpass<double>::samplesToP(sampleRate * decTime));
        }
      } break;

      case State::decay:
        value = lowpass.process(sustain);
        break;

      case State::release:
        value = lowpass.process(0);
        if (value < threshold) {
          value = threshold;
          state = State::tail;
          tailCounter = tailLength;
        }
        break;

      case State::tail:
        --tailCounter;
        value = threshold * tailCounter / float(tailLength);
        if (tailCounter <= 0) {
          state = State::terminated;
          lowpass.reset(0);
        } else {
          lowpass.reset(value);
        }
        break;

      default:
        return 0;
    }
    return value;
  }

private:
  enum class State { attack, decay, release, tail, terminated };
  const Sample threshold = 1e-5;

  int tailLength = 32;
  int tailCounter = tailLength;

  EmaLowpass<Sample> lowpass;
  State state = State::terminated;
  int atk = 0;
  Sample decTime = 0;
  Sample relTime = 0;
  Sample sampleRate = 44100;
  Sample sustain = 1;
  Sample value = 0;
};
