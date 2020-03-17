#include <algorithm>
#include <cstdint>

template<typename Sample> class LinearADSREnvelope {
public:
  void setup(Sample sampleRate) { this->sampleRate = sampleRate; }

  Sample adaptTime(Sample seconds, Sample noteFreq)
  {
    const Sample cycle = Sample(1) / noteFreq;
    return seconds >= cycle ? seconds : cycle > Sample(0.1) ? Sample(0.1) : cycle;
  }

  Sample secondToDelta(Sample seconds) { return Sample(1) / (sampleRate * seconds); }

  void trigger(
    Sample attackTime,
    Sample decayTime,
    Sample sustainLevel,
    Sample releaseTime,
    Sample noteFreq)
  {
    state = stateAttack;
    value = Sample(1);
    atkOffset = state == stateTerminated ? 0 : out;
    atkRange = Sample(1) - atkOffset;
    set(attackTime, decayTime, sustainLevel, releaseTime, noteFreq);
  }

  void set(
    Sample attackTime,
    Sample decayTime,
    Sample sustainLevel,
    Sample releaseTime,
    Sample noteFreq)
  {
    sus = std::clamp<Sample>(sustainLevel, Sample(0), Sample(1));
    atk = secondToDelta(adaptTime(attackTime, noteFreq));
    dec = secondToDelta(adaptTime(decayTime, noteFreq));
    rel = secondToDelta(adaptTime(releaseTime, noteFreq));
  }

  void release()
  {
    state = stateRelease;
    value = Sample(1);
    relRange = out;
  }

  Sample process()
  {
    if (value <= Sample(0)) {
      state = state + 1;
      value = Sample(1);
    }

    switch (state) {
      case stateAttack:
        value -= atk;
        out = atkOffset + atkRange * (Sample(1) - value);
        break;

      case stateDecay:
        value -= dec;
        out = (Sample(1) - sus) * value + sus;
        break;

      case stateSustain:
        out = sus;
        break;

      case stateRelease:
        value -= rel;
        out = relRange * value;
        break;

      default:
        return Sample(0);
    }
    return std::clamp<Sample>(out, Sample(0), Sample(1));
  }

protected:
  enum State : int32_t {
    stateAttack,
    stateDecay,
    stateSustain,
    stateRelease,
    stateTerminated
  };

  Sample sampleRate = 44100;
  Sample sus = 0.5;
  Sample atk = 0.01;
  Sample atkOffset = 0;
  Sample atkRange = 1;
  Sample dec = 0.01;
  Sample rel = 0.01;
  Sample relRange = 0.5;
  Sample value = 0;
  Sample out = 0;
  int32_t state = stateTerminated;
};
