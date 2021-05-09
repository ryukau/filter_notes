#include <cmath>

template<typename Sample> class TempoSynchronizer {
private:
  enum class State { free, steady, decelerating, accelerating };

  State state = State::free;
  Sample v1 = 0;
  Sample p1 = 0;
  Sample v2 = 0;
  Sample p2 = 0;
  Sample lastTempo = 0;
  Sample lastSync = 0;
  Sample lastElapsedBeats = 0;

  Sample midTime = 0; // In samples.
  Sample midVelocity = 0;
  Sample counter = 0;

  inline void setLastValues(Sample tempo, Sample sync, Sample elapsedBeats)
  {
    lastTempo = tempo;
    lastSync = sync;
    lastElapsedBeats = elapsedBeats;
  }

  inline Sample wrap(Sample value) { return value - std::floor(value); }

public:
  // Only use this when DAW is not playing and the synth resets LFO for each note-on.
  void reset()
  {
    p1 = 0;
    p2 = 0;
  }

  void reset(Sample sampleRate, Sample tempo, Sample sync)
  {
    v1 = std::fabs(sync) <= std::numeric_limits<Sample>::min()
      ? 0
      : tempo / (Sample(60) * sampleRate * sync);
    p1 = 0;
    v2 = v1;
    p2 = 0;
    lastTempo = tempo;
    lastSync = sync;
  }

  // Must call this method at the start of each DSP processing cycle.
  void
  prepare(Sample sampleRate, Sample tempo, Sample sync, Sample elapsedBeats, bool isFree)
  {
    if (std::fabs(sync) <= std::numeric_limits<Sample>::min()) {
      v1 = 0;
      setLastValues(tempo, sync, elapsedBeats);
      return;
    }

    v1 = tempo / (Sample(60) * sampleRate * sync);

    if (isFree) state = State::free;
    if (state == State::free) {
      if (!isFree) state = State::steady;
      setLastValues(tempo, sync, elapsedBeats);
      return;
    }

    auto p0 = p1;
    p1 = wrap(elapsedBeats / sync);

    if (lastTempo != tempo || lastSync != sync || lastElapsedBeats > elapsedBeats) {
      if (std::fabs(lastSync) <= std::numeric_limits<Sample>::min()) {
        v2 = 0;
        p2 = p1;
      } else if (lastElapsedBeats > elapsedBeats) {
        v2 = lastTempo / (Sample(60) * sampleRate * lastSync);
        p2 = p0;
      } else {
        v2 = lastTempo / (Sample(60) * sampleRate * lastSync);
        if (state == State::steady) p2 = wrap(elapsedBeats / lastSync);
      }

      state = State::decelerating;

      midTime = Sample(0.05) * sampleRate;
      auto distance = wrap(p1 + Sample(2) * v1 * midTime - p2);
      auto k = std::ceil((v1 + v2) * Sample(0.5) * midTime - distance);
      midVelocity = (distance + k) / midTime - (v1 + v2) * Sample(0.5);

      counter = 0;
    }

    setLastValues(tempo, sync, elapsedBeats);
  }

  Sample process()
  {
    Sample outPhase;
    switch (state) {
      default:
      case State::free:
      case State::steady: {
        outPhase = p1;
        p1 = wrap(p1 + v1);
      } break;

      case State::decelerating: {
        outPhase = p2;
        p2 = wrap(p2 + v2 + (midVelocity - v2) * counter / midTime);
        if (++counter >= midTime) {
          state = State::accelerating;
          counter = 0;
        }
      } break;

      case State::accelerating: {
        outPhase = p2;
        p2 = wrap(p2 + midVelocity + (v1 - midVelocity) * counter / midTime);
        if (++counter >= midTime) {
          state = State::steady;
          p1 = p2;
        }
      } break;
    }
    return outPhase;
  }
};

#include <array>
#include <iostream>

int main()
{
  float sampleRate = 48000.0f;
  float tempo = 120.0f;
  float sync = 0.25f;
  float elapsedBeats = 0.0f;

  std::array<float, 2048> buffer{};

  // Example usage of TempoSynchronizer.
  TempoSynchronizer<float> lfo;

  for (size_t idx = 0; idx < buffer.size();) {
    lfo.prepare(sampleRate, tempo, sync, elapsedBeats, false);

    auto beatsDelta = tempo / (60.0f * sampleRate);
    for (size_t j = 0; j < 512; ++j) {
      if (idx >= buffer.size()) break;
      buffer[idx++] = lfo.process();
      elapsedBeats += beatsDelta;

      std::cout << idx - 1 << ": " << buffer[idx - 1] << "\n";
    }
  }

  return 0;
}
