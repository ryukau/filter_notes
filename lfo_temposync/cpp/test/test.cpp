#include <cmath>

template<typename Sample> class TempoSynchronizer {
private:
  enum class State { steady, decelerating, accelerating };

  State state = State::steady;
  Sample v1 = 0;
  Sample p1 = 0;
  Sample v2 = 0;
  Sample p2 = 0;
  Sample lastTempo = 0;
  Sample lastSync = 0;

  Sample midTime = 0; // In samples.
  Sample midVelocity = 0;
  Sample counter = 0;

public:
  void reset(Sample sampleRate, Sample tempo, Sample sync)
  {
    v1 = tempo / (Sample(60) * sampleRate * sync);
    p1 = 0;
    v2 = v1;
    p2 = 0;
    lastTempo = tempo;
    lastSync = sync;
  }

  // Must call this method at the start of each DSP processing cycle.
  void prepare(Sample sampleRate, Sample tempo, Sample sync, Sample beatsElapsed)
  {
    v1 = tempo / (Sample(60) * sampleRate * sync);
    p1 = beatsElapsed / sync;
    p1 -= std::floor(p1);

    if (lastTempo != tempo || lastSync != sync) {
      state = State::decelerating;

      midTime = Sample(0.5) / v1;
      auto distance = p1 - p2;
      auto k = std::ceil((v1 + v2) * Sample(0.5) * midTime - distance);
      midVelocity = (distance + k) / midTime - (v1 + v2) * Sample(0.5);

      counter = 0;
    }
  }

  // Must call this method at the end of each DSP processing cycle.
  void postProcess(Sample tempo, Sample sync)
  {
    if (state == State::steady) {
      p2 = p1;
      v2 = v1;
    }
    lastTempo = tempo;
    lastSync = sync;
  }

  Sample process()
  {
    Sample outPhase;
    switch (state) {
      default:
      case State::steady: {
        outPhase = p1;

        p1 += v1;
        p1 -= std::floor(p1);
      } break;

      case State::decelerating: {
        outPhase = p2;

        p2 += v2 + (midVelocity - v2) * counter / midTime;
        p2 -= std::floor(p2);

        if (++counter >= midTime) {
          state = State::accelerating;
          counter = 0;
        }
      } break;

      case State::accelerating: {
        outPhase = p2;

        p2 += midVelocity + (v1 - midVelocity) * counter / midTime;
        p2 -= std::floor(p2);

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
  float beatsElapsed = 0.0f;

  std::array<float, 2048> buffer{};

  // Example usage of TempoSynchronizer.
  TempoSynchronizer<float> lfo;

  for (size_t idx = 0; idx < buffer.size();) {
    lfo.prepare(sampleRate, tempo, sync, beatsElapsed);

    auto beatsDelta = tempo / (60.0f * sampleRate);
    for (size_t j = 0; j < 512; ++j) {
      if (idx >= buffer.size()) break;
      buffer[idx++] = lfo.process();
      beatsElapsed += beatsDelta;
    }

    lfo.postProcess(tempo, sync);
  }

  for (const auto &v : buffer) std::cout << v << "\n";

  return 0;
}
