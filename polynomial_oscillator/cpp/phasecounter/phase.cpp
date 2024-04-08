/*
This code is written to verify various phase counter implementations. These are
probably not suitable for real use.
*/

#include <algorithm>
#include <array>
#include <cinttypes>
#include <cmath>
#include <format>
#include <fstream>
#include <limits>
#include <string>

template<typename Sample> class PhaseFloor {
public:
  Sample phase = 0;
  Sample process(Sample sampleRateHz, Sample frequencyHz)
  {
    phase += frequencyHz / sampleRateHz;
    return phase -= std::floor(phase);
  }
};

template<typename Sample> class PhaseFmodIf {
public:
  Sample phase = 0;
  Sample process(Sample sampleRateHz, Sample frequencyHz)
  {
    constexpr Sample upperLimit = Sample(1);
    phase = std::fmod(phase + frequencyHz / sampleRateHz, upperLimit);
    if (frequencyHz < 0) phase = -phase;
    return phase;
  }
};

template<typename Sample> class PhaseFmodCopysign {
public:
  Sample phase = 0;
  Sample process(Sample sampleRateHz, Sample frequencyHz)
  {
    constexpr Sample upperLimit = Sample(1);
    phase = std::fmod(phase + frequencyHz / sampleRateHz, upperLimit);
    return phase *= std::copysign(Sample(1), frequencyHz);
  }
};

template<typename Sample> class PhaseIf {
public:
  Sample phase = 0;
  Sample process(Sample sampleRateHz, Sample frequencyHz)
  {
    constexpr Sample upperLimit = Sample(1);
    phase += std::clamp(frequencyHz / sampleRateHz, Sample(0), Sample(0.5));
    if (phase >= upperLimit) phase -= upperLimit;
    return phase;
  }
};

template<typename Sample> class PhaseWhile {
public:
  Sample phase = 0;
  Sample process(Sample sampleRateHz, Sample frequencyHz)
  {
    constexpr Sample upperLimit = Sample(1);
    phase += std::max(frequencyHz / sampleRateHz, Sample(0));
    while (phase >= upperLimit) phase -= upperLimit;
    return phase;
  }
};

template<typename Sample> class PhaseUnsigned {
public:
  unsigned phase = 0;
  Sample process(Sample sampleRateHz, Sample frequencyHz)
  {
    constexpr auto scaler = Sample(std::numeric_limits<unsigned>::max());
    unsigned frequencyInt = unsigned(scaler * frequencyHz / sampleRateHz);
    phase += frequencyInt;
    return Sample(phase) / scaler;
  }
};

template<typename Sample> class PhaseBitmask {
public:
  unsigned phase = 0;
  Sample process(Sample sampleRateHz, Sample frequencyHz)
  {
    constexpr auto scaler = unsigned(0xfffff);
    unsigned frequencyInt = unsigned(Sample(scaler) * frequencyHz / sampleRateHz);
    phase += frequencyInt;
    phase &= scaler;
    return Sample(phase) / Sample(scaler);
  }
};

template<typename Sample> class PhaseFixedPoint {
public:
  int32_t phase = 0;
  Sample process(Sample sampleRateHz, Sample frequencyHz)
  {
    constexpr auto fractionBits = int32_t(16);
    constexpr auto fractionMult = int32_t(1) << fractionBits;
    constexpr auto fractionMask = fractionMult - int32_t(1);
    int32_t sampleRateFixed = int32_t(sampleRateHz * Sample(fractionMult));
    int32_t frequencyFixed = int32_t(frequencyHz * Sample(fractionMult));

    auto div = [](int32_t x, int32_t y) -> int32_t {
      return int32_t((int64_t(x) * fractionMult) / y);
    };

    phase += div(frequencyFixed, sampleRateFixed);
    phase &= fractionMask;

    return Sample(phase) / Sample(fractionMult);
  }
};

template<typename Sample, typename PhaseCounter> auto run(std::string name)
{
  constexpr auto sampleRateHz = Sample(1000);
  constexpr auto frequencyHz = Sample(5);

  std::array<Sample, size_t(sampleRateHz)> output{};
  PhaseCounter phaser;
  for (auto &x : output) x = phaser.process(sampleRateHz, frequencyHz);

  std::string data = std::format("\"{}\":[", name);
  for (const auto &x : output) data += std::format("{},", x);
  data.pop_back();
  data += "],";
  return data;
}

int main()
{
  using Sample = float;

  std::string jsonText = "{";

#define RUN(phi) run<Sample, phi<Sample>>(#phi);
  jsonText += RUN(PhaseFloor);
  jsonText += RUN(PhaseFmodIf);
  jsonText += RUN(PhaseFmodCopysign);
  jsonText += RUN(PhaseIf);
  jsonText += RUN(PhaseWhile);
  jsonText += RUN(PhaseUnsigned);
  jsonText += RUN(PhaseBitmask);
  jsonText += RUN(PhaseFixedPoint);
#undef RUN

  jsonText.pop_back();
  jsonText += "}";

  std::ofstream os("phase.json");
  os << jsonText;

  return 0;
}
