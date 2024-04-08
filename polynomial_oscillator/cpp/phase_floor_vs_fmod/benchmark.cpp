#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <string>
#include <vector>

template<typename Sample> class SineFloor {
public:
  std::string name{"SineFloor"};
  Sample phase = 0;
  Sample process(Sample normalizedFreq)
  {
    constexpr Sample twopi = Sample(2) * std::numbers::pi_v<Sample>;
    phase += normalizedFreq;
    phase -= std::floor(phase);
    return std::sin(twopi * phase);
  }
};

template<typename Sample> class SineFmod {
public:
  std::string name{"SineFmod"};
  Sample phase = 0;
  Sample process(Sample freqRadian)
  {
    constexpr Sample twopi = Sample(2) * std::numbers::pi_v<Sample>;
    phase = std::fmod(phase + freqRadian, twopi);
    return std::sin(phase);
  }
};

#define ROW_STR "|{:32}|{:20}|{:20}|\n"

template<typename Sample, typename Oscillator>
auto benchmark(Sample startFreq, Sample endFreq)
{
  constexpr Sample nSample = 100000;

  Oscillator osc;
  double sumElapsed = 0.0;
  std::vector<Sample> output(nSample);
  for (size_t i = 0; i < nSample; ++i) {
    auto freq = std::lerp(startFreq, endFreq, Sample(i) / Sample(nSample));

    auto start = std::chrono::steady_clock::now();

    output[i] = osc.process(freq);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }

  std::cout << std::format(ROW_STR, osc.name, sumElapsed, output.back());

  std::string data = std::format("\"{}\":[", osc.name);
  for (const auto &x : output) data += std::format("{},", x);
  data.pop_back();
  data += "],";
  return data;
}

int main()
{
  using Sample = float;
  constexpr Sample twopi = Sample(2) * std::numbers::pi_v<Sample>;

  std::cout << "--- Warm up ";
  std::mt19937 rng(std::random_device{}());
  rng.discard(700000);
  std::cout << std::uniform_int_distribution<int>{0, 2147483647}(rng) << "\n";

  std::string jsonText = "{";

  std::cout << std::format(ROW_STR, "Name", "Elapsed [ms]", "Last Output");
  constexpr Sample sampleRateHz = Sample(48000);
  constexpr Sample startFreq = Sample(10) / sampleRateHz;
  constexpr Sample endFreq = Sample(10000) / sampleRateHz;
  jsonText += benchmark<Sample, SineFloor<Sample>>(startFreq, endFreq);
  jsonText += benchmark<Sample, SineFmod<Sample>>(twopi * startFreq, twopi * endFreq);

  jsonText.pop_back();
  jsonText += "}";

  std::ofstream os("sine.json");
  os << jsonText;

  return 0;
}
