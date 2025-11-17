#include <algorithm>
#include <array>
#include <chrono>
#include <cinttypes>
#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <random>
#include <string>
#include <vector>

/*
`LowpassFir` is a linear-phase 2-band splitter. Cutoff can be chagned for each sample.

The compulation of windowed-sinc uses recursive sine algorithm (biquad oscillator) instead
of `std::sin`, to make it runs faster.
*/
template<typename Sample, int length = 63> class LowpassFir {
private:
  static_assert(length > 0);

  int wptr = 0;
  std::array<Sample, length> buf{};

public:
  static constexpr int latency = length / 2;

  std::string name() { return "ptr_rectangular"; }

  void reset()
  {
    wptr = 0;
    buf.fill({});
  }

  std::array<Sample, 2> process(Sample input, Sample cutoffNormalized)
  {
    if constexpr (length == 1) return Sample(2) * cutoffNormalized * input;

    // Write to buffer.
    if (++wptr >= length) wptr = 0;
    buf[wptr] = input;

    // Setup recursive sine oscillator.
    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample omega = Sample(2) * pi * cutoffNormalized;
    const Sample phi = -latency * omega;
    const Sample k = Sample(2) * std::cos(omega);
    Sample u1 = std::sin(phi - omega);
    Sample u2 = std::sin(phi - Sample(2) * omega);

    // Convolution.
    int rptr = wptr;
    Sample sum = 0;
    for (int i = 0; i < length; ++i) {
      const Sample u0 = k * u1 - u2;
      u2 = u1;
      u1 = u0;

      const Sample x = Sample(i) - latency;
      const Sample sinc = x == 0 ? Sample(2) * cutoffNormalized : u0 / (pi * x);

      if (++rptr >= length) rptr = 0;
      sum += sinc * buf[rptr];
    }

    int center = wptr - latency;
    if (center < 0) center += length;
    return {sum, buf[center] - sum};
  }
};

/*
`LowpassFirRotate` tries to save the computation of mirror symmetric part of the FIR
filter. For example, when the filter coeffiencts is [1, 2, 3, 2, 1], it only compute first
half [1, 2], and convolute the values to the buffer.

Slower than `LowpassFir`.
*/
template<typename Sample, int length = 63> class LowpassFirBidirection {
private:
  static_assert(length > 0);

  int wptr = 0;
  std::array<Sample, length> buf{};

public:
  static constexpr int latency = length / 2;

  std::string name() { return "bidi_rectangular"; }

  void reset()
  {
    wptr = 0;
    buf.fill({});
  }

  std::array<Sample, 2> process(Sample input, Sample cutoffNormalized)
  {
    if constexpr (length == 1) return Sample(2) * cutoffNormalized * input;

    // Write to buffer.
    if (++wptr >= length) wptr = 0;
    buf[wptr] = input;

    // Setup recursive sine oscillator.
    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample omega = Sample(2) * pi * cutoffNormalized;
    const Sample phi = -latency * omega;
    const Sample k = Sample(2) * std::cos(omega);
    Sample u1 = std::sin(phi - omega);
    Sample u2 = std::sin(phi - Sample(2) * omega);

    // Convolution.
    int rptr0 = wptr;
    int rptr1 = wptr + 1;
    Sample sum = 0;
    for (int i = 0; i < latency; ++i) {
      const Sample u0 = k * u1 - u2;
      u2 = u1;
      u1 = u0;

      const Sample sinc = u0 / (pi * Sample(i - latency));

      if (++rptr0 >= length) rptr0 = 0;
      if (--rptr1 < 0) rptr1 = length - 1;
      sum += sinc * (buf[rptr0] + buf[rptr1]);
    }

    int center = wptr - latency;
    if (center < 0) center += length;
    sum += Sample(2) * cutoffNormalized * buf[center];

    return {sum, buf[center] - sum};
  }
};

/*
`LowpassFirRotate` uses `std::rotate` instead of juggling indices with `wptr` and `rptr`.

Slower than `LowpassFir`.
*/
template<typename Sample, int length = 63> class LowpassFirRotate {
private:
  static_assert(length > 0);
  std::array<Sample, length> buf{};

public:
  static constexpr int latency = length / 2;

  std::string name() { return "rotate_rectangular"; }

  void reset() { buf.fill({}); }

  std::array<Sample, 2> process(Sample input, Sample cutoffNormalized)
  {
    if constexpr (length == 1) return Sample(2) * cutoffNormalized * input;

    std::rotate(buf.begin(), buf.begin() + 1, buf.end());
    buf[length - 1] = input;

    // Setup recursive sine oscillator.
    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample omega = Sample(2) * pi * cutoffNormalized;
    const Sample phi = -latency * omega;
    const Sample k = Sample(2) * std::cos(omega);
    Sample u1 = std::sin(phi - omega);
    Sample u2 = std::sin(phi - Sample(2) * omega);

    // Convolution.
    Sample sum = 0;
    for (int i = 0; i < length; ++i) {
      const Sample u0 = k * u1 - u2;
      u2 = u1;
      u1 = u0;

      const Sample x = Sample(i - latency);
      const Sample sinc = x == 0 ? Sample(2) * cutoffNormalized : u0 / (pi * x);

      sum += sinc * buf[i];
    }

    return {sum, buf[latency] - sum};
  }
};

template<typename Float>
std::string vectorToJsonList(const std::string &name, const std::vector<Float> &data)
{
  std::string text = std::format("\"{}\":[", name);
  for (const auto &value : data) text += std::format("{},", value);
  text.pop_back();
  text += "],";
  return text;
}

template<typename Filter, typename Sample> std::string testImpulseResponse()
{
  constexpr Sample cutoffNormalzied = Sample(0.25);
  Filter filter;

  std::string text;
  text += std::format("\"cutoff\":{},", cutoffNormalzied);
  text += std::format("\"latency\":{},", filter.latency);

  std::vector<Sample> impulse;
  impulse.resize(2048);
  impulse[0] = Sample(1);

  std::vector<Sample> low(impulse.size(), Sample(0));
  std::vector<Sample> high(impulse.size(), Sample(0));

  for (size_t i = 0; i < impulse.size(); ++i) {
    auto frame = filter.process(impulse[i], cutoffNormalzied);
    low[i] = frame[0];
    high[i] = frame[1];
  }
  text += vectorToJsonList("testIrLow", low);
  text += vectorToJsonList("testIrHigh", high);

  filter.reset();
  for (size_t i = 0; i < impulse.size(); ++i) {
    auto frame = filter.process(impulse[i], cutoffNormalzied);
    low[i] = frame[0];
    high[i] = frame[1];
  }
  text += vectorToJsonList("testResetLow", low);
  text += vectorToJsonList("testResetHigh", high);

  return text;
}

#define ROW_FORMAT_STR "|{:20}|{:20}|{:20}|{:20}|\n"

template<typename Filter, typename Sample> void benchmark()
{
  constexpr Sample sampleRate = Sample(48000);
  constexpr size_t nSample = 20 * size_t(sampleRate);

  Filter filter;
  auto rng = std::minstd_rand(3216748);
  std::uniform_real_distribution<Sample> dist(Sample(-0.25), Sample(0.25));
  double sumElapsed = 0.0;
  for (size_t i = 0; i < nSample; ++i) {
    const auto value = dist(rng);

    Sample cutoff = std::exp2(
      Sample(-3)
      + std::sin(Sample(4.01) * std::numbers::pi_v<Sample> * Sample(i) / sampleRate));

    auto start = std::chrono::steady_clock::now();
    filter.process(value, cutoff);
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }

  const auto frame = filter.process(0, 0);
  std::cout << std::format(
    ROW_FORMAT_STR, filter.name(), sumElapsed, nSample, frame[0] + frame[1]);
}

void testAll()
{
  std::string text{"{"};

  text += testImpulseResponse<LowpassFir<double, 15>, double>();
  // text += testImpulseResponse<LowpassFirBidirection<double, 15>, double>();
  // text += testImpulseResponse<LowpassFirRotate<double, 15>, double>();

  text.pop_back();
  text += "}";
  std::ofstream ofs("cpp_test_output.json");
  ofs << text;
}

void benchmarkAll()
{
  std::cout << "--- Warm up ";
  std::mt19937 rng(std::random_device{}());
  rng.discard(700000);
  std::cout << std::uniform_int_distribution<int>{0, 2147483647}(rng) << "\n";

  std::cout << std::format(
    ROW_FORMAT_STR, "Name", "Time [ms]", "Sample Count", "Last Output")
            << std::format(ROW_FORMAT_STR, "--", "--", "--", "--");

  benchmark<LowpassFir<double>, double>();
  benchmark<LowpassFirBidirection<double>, double>();
  benchmark<LowpassFirRotate<double>, double>();
}

int main()
{
  testAll();
  benchmarkAll();
  return 0;
}
