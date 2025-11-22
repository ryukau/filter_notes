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

template<typename Sample> inline Sample modifiedSinc(Sample x, Sample cutoff)
{
  constexpr Sample pi = std::numbers::pi_v<Sample>;
  return std::abs(x) < std::numeric_limits<Sample>::epsilon()
    ? Sample(2) * cutoff
    : std::sin(Sample(2) * pi * cutoff * x) / (pi * x);
}

// `length` may be shorter than `size`. `0 <= length < size`.
template<typename Sample, size_t size>
inline void
fillFirEven(std::array<Sample, size> &fir, int length, Sample cutoff, Sample fraction)
{
  const Sample mid = Sample(length / 2) - fraction;
  for (int i = 0; i < length; ++i) {
    const Sample x = Sample(i) - mid;
    fir[i] = modifiedSinc(x, cutoff);
  }
}

// No interplation.
// `maxTap` is unused. It is only there to match the interface.
template<typename Sample, int maxTap = 256> class DelayInt {
private:
  int wptr = 0;
  std::vector<Sample> buf{Sample(0), Sample(0)};

public:
  std::string name() { return "int"; }

  void setup(Sample maxTimeSamples)
  {
    buf.resize(std::max(size_t(1), size_t(maxTimeSamples) + 1));
    reset();
  }

  void reset() { std::fill(buf.begin(), buf.end(), Sample(0)); }

  Sample process(Sample input, Sample timeInSamples)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const Sample clamped = std::clamp(timeInSamples, Sample(0), Sample(size - 1));
    auto rptr0 = wptr - int(clamped);
    if (rptr0 < 0) rptr0 += size;
    return buf[rptr0];
  }
};

// Linear interplation.
// `maxTap` is unused. It is only there to match the interface.
template<typename Sample, int maxTap = 256> class DelayLinear {
private:
  int wptr = 0;
  std::vector<Sample> buf{Sample(0), Sample(0)};

public:
  std::string name() { return "linear"; }

  void setup(Sample maxTimeSamples)
  {
    buf.resize(std::max(size_t(2), size_t(maxTimeSamples) + 2));
    reset();
  }

  void reset() { std::fill(buf.begin(), buf.end(), Sample(0)); }

  Sample process(Sample input, Sample timeInSamples)
  {
    const int size = int(buf.size());
    const Sample clamped = std::clamp(timeInSamples, Sample(0), Sample(size - 1));
    const int timeInt = int(clamped);
    const Sample rFraction = clamped - Sample(timeInt);

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    auto rptr0 = wptr - timeInt;
    auto rptr1 = rptr0 - 1;
    if (rptr0 < 0) rptr0 += size;
    if (rptr1 < 0) rptr1 += size;
    return std::lerp(buf[rptr0], buf[rptr1], rFraction);
  }
};

// Order 3 Lagrange interplation.
// `maxTap` is unused. It is only there to match the interface.
//
// `process1` can't set delay time below 1 sample. It may sounds smoother than `process`.
template<typename Sample, int maxTap = 256> class DelayLagrange3 {
private:
  int wptr = 0;
  std::vector<Sample> buf{Sample(0), Sample(0)};

  inline Sample lagrange3Interp(Sample y0, Sample y1, Sample y2, Sample y3, Sample t)
  {
    auto u = Sample(1) + t;
    auto d0 = y0 - y1;
    auto d1 = d0 - (y1 - y2);
    auto d2 = d1 - ((y1 - y2) - (y2 - y3));
    return y0
      - u * (d0 + (Sample(1) - u) / Sample(2) * (d1 + (Sample(2) - u) / Sample(3) * d2));
  }

public:
  std::string name() { return "lagrange3"; }

  void setup(Sample maxTimeSamples)
  {
    buf.resize(std::max(size_t(4), size_t(maxTimeSamples) + 4));
    reset();
  }

  void reset() { std::fill(buf.begin(), buf.end(), Sample(0)); }

  Sample process(Sample input, Sample timeInSamples)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const Sample clamped
      = std::clamp(timeInSamples - Sample(1), Sample(-1), Sample(size - 4));
    const Sample rFraction = clamped - std::floor(clamped);

    if (clamped < 0) { // Linear interpolation when delay time is less than 1 sample.
      auto rptr1 = wptr - 1;
      if (rptr1 < 0) rptr1 += size;
      return std::lerp(buf[wptr], buf[rptr1], rFraction);
    }

    auto rptr0 = wptr - int(clamped);
    auto rptr1 = rptr0 - 1;
    auto rptr2 = rptr0 - 2;
    auto rptr3 = rptr0 - 3;
    if (rptr0 < 0) rptr0 += size;
    if (rptr1 < 0) rptr1 += size;
    if (rptr2 < 0) rptr2 += size;
    if (rptr3 < 0) rptr3 += size;
    return lagrange3Interp(buf[rptr0], buf[rptr1], buf[rptr2], buf[rptr3], rFraction);
  }

  Sample process1(Sample input, Sample timeInSamples)
  {
    const int size = int(buf.size());
    const Sample clamped
      = std::clamp(timeInSamples - Sample(1), Sample(0), Sample(size - 4));
    const int timeInt = int(clamped);
    const Sample rFraction = clamped - Sample(timeInt);

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    auto rptr0 = wptr - timeInt;
    auto rptr1 = rptr0 - 1;
    auto rptr2 = rptr0 - 2;
    auto rptr3 = rptr0 - 3;
    if (rptr0 < 0) rptr0 += size;
    if (rptr1 < 0) rptr1 += size;
    if (rptr2 < 0) rptr2 += size;
    if (rptr3 < 0) rptr3 += size;
    return lagrange3Interp(buf[rptr0], buf[rptr1], buf[rptr2], buf[rptr3], rFraction);
  }
};

// Reference implementation. Slow but simple and accurate. Rectangular window.
template<typename Sample, int maxTap = 256> class DelayAntialiasedReference {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf;
  std::array<Sample, maxTap> fir{};

public:
  std::string name() { return "naive"; }

  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
    fir.fill(Sample(0));
  }

  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);
    if (timeInSample <= 0) return input * modifiedSinc(Sample(0), cutoff);

    const int timeInt = int(clamped);
    const Sample fraction = clamped - timeInt;
    fillFirEven(fir, localTap, cutoff, fraction);

    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;

    Sample sig = 0;
    for (int idx = 0; idx < localTap; ++idx) {
      sig += fir[idx] * buf[rptr];
      if (++rptr >= size) rptr = 0;
    }
    return sig;
  }
};

// An attempt to speed up the reference implementation.
template<typename Sample, int maxTap = 256> class DelayAntialiasedInnerProduct {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf;
  std::array<Sample, maxTap> fir{};
  std::array<Sample, maxTap> read{};

public:
  std::string name() { return "inner_product"; }

  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
    fir.fill(Sample(0));
    read.fill(Sample(0));
  }

  // Using `std::inner_product`.
  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const int localTap = std::clamp(2 * int(timeInSample), 2, maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);
    if (timeInSample <= 0) return input * Sample(2) * cutoff;

    const int timeInt = int(clamped);
    const Sample fraction = clamped - timeInt;
    fillFirEven(fir, localTap, cutoff, fraction);

    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;
    for (int idx = 0; idx < localTap; ++idx) {
      read[idx] = buf[rptr];
      if (++rptr >= size) rptr = 0;
    }

    return std::inner_product(read.begin(), read.end(), fir.begin(), Sample(0));
  }
};

// Another attempt to speed up the reference implementation.
template<typename Sample, int maxTap = 256> class DelayAntialiasedInlined {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf;

public:
  std::string name() { return "inlined"; }

  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
  }

  // FIR coefficient computations are inlined into a loop.
  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);
    if (timeInSample <= 0) return input * Sample(2) * cutoff;

    const int timeInt = int(clamped);
    const Sample mid = halfTap + timeInt - clamped;

    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;

    Sample sum = 0;
    for (int i = 0; i < localTap; ++i) {
      const Sample sinc = modifiedSinc(Sample(i) - mid, cutoff);
      sum += sinc * buf[rptr];
      if (++rptr >= size) rptr = 0;
    }
    return sum;
  }
};

// Fast implementation with severe inaccuracy when `x` is around 0.
template<typename Sample, int maxTap = 256> class DelayAntialiasedBiquadNaive {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf;

public:
  std::string name() { return "biquad_sine"; }

  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
  }

  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);

    if (timeInSample <= 0) return input * Sample(2) * cutoff;

    // Setup recursive sine oscillator.
    const int timeInt = int(clamped);
    Sample fraction = clamped - Sample(timeInt);
    const Sample mid = fraction - halfTap;

    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample omega = Sample(2) * pi * cutoff;
    const Sample phi = mid * omega;
    const Sample k = Sample(2) * std::cos(omega);
    Sample u1 = std::sin(phi - omega);
    Sample u2 = std::sin(phi - Sample(2) * omega);

    // Convolution.
    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;

    Sample sum = 0;
    for (int i = 0; i < localTap; ++i) {
      const Sample u0 = k * u1 - u2;
      u2 = u1;
      u1 = u0;

      const Sample x = Sample(i) + mid;
      const Sample sinc = x == 0 ? Sample(2) * cutoff : u0 / (pi * x);
      sum += sinc * buf[rptr];
      if (++rptr >= size) rptr = 0;
    }
    return sum;
  }
};

// Fast implementation with usable accuracy.
template<typename Sample, int maxTap = 256> class DelayAntialiasedBiquad2 {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf;

public:
  std::string name() { return "biquad_less_error"; }

  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
  }

  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);

    if (timeInSample <= 0) return input * Sample(2) * cutoff;

    // Setup recursive sine oscillator.
    const int timeInt = int(clamped);
    Sample fraction = clamped - Sample(timeInt);
    const Sample mid = fraction - halfTap;

    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample omega = Sample(2) * pi * cutoff;
    const Sample phi = mid * omega;
    const Sample cs_o = std::cos(omega);
    const Sample sn_o = std::sin(omega);
    const Sample cs_p = std::cos(phi);
    const Sample sn_p = std::sin(phi);
    const Sample k = Sample(2) * cs_o;
    Sample u1 = sn_p * cs_o - cs_p * sn_o;
    const Sample cs_o2 = Sample(2) * cs_o * cs_o - Sample(1);
    const Sample sn_o2 = Sample(2) * sn_o * cs_o;
    Sample u2 = sn_p * cs_o2 - cs_p * sn_o2;

    // Convolution.
    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;

    Sample sum = 0;
    for (int i = 0; i < localTap; ++i) {
      const Sample u0 = k * u1 - u2;
      u2 = u1;
      u1 = u0;

      const Sample x = Sample(i) + mid;
      Sample sinc;
      if (std::abs(x) <= Sample(0.1)) {
        Sample q = pi * cutoff * x;
        q *= q;
        sinc = Sample(2) / Sample(3) * cutoff * (Sample(15) - Sample(7) * q)
          / (Sample(5) + q);
      } else {
        sinc = u0 / (pi * x);
      }
      sum += sinc * buf[rptr];
      if (++rptr >= size) rptr = 0;
    }
    return sum;
  }
};

// Fast implementation with Reinsch oscillator.
template<typename Sample, int maxTap = 256> class DelayAntialiasedReinsch {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf;

public:
  std::string name() { return "reinsch_sine"; }

  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
  }

  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);

    if (timeInSample <= 0) return input * Sample(2) * cutoff;

    // Setup recursive sine oscillator.
    const int timeInt = int(clamped);
    Sample fraction = clamped - Sample(timeInt);
    const Sample mid = fraction - halfTap;

    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample omega = Sample(2) * pi * cutoff;
    const Sample phi = mid * omega;
    const Sample A = Sample(2) * std::sin(omega / Sample(2));
    const Sample k = A * A;
    Sample u = std::sin(phi - omega);
    Sample v = A * std::cos(phi - omega / Sample(2));

    // Convolution.
    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;

    Sample sum = 0;
    for (int i = 0; i < localTap; ++i) {
      u = u + v;
      v = v - k * u;

      const Sample x = Sample(i) + mid;
      Sample sinc;
      if (std::abs(x) <= Sample(0.1)) {
        Sample q = pi * cutoff * x;
        q *= q;
        sinc = Sample(2) / Sample(3) * cutoff * (Sample(15) - Sample(7) * q)
          / (Sample(5) + q);
      } else {
        sinc = u / (pi * x);
      }

      sum += sinc * buf[rptr];
      if (++rptr >= size) rptr = 0;
    }
    return sum;
  }
};

// Fast implementation with magic circle oscillator.
template<typename Sample, int maxTap = 256> class DelayAntialiasedMagic {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf;

public:
  std::string name() { return "magic_circle"; }

  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
  }

  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);

    if (timeInSample <= 0) return input * Sample(2) * cutoff;

    // Setup recursive sine oscillator.
    const int timeInt = int(clamped);
    Sample fraction = clamped - Sample(timeInt);
    const Sample mid = fraction - halfTap;

    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample omega = Sample(2) * pi * cutoff;
    const Sample phi = mid * omega;
    const Sample k = Sample(2) * std::sin(omega / Sample(2));
    Sample u = std::cos(phi - omega * Sample(3) / Sample(2));
    Sample v = std::sin(phi - omega);

    // Convolution.
    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;

    Sample sum = 0;
    for (int i = 0; i < localTap; ++i) {
      u -= k * v;
      v += k * u;

      const Sample x = Sample(i) + mid;
      Sample sinc;
      if (std::abs(x) <= Sample(0.1)) {
        Sample q = pi * cutoff * x;
        q *= q;
        sinc = Sample(2) / Sample(3) * cutoff * (Sample(15) - Sample(7) * q)
          / (Sample(5) + q);
      } else {
        sinc = v / (pi * x);
      }

      sum += sinc * buf[rptr];
      if (++rptr >= size) rptr = 0;
    }
    return sum;
  }
};

// Fast implementation with stable quad oscillator.
//
// For this application, coupled form (2D rotation matrix) might be better because of the
// simpler initial phase computation that doesn't involve `tan`.
template<typename Sample, int maxTap = 256> class DelayAntialiasedStableQuad {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf;

public:
  std::string name() { return "stable_quad_sine"; }

  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
  }

  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);

    if (timeInSample <= 0) return input * Sample(2) * cutoff;

    // Setup recursive sine oscillator.
    const int timeInt = int(clamped);
    Sample fraction = clamped - Sample(timeInt);
    const Sample mid = fraction - halfTap;

    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample omega = Sample(2) * pi * cutoff;
    const Sample phi = mid * omega;
    const Sample k1 = std::tan(omega / Sample(2));
    const Sample k2 = std::sin(omega);
    Sample u = std::cos(phi - omega);
    Sample v = std::sin(phi - omega);

    // Convolution.
    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;

    Sample sum = 0;
    for (int i = 0; i < localTap; ++i) {
      const Sample w = u - k1 * v;
      v += k2 * w;
      u = w - k1 * v;

      const Sample x = Sample(i) + mid;
      Sample sinc;
      if (std::abs(x) <= Sample(0.1)) {
        Sample q = pi * cutoff * x;
        q *= q;
        sinc = Sample(2) / Sample(3) * cutoff * (Sample(15) - Sample(7) * q)
          / (Sample(5) + q);
      } else {
        sinc = v / pi / x;
      }

      sum += sinc * buf[rptr];
      if (++rptr >= size) rptr = 0;
    }
    return sum;
  }
};

// Fast implementation with Lanczos window that adapts to cutoff frequency.
// Low quality, too much noise is added.
template<typename Sample, int maxTap = 256> class DelayLanczosABiquadSine {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf;

public:
  std::string name() { return "fast_lanczosA"; }

  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
  }

  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);

    if (timeInSample <= 0) return input * Sample(2) * cutoff;

    // Setup convolution filter parameters.
    const int timeInt = int(clamped);
    Sample fraction = clamped - Sample(timeInt);
    const Sample mid = fraction - halfTap;

    // Setup oscillator 1 (o1). Windowed sinc lowpass.
    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample o1_omega = Sample(2) * pi * cutoff;
    const Sample o1_phi = mid * o1_omega;
    const Sample o1_k = Sample(2) * std::cos(o1_omega);
    Sample o1_u1 = std::sin(o1_phi - o1_omega);
    Sample o1_u2 = std::sin(o1_phi - Sample(2) * o1_omega);

    // Setup oscillator 2 (o2). Lanczos window. `a` is sidelobe parameter.
    // This Lanczos window adapts to cutoff which is not following usual definition.
    const Sample lanczos_a = std::max(Sample(1), std::sqrt(halfTap));
    const Sample o2_omega = Sample(2) * pi * cutoff / lanczos_a;
    const Sample o2_phi = mid * o2_omega;
    const Sample o2_k = Sample(2) * std::cos(o2_omega);
    Sample o2_u1 = std::sin(o2_phi - o2_omega);
    Sample o2_u2 = std::sin(o2_phi - Sample(2) * o2_omega);

    // Convolution.
    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;

    Sample sum = 0;
    const Sample A = lanczos_a / (Sample(2) * cutoff * pi * pi);
    for (int i = 0; i < localTap; ++i) {
      const Sample o1_u0 = o1_k * o1_u1 - o1_u2;
      o1_u2 = o1_u1;
      o1_u1 = o1_u0;

      const Sample o2_u0 = o2_k * o2_u1 - o2_u2;
      o2_u2 = o2_u1;
      o2_u1 = o2_u0;

      const Sample x = Sample(i) + mid;
      Sample sinc;
      if (std::abs(x) <= Sample(0.1)) {
        Sample q = pi * x * cutoff;
        q *= q;
        sinc = Sample(2) / Sample(3) * cutoff * (Sample(15) - Sample(7) * q)
          / (Sample(5) + q);

        Sample p = pi * x * cutoff / lanczos_a;
        p *= p;
        sinc *= Sample(2) / Sample(3) * cutoff / lanczos_a * (Sample(15) - Sample(7) * p)
          / (Sample(5) + p);
      } else {
        sinc = o1_u0 * o2_u0 * A / x / x;
      }

      sum += sinc * buf[rptr];
      if (++rptr >= size) rptr = 0;
    }
    return sum;
  }
};

// Fast implementation with Lanczos window with sidelobe parameter `a` fixed to 1.
template<typename Sample, int maxTap = 256> class DelayLanczos1BiquadSine {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf;

public:
  std::string name() { return "fast_lanczos1"; }

  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
  }

  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);

    if (timeInSample <= 0) return input * Sample(2) * cutoff;

    // Setup convolution filter parameters.
    const int timeInt = int(clamped);
    Sample fraction = clamped - Sample(timeInt);
    const Sample mid = fraction - halfTap;

    // Setup oscillator 1 (o1). Windowed sinc lowpass.
    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample o1_omega = Sample(2) * pi * cutoff;
    const Sample o1_phi = mid * o1_omega;
    const Sample o1_k = Sample(2) * std::cos(o1_omega);
    Sample o1_u1 = std::sin(o1_phi - o1_omega);
    Sample o1_u2 = std::sin(o1_phi - Sample(2) * o1_omega);

    // Setup oscillator 2 (o2). Lanczos window.
    constexpr Sample lanczos_a = Sample(1);
    const Sample o2_omega = Sample(2) * pi / lanczos_a;
    const Sample o2_phi = mid * o2_omega;
    const Sample o2_k = Sample(2) * std::cos(o2_omega);
    Sample o2_u1 = std::sin(o2_phi - o2_omega);
    Sample o2_u2 = std::sin(o2_phi - Sample(2) * o2_omega);

    // Convolution.
    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;

    Sample sum = 0;
    const Sample A = lanczos_a / (Sample(2) * pi * pi);
    for (int i = 0; i < localTap; ++i) {
      const Sample o1_u0 = o1_k * o1_u1 - o1_u2;
      o1_u2 = o1_u1;
      o1_u1 = o1_u0;

      const Sample o2_u0 = o2_k * o2_u1 - o2_u2;
      o2_u2 = o2_u1;
      o2_u1 = o2_u0;

      const Sample x = Sample(i) + mid;
      Sample sinc;
      if (std::abs(x) <= Sample(0.1)) {
        Sample q = pi * x * cutoff;
        q *= q;
        sinc = Sample(2) / Sample(3) * cutoff * (Sample(15) - Sample(7) * q)
          / (Sample(5) + q);

        Sample p = pi * x / lanczos_a;
        p *= p;
        sinc *= Sample(2) / Sample(3) * cutoff / lanczos_a * (Sample(15) - Sample(7) * p)
          / (Sample(5) + p);
      } else {
        sinc = o1_u0 * o2_u0 * A / x / x;
      }

      sum += sinc * buf[rptr];
      if (++rptr >= size) rptr = 0;
    }
    return sum;
  }
};

#define DELAY_COSINE_SUM(WINDOW_NAME, WINDOW_POLYNOMIAL)                                 \
  template<typename Sample, int maxTap = 256> class Delay##WINDOW_NAME##BiquadSine {     \
  private:                                                                               \
    static_assert(maxTap > 0 && maxTap % 2 == 0);                                        \
                                                                                         \
    static constexpr int minTimeSample = maxTap / 2 - 1;                                 \
                                                                                         \
    Sample maxTime = 0;                                                                  \
    Sample prevTime = 0;                                                                 \
    int wptr = 0;                                                                        \
    std::vector<Sample> buf;                                                             \
                                                                                         \
  public:                                                                                \
    std::string name() { return "fast_" #WINDOW_NAME; }                                  \
                                                                                         \
    void setup(size_t maxTimeSample)                                                     \
    {                                                                                    \
      maxTime = Sample(maxTimeSample);                                                   \
      buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));              \
    }                                                                                    \
                                                                                         \
    void reset()                                                                         \
    {                                                                                    \
      prevTime = 0;                                                                      \
      wptr = 0;                                                                          \
      std::fill(buf.begin(), buf.end(), Sample(0));                                      \
    }                                                                                    \
                                                                                         \
    Sample process(Sample input, Sample timeInSample)                                    \
    {                                                                                    \
      const int size = int(buf.size());                                                  \
                                                                                         \
      if (++wptr >= size) wptr = 0;                                                      \
      buf[wptr] = input;                                                                 \
                                                                                         \
      const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);            \
      const int halfTap = localTap / 2;                                                  \
      const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);     \
                                                                                         \
      const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));                  \
      prevTime = clamped;                                                                \
      Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);        \
                                                                                         \
      if (timeInSample <= 0) return input * Sample(2) * cutoff;                          \
                                                                                         \
      const int timeInt = int(clamped);                                                  \
      Sample fraction = clamped - Sample(timeInt);                                       \
      const Sample mid = fraction - halfTap;                                             \
                                                                                         \
      constexpr Sample pi = std::numbers::pi_v<Sample>;                                  \
      const Sample o1_omega = Sample(2) * pi * cutoff;                                   \
      const Sample o1_phi = mid * o1_omega;                                              \
      const Sample o1_k = Sample(2) * std::cos(o1_omega);                                \
      Sample o1_u1 = std::sin(o1_phi - o1_omega);                                        \
      Sample o1_u2 = std::sin(o1_phi - Sample(2) * o1_omega);                            \
                                                                                         \
      const Sample o2_omega = Sample(2) * pi / Sample(localTap + 1);                     \
      const Sample o2_phi = pi / Sample(2);                                              \
      const Sample o2_k = Sample(2) * std::cos(o2_omega);                                \
      Sample o2_u1 = Sample(1);                                                          \
      Sample o2_u2 = std::sin(o2_phi - o2_omega);                                        \
                                                                                         \
      int rptr = wptr - timeInt - halfTap;                                               \
      if (rptr < 0) rptr += size;                                                        \
                                                                                         \
      Sample sum = 0;                                                                    \
      for (int i = 0; i < localTap; ++i) {                                               \
        const Sample o1_u0 = o1_k * o1_u1 - o1_u2;                                       \
        o1_u2 = o1_u1;                                                                   \
        o1_u1 = o1_u0;                                                                   \
                                                                                         \
        const Sample o2_u0 = o2_k * o2_u1 - o2_u2;                                       \
        o2_u2 = o2_u1;                                                                   \
        o2_u1 = o2_u0;                                                                   \
                                                                                         \
        const Sample window = WINDOW_POLYNOMIAL;                                         \
                                                                                         \
        const Sample x = Sample(i) + mid;                                                \
        Sample sinc;                                                                     \
        if (std::abs(x) <= Sample(0.1)) {                                                \
          Sample q = pi * cutoff * x;                                                    \
          q *= q;                                                                        \
          sinc = Sample(2) / Sample(3) * cutoff * (Sample(15) - Sample(7) * q)           \
            / (Sample(5) + q);                                                           \
        } else {                                                                         \
          sinc = o1_u0 / (pi * x);                                                       \
        }                                                                                \
                                                                                         \
        sum += window * sinc * buf[rptr];                                                \
        if (++rptr >= size) rptr = 0;                                                    \
      }                                                                                  \
      return sum;                                                                        \
    }                                                                                    \
  };

#define DELAY_COSINE_SUM_SMOOTH(WINDOW_NAME, WINDOW_POLYNOMIAL)                          \
  template<typename Sample, int maxTap = 256>                                            \
  class Delay##WINDOW_NAME##SmoothBiquadSine {                                           \
  private:                                                                               \
    static_assert(maxTap > 0 && maxTap % 2 == 0);                                        \
                                                                                         \
    static constexpr int minTimeSample = maxTap / 2 - 1;                                 \
                                                                                         \
    Sample maxTime = 0;                                                                  \
    Sample prevTime = 0;                                                                 \
    int wptr = 0;                                                                        \
    std::vector<Sample> buf;                                                             \
                                                                                         \
  public:                                                                                \
    std::string name() { return "smoo_" #WINDOW_NAME; }                                  \
                                                                                         \
    void setup(size_t maxTimeSample)                                                     \
    {                                                                                    \
      maxTime = Sample(maxTimeSample);                                                   \
      buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));              \
    }                                                                                    \
                                                                                         \
    void reset()                                                                         \
    {                                                                                    \
      prevTime = 0;                                                                      \
      wptr = 0;                                                                          \
      std::fill(buf.begin(), buf.end(), Sample(0));                                      \
    }                                                                                    \
                                                                                         \
    Sample process(Sample input, Sample timeInSample)                                    \
    {                                                                                    \
      const int size = int(buf.size());                                                  \
                                                                                         \
      if (++wptr >= size) wptr = 0;                                                      \
      buf[wptr] = input;                                                                 \
                                                                                         \
      const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);            \
      const int halfTap = localTap / 2;                                                  \
      const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);     \
                                                                                         \
      const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));                  \
      prevTime = clamped;                                                                \
      Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);        \
                                                                                         \
      if (timeInSample <= 0) return input * Sample(2) * cutoff;                          \
                                                                                         \
      const int timeInt = int(clamped);                                                  \
      Sample fraction = clamped - Sample(timeInt);                                       \
      const Sample mid = fraction - halfTap;                                             \
                                                                                         \
      constexpr Sample pi = std::numbers::pi_v<Sample>;                                  \
      const Sample o1_omega = Sample(2) * pi * cutoff;                                   \
      const Sample o1_phi = mid * o1_omega;                                              \
      const Sample o1_k = Sample(2) * std::cos(o1_omega);                                \
      Sample o1_u1 = std::sin(o1_phi - o1_omega);                                        \
      Sample o1_u2 = std::sin(o1_phi - Sample(2) * o1_omega);                            \
                                                                                         \
      const Sample o2_omega = Sample(2) * pi / Sample(maxTap + 1);                       \
      const Sample o2_phi = pi / Sample(2) + o2_omega * Sample(maxTap / 2 - halfTap);    \
      const Sample o2_k = Sample(2) * std::cos(o2_omega);                                \
      Sample o2_u1 = std::sin(o2_phi);                                                   \
      Sample o2_u2 = std::sin(o2_phi - o2_omega);                                        \
                                                                                         \
      int rptr = wptr - timeInt - halfTap;                                               \
      if (rptr < 0) rptr += size;                                                        \
                                                                                         \
      Sample sum = 0;                                                                    \
      for (int i = 0; i < localTap; ++i) {                                               \
        const Sample o1_u0 = o1_k * o1_u1 - o1_u2;                                       \
        o1_u2 = o1_u1;                                                                   \
        o1_u1 = o1_u0;                                                                   \
                                                                                         \
        const Sample o2_u0 = o2_k * o2_u1 - o2_u2;                                       \
        o2_u2 = o2_u1;                                                                   \
        o2_u1 = o2_u0;                                                                   \
                                                                                         \
        const Sample window = WINDOW_POLYNOMIAL;                                         \
                                                                                         \
        const Sample x = Sample(i) + mid;                                                \
        Sample sinc;                                                                     \
        if (std::abs(x) <= Sample(0.1)) {                                                \
          Sample q = pi * cutoff * x;                                                    \
          q *= q;                                                                        \
          sinc = Sample(2) / Sample(3) * cutoff * (Sample(15) - Sample(7) * q)           \
            / (Sample(5) + q);                                                           \
        } else {                                                                         \
          sinc = o1_u0 / (pi * x);                                                       \
        }                                                                                \
                                                                                         \
        sum += window * sinc * buf[rptr];                                                \
        if (++rptr >= size) rptr = 0;                                                    \
      }                                                                                  \
      return sum;                                                                        \
    }                                                                                    \
  };

#define BLACKMAN_POLY                                                                    \
  Sample(0.349742046431642)                                                              \
    + o2_u0 *(Sample(-0.496560619088564) + Sample(0.153697334479794) * o2_u0)

#define NUTTALL_POLY                                                                     \
  Sample(0.211536)                                                                       \
    + o2_u0 *(                                                                           \
      Sample(-0.449584) + o2_u0 * (Sample(0.288464) + o2_u0 * Sample(-0.050416)))

#define BLACKMANHARRIS_POLY                                                              \
  Sample(0.21747)                                                                        \
    + o2_u0 *(Sample(-0.45325) + o2_u0 * (Sample(0.28256) + o2_u0 * Sample(-0.04672)))

#define BLACKMANNUTTALL_POLY                                                             \
  Sample(0.2269824)                                                                      \
    + o2_u0 *(                                                                           \
      Sample(-0.4572542) + o2_u0 * (Sample(0.273199) + o2_u0 * Sample(-0.0425644)))

#define FLATTOP_POLY                                                                     \
  Sample(-0.05473684)                                                                    \
    + o2_u0 *(                                                                           \
      Sample(-0.165894739)                                                               \
      + o2_u0                                                                            \
        * (Sample(0.498947372) + o2_u0 * (Sample(-0.334315788) + o2_u0 * (Sample(0.055578944)))))

// Fast implementation with Blackman window.
DELAY_COSINE_SUM(Blackman, BLACKMAN_POLY);
DELAY_COSINE_SUM_SMOOTH(Blackman, BLACKMAN_POLY);

// Fast implementation with Nuttall window.
DELAY_COSINE_SUM(Nuttall, NUTTALL_POLY);
DELAY_COSINE_SUM_SMOOTH(Nuttall, NUTTALL_POLY);

// Fast implementation with Blackman-Harris window.
DELAY_COSINE_SUM(BlackmanHarris, BLACKMANHARRIS_POLY);
DELAY_COSINE_SUM_SMOOTH(BlackmanHarris, BLACKMANHARRIS_POLY);

// Fast implementation with Blackman-Nuttall window.
DELAY_COSINE_SUM(BlackmanNuttall, BLACKMANNUTTALL_POLY);
DELAY_COSINE_SUM_SMOOTH(BlackmanNuttall, BLACKMANNUTTALL_POLY);

// Fast implementation with flat top window.
DELAY_COSINE_SUM(Flattop, FLATTOP_POLY);
DELAY_COSINE_SUM_SMOOTH(Flattop, FLATTOP_POLY);

// Fast implementation with triangle (Bartlett) window.
template<typename Sample, int maxTap = 256> class DelayTriangleBiquadSine {
private:
  static_assert(maxTap > 0 && maxTap % 2 == 0);

  static constexpr int minTimeSample = maxTap / 2 - 1;

  Sample maxTime = 0;
  Sample prevTime = 0;
  int wptr = 0;
  std::vector<Sample> buf;

public:
  std::string name() { return "fast_triangle"; }

  void setup(size_t maxTimeSample)
  {
    maxTime = Sample(maxTimeSample);
    buf.resize(std::max(size_t(maxTap), maxTimeSample + maxTap / 2 + 1));
  }

  void reset()
  {
    prevTime = 0;
    wptr = 0;
    std::fill(buf.begin(), buf.end(), Sample(0));
  }

  Sample process(Sample input, Sample timeInSample)
  {
    const int size = int(buf.size());

    // Write to buffer.
    if (++wptr >= size) wptr = 0;
    buf[wptr] = input;

    // Read from buffer.
    const int localTap = std::clamp(2 * int(timeInSample), int(2), maxTap);
    const int halfTap = localTap / 2;
    const Sample clamped = std::clamp(timeInSample, Sample(halfTap - 1), maxTime);

    const Sample timeDiff = std::abs(prevTime - clamped + Sample(1));
    prevTime = clamped;
    Sample cutoff = timeDiff <= Sample(1) ? Sample(0.5) : std::exp2(-timeDiff);

    if (timeInSample <= 0) return input * Sample(2) * cutoff;

    // Setup convolution filter parameters.
    const int timeInt = int(clamped);
    Sample fraction = clamped - Sample(timeInt);
    const Sample mid = fraction - halfTap;

    // Setup oscillator 1 (o1). Windowed sinc lowpass.
    constexpr Sample pi = std::numbers::pi_v<Sample>;
    const Sample omega = Sample(2) * pi * cutoff;
    const Sample phi = mid * omega;
    const Sample k = Sample(2) * std::cos(omega);
    Sample u1 = std::sin(phi - omega);
    Sample u2 = std::sin(phi - Sample(2) * omega);

    // Convolution.
    int rptr = wptr - timeInt - halfTap;
    if (rptr < 0) rptr += size;

    Sample sum = 0;
    const int tri_N = std::max(2, localTap + 1);
    const Sample inv_N = Sample(1) / Sample(tri_N);
    for (int i = 0; i < localTap; ++i) {
      const Sample u0 = k * u1 - u2;
      u2 = u1;
      u1 = u0;

      const Sample x = Sample(i) + mid;
      Sample sinc;
      if (std::abs(x) <= Sample(0.1)) {
        Sample q = pi * cutoff * x;
        q *= q;
        sinc = Sample(2) / Sample(3) * cutoff * (Sample(15) - Sample(7) * q)
          / (Sample(5) + q);
      } else {
        sinc = u0 / (pi * x);
      }

      const Sample window = Sample(1) - std::abs(Sample(2 * (i + 1) - tri_N) * inv_N);
      sum += window * sinc * buf[rptr];
      if (++rptr >= size) rptr = 0;
    }
    return sum;
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

template<typename Delay, typename Sample> std::string testDelayTime()
{
  Delay delay;
  delay.setup(8);

  std::vector<Sample> impulse;
  impulse.resize(32);
  impulse[0] = Sample(1);

  std::vector<Sample> sig;
  sig.resize(impulse.size());

  std::string text{"\"testDelayTime\":{"};

  constexpr Sample fraction = Sample(0.25);
  for (Sample time = fraction; time < Sample(4); ++time) {
    for (size_t i = 0; i < sig.size(); ++i) {
      sig[i] = delay.process(impulse[i], time);
    }
    delay.reset();

    text += vectorToJsonList(std::to_string(int(time)), sig);
  }
  text.pop_back();

  return text + "},";
}

template<typename Sample>
std::vector<Sample>
generateSawtooth(Sample sampleRate, Sample frequencyHz, Sample durationSecond)
{
  int period = int(sampleRate / frequencyHz);
  int durationSample = int(sampleRate * durationSecond);
  if (period == 0) {
    std::vector<Sample> zeros(durationSample);
    return zeros;
  }

  std::vector<Sample> sig(durationSample);
  Sample gain = 2 / Sample(period);
  constexpr Sample dcOffset = Sample(1);
  for (int i = 0; i < durationSample; ++i) {
    sig[i] = Sample(i % period) * gain - dcOffset;
  }
  return sig;
}

template<typename Delay, typename Sample> std::string testAntialiasing()
{
  constexpr Sample sampleRate{48000};

  Delay delay;
  delay.setup(48000);
  delay.reset();

  auto sawototh = generateSawtooth<Sample>(sampleRate, 4000, 1);

  std::vector<Sample> sig;
  sig.resize(sawototh.size());

  std::string text(std::format("\"testAntialiasing\":{{\"sampleRate\":{},", sampleRate));
  text += vectorToJsonList("source", sawototh);

  constexpr Sample pi = std::numbers::pi_v<Sample>;
  for (size_t i = 0; i < sig.size(); ++i) {
    Sample time = 16 * std::exp2(4 * std::sin(4 * pi * Sample(i) / sampleRate));
    // Sample time = 200 * std::exp2(8 * std::sin(4 * pi * Sample(i) / sampleRate));
    sig[i] = delay.process(sawototh[i], time);
  }
  text += vectorToJsonList("sinMod", sig);
  text.pop_back();

  return text + "},";
}

#define ROW_FORMAT_STR "|{:20}|{:20}|{:20}|{:20}|\n"

template<typename Delay, typename Sample> void benchmark()
{
  constexpr Sample pi = std::numbers::pi_v<Sample>;

  constexpr Sample sampleRate = Sample(48000);
  constexpr size_t nSample = 20 * size_t(sampleRate);

  Delay delay;
  delay.setup(size_t(sampleRate));
  auto rng = std::minstd_rand(3216748);
  std::uniform_real_distribution<Sample> dist(Sample(-0.25), Sample(0.25));
  double sumElapsed = 0.0;
  for (size_t i = 0; i < nSample; ++i) {
    const auto value = dist(rng);
    Sample time = Sample(200)
      * std::exp2(Sample(7.98) * std::sin(Sample(4.01) * pi * Sample(i) / sampleRate));
    auto start = std::chrono::steady_clock::now();

    delay.process(value, time);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }

  std::cout << std::format(
    ROW_FORMAT_STR, delay.name(), sumElapsed, nSample, delay.process(0, 0));
}

void testAll()
{
  std::string text{"{"};

#define ADD_DATA(DELAY_TYPE)                                                             \
  text += "\"" #DELAY_TYPE "\":{";                                                       \
  text += testDelayTime<DELAY_TYPE<double, 8>, double>();                                \
  text += testAntialiasing<DELAY_TYPE<double, 256>, double>();                           \
  text.pop_back();                                                                       \
  text += "},";

  ADD_DATA(DelayInt);
  ADD_DATA(DelayLinear);
  ADD_DATA(DelayLagrange3);
  ADD_DATA(DelayAntialiasedReference);
  ADD_DATA(DelayAntialiasedInnerProduct);
  ADD_DATA(DelayAntialiasedInlined);
  ADD_DATA(DelayAntialiasedBiquadNaive);
  ADD_DATA(DelayAntialiasedBiquad2);
  ADD_DATA(DelayAntialiasedReinsch);
  ADD_DATA(DelayAntialiasedMagic);
  ADD_DATA(DelayAntialiasedStableQuad);
  ADD_DATA(DelayLanczosABiquadSine);
  ADD_DATA(DelayLanczos1BiquadSine);
  ADD_DATA(DelayBlackmanBiquadSine);
  ADD_DATA(DelayNuttallBiquadSine);
  ADD_DATA(DelayBlackmanHarrisBiquadSine);
  ADD_DATA(DelayBlackmanNuttallBiquadSine);
  ADD_DATA(DelayFlattopBiquadSine);
  ADD_DATA(DelayBlackmanSmoothBiquadSine);
  ADD_DATA(DelayNuttallSmoothBiquadSine);
  ADD_DATA(DelayBlackmanHarrisSmoothBiquadSine);
  ADD_DATA(DelayBlackmanNuttallSmoothBiquadSine);
  ADD_DATA(DelayFlattopSmoothBiquadSine);
  ADD_DATA(DelayTriangleBiquadSine);
#undef ADD_DATA

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

  benchmark<DelayInt<double>, double>();
  benchmark<DelayLinear<double>, double>();
  benchmark<DelayLagrange3<double>, double>();
  benchmark<DelayAntialiasedReference<double>, double>();
  benchmark<DelayAntialiasedInnerProduct<double>, double>();
  benchmark<DelayAntialiasedInlined<double>, double>();
  benchmark<DelayAntialiasedBiquadNaive<double>, double>();
  benchmark<DelayAntialiasedBiquad2<double>, double>();
  benchmark<DelayAntialiasedReinsch<double>, double>();
  benchmark<DelayAntialiasedMagic<double>, double>();
  benchmark<DelayAntialiasedStableQuad<double>, double>();
  benchmark<DelayLanczosABiquadSine<double>, double>();
  benchmark<DelayLanczos1BiquadSine<double>, double>();
  benchmark<DelayBlackmanBiquadSine<double>, double>();
  benchmark<DelayNuttallBiquadSine<double>, double>();
  benchmark<DelayBlackmanHarrisBiquadSine<double>, double>();
  benchmark<DelayBlackmanNuttallBiquadSine<double>, double>();
  benchmark<DelayFlattopBiquadSine<double>, double>();
  benchmark<DelayBlackmanSmoothBiquadSine<double>, double>();
  benchmark<DelayNuttallSmoothBiquadSine<double>, double>();
  benchmark<DelayBlackmanHarrisSmoothBiquadSine<double>, double>();
  benchmark<DelayBlackmanNuttallSmoothBiquadSine<double>, double>();
  benchmark<DelayFlattopSmoothBiquadSine<double>, double>();
  benchmark<DelayTriangleBiquadSine<double>, double>();
}

int main()
{
  testAll();
  benchmarkAll();
  return 0;
}
