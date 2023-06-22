#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <numeric>
#include <random>
#include <string>
#include <vector>

// Lowpass and allpass. Highpass = input - lowpass.
template<typename Sample> class TptBiquad {
private:
  Sample g = 0;
  Sample R = Sample(1);
  Sample d = 0;

  Sample s1 = 0;
  Sample s2 = 0;

public:
  void reset()
  {
    s1 = 0;
    s2 = 0;
  }

  void prepare(Sample normalizedFreq, Sample Q, Sample value = 0)
  {
    s1 = value;
    s2 = value;

    g = std::tan(std::numbers::pi_v<Sample> * normalizedFreq);
    R = Sample(1) / Q;
    d = Sample(1) / (Sample(1) + g * g + g * R);
  }

#define TICK                                                                             \
  Sample v1 = (s1 + g * (v0 - s2)) * d;                                                  \
  Sample v2 = s2 + g * v1;                                                               \
  s1 = Sample(2) * v1 - s1;                                                              \
  s2 = Sample(2) * v2 - s2;

  Sample processLowpass(Sample v0)
  {
    TICK;
    return v2;
  }

  Sample processHighpass(Sample v0)
  {
    TICK;
    return v0 - R * v1 - v2;
  }

  Sample processAllpass(Sample v0)
  {
    TICK;
    return v0 - Sample(2) * R * v1;
  }
#undef TICK
};

enum class LinkwitzRileyIIRType { lowpass, highpass, allpass };

// if `isLowpass == false`, it becomes allpass.
template<typename Sample, size_t order, LinkwitzRileyIIRType filterType>
class LinkwitzRileyIIR {
private:
  static constexpr size_t minOrder = filterType == LinkwitzRileyIIRType::allpass ? 4 : 2;
  static constexpr size_t nSection = order / minOrder;
  std::array<TptBiquad<Sample>, nSection> section;

public:
  LinkwitzRileyIIR() { static_assert(order % minOrder == 0 && order >= minOrder); }

  void reset()
  {
    for (auto &x : section) x.reset();
  }

  void prepare(Sample normalizedCrossover)
  {
    constexpr Sample pi = std::numbers::pi_v<Sample>;
    for (size_t idx = 0; idx < nSection; ++idx) {
      section[idx].prepare(
        normalizedCrossover,
        Sample(0.5) / std::sin((Sample(idx) + Sample(0.5)) * pi / Sample(nSection)));
    }
  }

  Sample process(Sample x0)
  {
    if constexpr (filterType == LinkwitzRileyIIRType::lowpass) {
      for (auto &lp : section) x0 = lp.processLowpass(x0);
    } else if constexpr (filterType == LinkwitzRileyIIRType::highpass) {
      for (auto &lp : section) x0 = lp.processHighpass(x0);
    } else {
      for (auto &lp : section) x0 = lp.processAllpass(x0);
    }
    return x0;
  }
};

template<typename Sample, size_t order> class LinkwitzRileyIIR2Band4n {
private:
  LinkwitzRileyIIR<Sample, order, LinkwitzRileyIIRType::lowpass> lowpass;
  LinkwitzRileyIIR<Sample, order, LinkwitzRileyIIRType::highpass> highpass;

public:
  const std::string name{"LinkwitzRileyIIR2Band4n"};
  std::array<Sample, 2> output{}; // 0: low, 1: high.

  LinkwitzRileyIIR2Band4n() { static_assert(order % 4 == 0 && order >= 4); }

  void reset()
  {
    lowpass.reset();
    highpass.reset();
  }

  void prepare(Sample normalizedCrossover)
  {
    lowpass.prepare(normalizedCrossover);
    highpass.prepare(normalizedCrossover);
  }

  void process(Sample x0)
  {
    output[0] = lowpass.process(x0);
    output[1] = highpass.process(x0);
  }
};

// Slow.
template<typename Sample, size_t length> class RotateDelay {
private:
  std::array<Sample, length> buf{};

public:
  void reset(Sample value = 0) { buf.fill(value); }

  Sample process(Sample x0)
  {
    std::rotate(buf.begin(), buf.begin() + 1, buf.end());
    Sample output = buf[0];
    buf[0] = x0;
    return output;
  }
};

template<typename Sample, int length> class FixedIntDelay {
private:
  int ptr = 0;
  std::array<Sample, length> buf{};

public:
  void reset(Sample value = 0)
  {
    ptr = 0;
    buf.fill(value);
  }

  Sample process(Sample x0)
  {
    if (++ptr >= length) ptr = 0;
    Sample output = buf[ptr];
    buf[ptr] = x0;
    return output;
  }
};

template<typename Sample, size_t fullStage, size_t recStage> struct ComplexIIRDelay;

template<typename Sample, size_t fullStage> struct ComplexIIRDelay<Sample, fullStage, 0> {
  void reset() {}

  using C = std::complex<Sample>;
  inline C process1PoleForward(C x0, const std::array<C, fullStage> &) { return x0; }
  inline C process1PoleReversed(C x0, const std::array<C, fullStage> &) { return x0; }
};

template<typename Sample, size_t fullStage, size_t recStage> struct ComplexIIRDelay {
  static constexpr size_t index = fullStage - recStage;
  FixedIntDelay<std::complex<Sample>, size_t(1) << index> delay;
  ComplexIIRDelay<Sample, fullStage, recStage - 1> recursion;

  void reset()
  {
    delay.reset();
    recursion.reset();
  }

  inline std::complex<Sample> process1PoleForward(
    std::complex<Sample> x0, const std::array<std::complex<Sample>, fullStage> &poles)
  {
    return recursion.process1PoleForward(x0 + poles[index] * delay.process(x0), poles);
  }

  inline std::complex<Sample> process1PoleReversed(
    std::complex<Sample> x0, const std::array<std::complex<Sample>, fullStage> &poles)
  {
    return recursion.process1PoleReversed(poles[index] * x0 + delay.process(x0), poles);
  }
};

template<typename Sample, size_t stage = 8> class ComplexIIR {
private:
  static constexpr size_t nPoles = stage - 1;

  Sample a_per_b = 0;
  std::array<std::complex<Sample>, stage> poles;

  std::complex<Sample> x1 = 0;
  ComplexIIRDelay<Sample, stage, stage - 1> delay;

public:
  void reset()
  {
    x1 = 0;
    delay.reset();
  }

  void prepare(std::complex<Sample> pole)
  {
    a_per_b = pole.real() / pole.imag();
    for (auto &value : poles) {
      value = pole;
      pole *= pole;
    }
  }

  std::complex<Sample> process1PoleForward(Sample x0)
  {
    std::complex<Sample> sig = x0 + poles[0] * x1;
    x1 = x0;
    return delay.process1PoleForward(sig, poles);
  }

  std::complex<Sample> process1PoleReversed(Sample x0)
  {
    std::complex<Sample> sig = poles[0] * x0 + x1;
    x1 = x0;
    return delay.process1PoleReversed(sig, poles);
  }

  Sample process2PoleForward(Sample x0)
  {
    std::complex<Sample> sig = process1PoleForward(x0);
    return sig.real() + a_per_b * sig.imag();
  }

  Sample process2PoleReversed(Sample x0)
  {
    std::complex<Sample> sig = process1PoleReversed(x0);
    return sig.real() + a_per_b * sig.imag();
  }
};

template<typename Sample, size_t order, size_t stage = 8> class LinkwitzRileyFIR {
private:
  static constexpr size_t nSection = order / 4;

  std::array<ComplexIIR<Sample, stage>, nSection> reverse;
  std::array<ComplexIIR<Sample, stage>, nSection> forward;

  std::array<Sample, nSection> u1{};
  std::array<Sample, nSection> u2{};
  std::array<Sample, nSection> v1{};
  std::array<Sample, nSection> v2{};

  Sample gain = Sample(1);

public:
  static constexpr size_t latency = nSection * (size_t(1) << stage) + 1;

  LinkwitzRileyFIR() { static_assert(order % 4 == 0 && order >= 4); }

  void reset()
  {
    for (auto &x : reverse) x.reset();
    for (auto &x : forward) x.reset();

    u1.fill({});
    u2.fill({});
    v1.fill({});
    v2.fill({});
  }

  void prepare(Sample normalizedCrossover)
  {
    constexpr Sample pi = std::numbers::pi_v<Sample>;
    constexpr size_t N = 2 * nSection; // Butterworth order.

    gain = Sample(1);

    auto cutoffRadian = Sample(2) * pi * normalizedCrossover;
    for (size_t idx = 0; idx < nSection; ++idx) {
      auto m = Sample(2 * idx) - Sample(N) + Sample(1);
      auto analogPole = cutoffRadian * std::polar(Sample(-1), pi * m / Sample(2 * N));
      auto pole = (Sample(2) + analogPole) / (Sample(2) - analogPole);
      reverse[idx].prepare(pole);
      forward[idx].prepare(pole);
      gain *= (Sample(1) + Sample(-2) * pole.real() + std::norm(pole)) / Sample(4);
    }

    gain = std::pow(gain, Sample(1) / Sample(nSection));
  }

  Sample process(Sample x0)
  {
    constexpr Sample a1 = Sample(2); // -2 for highpass.

    for (size_t i = 0; i < nSection; ++i) {
      Sample u0 = reverse[i].process2PoleReversed(x0 * gain);
      x0 = u0 + a1 * u1[i] + u2[i];
      u2[i] = u1[i];
      u1[i] = u0;

      Sample v0 = forward[i].process2PoleForward(x0 * gain);
      x0 = v0 + a1 * v1[i] + v2[i];
      v2[i] = v1[i];
      v1[i] = v0;
    }
    return x0;
  }
};

template<typename Sample, size_t order, size_t stage = 8> class LinkwitzRileyFIR2Band4n {
private:
  LinkwitzRileyFIR<Sample, order, stage> lowpass;
  FixedIntDelay<Sample, LinkwitzRileyFIR<Sample, order, stage>::latency> highpassDelay;

public:
  const std::string name{"LinkwitzRileyFIR2Band4n"};
  std::array<Sample, 2> output{}; // 0: low, 1: high.

  LinkwitzRileyFIR2Band4n() { static_assert(order % 4 == 0 && order >= 4); }

  void reset()
  {
    lowpass.reset();
    highpassDelay.reset();
  }

  void prepare(Sample normalizedCrossover) { lowpass.prepare(normalizedCrossover); }

  void process(Sample x0)
  {
    output[0] = lowpass.process(x0);
    output[1] = highpassDelay.process(x0) - output[0];
  }
};

template<typename Sample, size_t length> class WindowedFIR {
private:
  size_t ptr = 0;
  std::array<Sample, 2 * length> buf{};
  std::array<Sample, length> fir{};

public:
  static constexpr bool isOdd = length % 2 == 1;
  static constexpr size_t latency
    = isOdd ? length / 2 + length % 2 : (length - 1) / 2 + 1;

  void reset()
  {
    ptr = 0;
    buf.fill({});
  }

  void prepare(Sample normalizedCrossover, Sample kaiserAlpha = Sample(6))
  {
    constexpr Sample pi = std::numbers::pi_v<Sample>;
    constexpr Sample K = Sample(length / 2 + length % 2);
    auto denom = std::cosh(kaiserAlpha);
    for (int idx = 0; idx < length; ++idx) {
      // Kaiser window approximation. (https://dsp.stackexchange.com/a/37720)
      Sample pos = (Sample(idx) + Sample(1) - K) / K;
      Sample window = std::cosh(kaiserAlpha * std::sqrt(Sample(1) - pos * pos)) / denom;

      // Lowpass.
      Sample x = pi * (idx + Sample(1) - K);
      Sample lp
        = x != 0 ? std::sin(2 * normalizedCrossover * x) / x : normalizedCrossover;
      fir[idx] = window * lp;
    }

    Sample sum = Sample(1) / std::accumulate(fir.begin(), fir.end(), Sample(0));
    for (auto &v : fir) v *= sum;
  }

  Sample process(Sample x0)
  {
    if (++ptr >= length) ptr = 0;
    buf[ptr] = x0;
    buf[ptr + length] = x0;
    return std::inner_product(fir.begin(), fir.end(), buf.begin() + ptr, Sample(0));
  }

  Sample highpassDelay() { return buf[ptr + latency - 1]; }
};

template<typename Sample, size_t length> class WindowedFIR2Band {
private:
  WindowedFIR<Sample, length> lowpass;

public:
  const std::string name{"WindowedFIR2Band"};
  std::array<Sample, 2> output{}; // 0: low, 1: high.

  void reset() { lowpass.reset(); }

  void prepare(Sample normalizedCrossover) { lowpass.prepare(normalizedCrossover); }

  void process(Sample x0)
  {
    output[0] = lowpass.process(x0);
    output[1] = lowpass.highpassDelay() - output[0];
  }
};

template<typename Float, size_t bands>
std::string signalToJsonString(const std::vector<std::array<Float, bands>> &sig)
{
  std::string text{"["};
  for (size_t band = 0; band < bands; ++band) {
    text += std::format("[", band);
    for (size_t idx = 0; idx < sig.size(); ++idx) {
      text += std::format("{},", sig[idx][band]);
    }
    text.pop_back();
    text += "],";
  }
  text.pop_back();
  return text + "]";
}

#define ROW_STR "|{:32}|{:20}|\n"

template<typename Filter, typename Float> auto test2band()
{
  constexpr bool isRandom = false;
  constexpr int seed = 3216741;
  constexpr Float sampleRate = 48000;
  constexpr Float cutoffHz = 1000;

  Filter flt;
  flt.reset();
  flt.prepare(cutoffHz / sampleRate);

  std::vector<std::array<Float, 2>> sig;
  sig.resize(size_t(sampleRate));
  sig[0] = {Float(1), Float(0)};

  if (isRandom) {
    auto rng = std::minstd_rand(seed);
    std::uniform_real_distribution<Float> dist(Float(-0.25), Float(0.25));
    for (auto &v : sig) v = {dist(rng), Float(0)};
  }

  double sumElapsed = 0.0;
  for (size_t i = 0; i < size_t(sampleRate); ++i) {
    auto start = std::chrono::steady_clock::now();

    flt.process(sig[i][0]);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();

    sig[i] = flt.output;
  }
  std::cout << std::format(ROW_STR, flt.name, sumElapsed);

  return std::format("\"{}\": {},", flt.name, signalToJsonString(sig));
}

int main()
{
  using Float = double;

  std::cout << "--- Warm up ";
  std::mt19937 rng(std::random_device{}());
  rng.discard(700000);
  std::cout << std::uniform_int_distribution<int>{0, 2147483647}(rng) << "\n";

  std::string text{"{"};

  std::cout << std::format(ROW_STR, "Name", "Elapsed [ms]");
  text += test2band<LinkwitzRileyIIR2Band4n<Float, 4>, Float>();
  text += test2band<LinkwitzRileyFIR2Band4n<Float, 4, 8>, Float>();
  text += test2band<WindowedFIR2Band<Float, 255>, Float>();

  text.pop_back();
  text += "}";

  std::ofstream os("ir.json");
  os << text;

  return 0;
}
