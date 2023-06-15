#include <array>
#include <chrono>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <numbers>
#include <numeric>
#include <random>
#include <string>
#include <type_traits>

template<typename Sample> struct ExpSmoother {
  Sample value = 0;
  Sample target = 0;
  Sample kp = 1; // Exponential moving average coefficient.

  void setCutoff(Sample sampleRate, Sample cutoffHz)
  {
    // `double` is used for accuracy.
    double y = double(1)
      - std::cos(double(2) * std::numbers::pi_v<double> * cutoffHz / sampleRate);
    kp = Sample(std::sqrt((y + double(2)) * y) - y);
  }

  void reset(Sample newTarget = 0)
  {
    value = newTarget;
    target = newTarget;
  }

  void push(Sample newTarget = 0) { target = newTarget; }
  inline Sample process() { return value += kp * (target - value); }
};

template<typename Sample> class RbjBiquad {
  Sample x1 = 0;
  Sample x2 = 0;
  Sample y1 = 0;
  Sample y2 = 0;
  std::array<ExpSmoother<Sample>, 5> co; // Coefficients {b0, b1, b2, -a1, -a2}.

public:
  const std::string name{"RbjBiquad"};

  RbjBiquad()
  {
    for (auto &x : co) x.setCutoff(Sample(1), Sample(0.02));
  }

#define ASSIGN_COEFFICINETS(method)                                                      \
  /* normalizedFreq = cutoffHz / sampleRate. */                                          \
  void method##Lowpass(Sample normalizedFreq, Sample Q)                                  \
  {                                                                                      \
    auto omega = std::numbers::pi_v<Sample> * Sample(2) * normalizedFreq;                \
    auto cos = std::cos(omega);                                                          \
    auto sin = std::sin(omega);                                                          \
    auto alpha = sin / (Sample(2) * Q);                                                  \
    auto a0_inv = Sample(1) / (Sample(1) + alpha);                                       \
    co[0].method(a0_inv *((Sample(1) - cos) / Sample(2))); /* b0 */                      \
    co[1].method(a0_inv *(Sample(1) - cos));               /* b1 */                      \
    co[2].method(a0_inv *((Sample(1) - cos) / Sample(2))); /* b2 */                      \
    co[3].method(-a0_inv *(Sample(-2) * cos));             /* a1 */                      \
    co[4].method(-a0_inv *(Sample(1) - alpha));            /* a2 */                      \
  }                                                                                      \
                                                                                         \
  void method##Allpass(Sample normalizedFreq, Sample Q)                                  \
  {                                                                                      \
    auto omega = std::numbers::pi_v<Sample> * Sample(2) * normalizedFreq;                \
    auto cos = std::cos(omega);                                                          \
    auto sin = std::sin(omega);                                                          \
    auto alpha = sin / (Sample(2) * Q);                                                  \
    auto a0_inv = Sample(1) / (Sample(1) + alpha);                                       \
    co[0].method(a0_inv *(Sample(1) - alpha));  /* b0 */                                 \
    co[1].method(a0_inv *(Sample(-2) * cos));   /* b1 */                                 \
    co[2].method(Sample(1));                    /* b2 */                                 \
    co[3].method(-a0_inv *(Sample(-2) * cos));  /* a1 */                                 \
    co[4].method(-a0_inv *(Sample(1) - alpha)); /* a2 */                                 \
  }

  ASSIGN_COEFFICINETS(reset);
  ASSIGN_COEFFICINETS(push);
#undef ASSIGN_COEFFICINETS

  void reset(Sample normalizedFreq, Sample Q, Sample value = 0)
  {
    x1 = value;
    x2 = value;
    y1 = value;
    y2 = value;

    resetLowpass(normalizedFreq, Q);
  }

  void prepare(Sample normalizedFreq, Sample Q) { pushLowpass(normalizedFreq, Q); }

  Sample process(Sample x0)
  {
    for (auto &x : co) x.process();

    auto y0 = co[0].value * x0 + co[1].value * x1 + co[2].value * x2 + co[3].value * y1
      + co[4].value * y2;

    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;

    return y0;
  }
};

template<typename Sample> class TptBiquad {
private:
  ExpSmoother<Sample> g;
  ExpSmoother<Sample> R;
  ExpSmoother<Sample> d;

  Sample s1 = 0;
  Sample s2 = 0;

public:
  const std::string name{"TptBiquad"};

  TptBiquad()
  {
    constexpr Sample cutoff = Sample(0.02);

    g.setCutoff(Sample(1), cutoff);
    R.setCutoff(Sample(1), cutoff);
    d.setCutoff(Sample(1), cutoff);
  }

#define ASSIGN_COEFFICINETS(method)                                                      \
  /* normalizedFreq = cutoffHz / sampleRate. */                                          \
  void method##Coefficient(Sample normalizedFreq, Sample Q)                              \
  {                                                                                      \
    g.method(std::tan(std::numbers::pi_v<Sample> *normalizedFreq));                      \
    R.method(Sample(1) / Q);                                                             \
    d.method(Sample(1) / (Sample(1) + g.target * g.target + g.target * R.target));       \
  }

  ASSIGN_COEFFICINETS(push);
  ASSIGN_COEFFICINETS(reset);
#undef ASSIGN_COEFFICINETS

  void reset(Sample normalizedFreq, Sample Q, Sample value = 0)
  {
    s1 = value;
    s2 = value;
    resetCoefficient(normalizedFreq, Q);
  }

  void prepare(Sample normalizedFreq, Sample Q) { pushCoefficient(normalizedFreq, Q); }

  // Lowpass.
  Sample process(Sample x0)
  {
    g.process();
    d.process();

    auto v1 = (s1 + g.value * (x0 - s2)) * d.value;
    auto v2 = s2 + g.value * v1;
    s1 = Sample(2) * v1 - s1;
    s2 = Sample(2) * v2 - s2;

    return v2;
  }
};

// Alternative implementation of TptBiquad.
// This might be faster if lp, bp, hp in `process()` are all used.
template<typename Sample> class TptSvf {
private:
  ExpSmoother<Sample> g;
  ExpSmoother<Sample> g1;
  ExpSmoother<Sample> d;

  Sample s1 = 0;
  Sample s2 = 0;

public:
  const std::string name{"TptSvf"};

  TptSvf()
  {
    constexpr Sample cutoff = Sample(0.02);

    g.setCutoff(Sample(1), cutoff);
    g1.setCutoff(Sample(1), cutoff);
    d.setCutoff(Sample(1), cutoff);
  }

#define ASSIGN_COEFFICINETS(method)                                                      \
  /* normalizedFreq = cutoffHz / sampleRate. */                                          \
  void method##Coefficient(Sample normalizedFreq, Sample Q)                              \
  {                                                                                      \
    g.method(std::tan(std::numbers::pi_v<Sample> *normalizedFreq));                      \
    auto R = Sample(1) / Q;                                                              \
    g1.method(R + g.target);                                                             \
    d.method(Sample(1) / (Sample(1) + g1.target * g.target));                            \
  }

  ASSIGN_COEFFICINETS(push);
  ASSIGN_COEFFICINETS(reset);
#undef ASSIGN_COEFFICINETS

  void reset(Sample normalizedFreq, Sample Q, Sample value = 0)
  {
    s1 = value;
    s2 = value;
    resetCoefficient(normalizedFreq, Q);
  }

  void prepare(Sample normalizedFreq, Sample Q) { pushCoefficient(normalizedFreq, Q); }

  // Lowpass.
  Sample process(Sample x)
  {
    g.process();

    auto hp = (x - g1.process() * s1 - s2) * d.process();

    auto v1 = g.value * hp;
    auto bp = v1 + s1;
    s1 = bp + v1;

    auto v2 = g.value * bp;
    auto lp = v2 + s2;
    s2 = lp + v2;

    return lp;
  }
};

template<typename Sample, size_t length> struct ExpSmootherParallel {
  std::array<Sample, length> value{};
  std::array<Sample, length> output{};
  std::array<Sample, length> target{};
  Sample kp = 1; // Exponential moving average coefficient.

  void setCutoff(Sample sampleRate, Sample cutoffHz)
  {
    // `double` is used for accuracy.
    double y = double(1)
      - std::cos(double(2) * std::numbers::pi_v<double> * cutoffHz / sampleRate);
    kp = Sample(std::sqrt((y + double(2)) * y) - y);
  }

  void reset()
  {
    value.fill(0);
    target.fill(0);
  }

  void reset(const std::array<Sample, length> &newTarget)
  {
    for (size_t i = 0; i < length; ++i) {
      value[i] = newTarget[i];
      target[i] = newTarget[i];
    }
  }

  void push(const std::array<Sample, length> &newTarget)
  {
    for (size_t i = 0; i < length; ++i) target[i] = newTarget[i];
  }

  void resetAt(size_t index, Sample newTarget = 0)
  {
    value[index] = newTarget;
    target[index] = newTarget;
  }

  void pushAt(size_t index, Sample newTarget = 0) { target[index] = newTarget; }

  void process()
  {
    for (size_t i = 0; i < length; ++i) value[i] += kp * (target[i] - value[i]);
  }
};

template<typename Sample, size_t length> class RbjBiquadParallel {
  std::array<Sample, length> x1{};
  std::array<Sample, length> x2{};
  std::array<Sample, length> y1{};
  std::array<Sample, length> y2{};

  // Coefficients {b0, b1, b2, -a1, -a2}.
  std::array<ExpSmootherParallel<Sample, length>, 5> co;

public:
  const std::string name{"RbjBiquadParallel"};

  RbjBiquadParallel()
  {
    for (auto &x : co) x.setCutoff(Sample(1), Sample(0.02));
  }

#define ASSIGN_COEFFICINETS(method)                                                      \
  /* normalizedFreq = cutoffHz / sampleRate. */                                          \
  void method##Lowpass(                                                                  \
    const std::array<Sample, length> &normalizedFreq,                                    \
    const std::array<Sample, length> &Q)                                                 \
  {                                                                                      \
    for (size_t idx = 0; idx < length; ++idx) {                                          \
      auto omega = std::numbers::pi_v<Sample> * Sample(2) * normalizedFreq[idx];         \
      auto cos = std::cos(omega);                                                        \
      auto sin = std::sin(omega);                                                        \
      auto alpha = sin / (Sample(2) * Q[idx]);                                           \
      auto a0_inv = Sample(1) / (Sample(1) + alpha);                                     \
      co[0].method##At(idx, a0_inv *((Sample(1) - cos) / Sample(2))); /* b0 */           \
      co[1].method##At(idx, a0_inv *(Sample(1) - cos));               /* b1 */           \
      co[2].method##At(idx, a0_inv *((Sample(1) - cos) / Sample(2))); /* b2 */           \
      co[3].method##At(idx, -a0_inv *(Sample(-2) * cos));             /* a1 */           \
      co[4].method##At(idx, -a0_inv *(Sample(1) - alpha));            /* a2 */           \
    }                                                                                    \
  }                                                                                      \
                                                                                         \
  void method##Allpass(                                                                  \
    const std::array<Sample, length> &normalizedFreq,                                    \
    const std::array<Sample, length> &Q)                                                 \
  {                                                                                      \
    for (size_t idx = 0; idx < length; ++idx) {                                          \
      auto omega = std::numbers::pi_v<Sample> * Sample(2) * normalizedFreq[idx];         \
      auto cos = std::cos(omega);                                                        \
      auto sin = std::sin(omega);                                                        \
      auto alpha = sin / (Sample(2) * Q[idx]);                                           \
      auto a0_inv = Sample(1) / (Sample(1) + alpha);                                     \
      co[0].method##At(idx, a0_inv *(Sample(1) - alpha));  /* b0 */                      \
      co[1].method##At(idx, a0_inv *(Sample(-2) * cos));   /* b1 */                      \
      co[2].method##At(idx, Sample(1));                    /* b2 */                      \
      co[3].method##At(idx, -a0_inv *(Sample(-2) * cos));  /* a1 */                      \
      co[4].method##At(idx, -a0_inv *(Sample(1) - alpha)); /* a2 */                      \
    }                                                                                    \
  }

  ASSIGN_COEFFICINETS(reset);
  ASSIGN_COEFFICINETS(push);
#undef ASSIGN_COEFFICINETS

  void reset(
    const std::array<Sample, length> &normalizedFreq,
    const std::array<Sample, length> &Q,
    Sample value = 0)
  {
    x1.fill(value);
    x2.fill(value);
    y1.fill(value);
    y2.fill(value);

    resetLowpass(normalizedFreq, Q);
  }

  void prepare(
    const std::array<Sample, length> &normalizedFreq, const std::array<Sample, length> &Q)
  {
    pushLowpass(normalizedFreq, Q);
  }

  void process(std::array<Sample, length> &x0)
  {
    for (size_t i = 0; i < length; ++i) {
      for (auto &x : co) x.process();

      auto y0 = co[0].value[i] * x0[i] + co[1].value[i] * x1[i] + co[2].value[i] * x2[i]
        + co[3].value[i] * y1[i] + co[4].value[i] * y2[i];

      x2[i] = x1[i];
      x1[i] = x0[i];
      y2[i] = y1[i];
      y1[i] = y0;

      x0[i] = y0;
    }
  }
};

template<typename Sample, size_t length> class TptBiquadParallel {
private:
  ExpSmootherParallel<Sample, length> g;
  ExpSmootherParallel<Sample, length> R;
  ExpSmootherParallel<Sample, length> d;

  std::array<Sample, length> s1{};
  std::array<Sample, length> s2{};

public:
  const std::string name{"TptBiquadParallel"};

  TptBiquadParallel()
  {
    constexpr Sample cutoff = Sample(0.02);
    g.setCutoff(Sample(1), cutoff);
    R.setCutoff(Sample(1), cutoff);
    d.setCutoff(Sample(1), cutoff);
  }

#define ASSIGN_COEFFICINETS(method)                                                      \
  /* normalizedFreq = cutoffHz / sampleRate. */                                          \
  void method##Coefficient(                                                              \
    const std::array<Sample, length> &normalizedFreq,                                    \
    const std::array<Sample, length> &Q)                                                 \
  {                                                                                      \
    for (size_t idx = 0; idx < length; ++idx) {                                          \
      auto gv = std::tan(std::numbers::pi_v<Sample> * normalizedFreq[idx]);              \
      g.method##At(idx, gv);                                                             \
      R.method##At(idx, Sample(1) / Q[idx]);                                             \
      d.method##At(idx, Sample(1) / (Sample(1) + gv * gv + gv * R.target[idx]));         \
    }                                                                                    \
  }

  ASSIGN_COEFFICINETS(push);
  ASSIGN_COEFFICINETS(reset);
#undef ASSIGN_COEFFICINETS

  void reset(
    const std::array<Sample, length> &normalizedFreq,
    const std::array<Sample, length> &Q,
    Sample value = 0)
  {
    for (size_t i = 0; i < length; ++i) {
      s1[i] = value;
      s2[i] = value;
      resetCoefficient(normalizedFreq, Q);
    }
  }

  void prepare(
    const std::array<Sample, length> &normalizedFreq, const std::array<Sample, length> &Q)
  {
    pushCoefficient(normalizedFreq, Q);
  }

  // Lowpass.
  void process(std::array<Sample, length> &x0)
  {
    for (size_t i = 0; i < length; ++i) {
      g.process();
      d.process();

      auto v1 = (s1[i] + g.value[i] * (x0[i] - s2[i])) * d.value[i];
      auto v2 = s2[i] + g.value[i] * v1;
      s1[i] = Sample(2) * v1 - s1[i];
      s2[i] = Sample(2) * v2 - s2[i];

      x0[i] = v2;
    }
  }
};

#define ROW_FORMAT_STR "|{:20}|{:20}|{:20}|{:20}|\n"

template<typename Filter, typename Sample> void benchSerial()
{
  constexpr Sample sampleRate = Sample(48000);
  constexpr Sample cutoffHz = Sample(1000);
  constexpr Sample Q = Sample(1) / std::numbers::sqrt2_v<Sample>;

  constexpr size_t nSample = 20 * size_t(sampleRate);

  Filter flt;
  flt.reset(cutoffHz / sampleRate, Q);
  auto rng = std::minstd_rand(3216748);
  std::uniform_real_distribution<Sample> dist(Sample(-0.25), Sample(0.25));
  double sumElapsed = 0.0;
  for (size_t i = 0; i < nSample; ++i) {
    if ((i & (512 - 1)) == 0) {
      auto scalar = Sample(nSample + i) / Sample(nSample);
      flt.prepare(scalar * cutoffHz / sampleRate, Q);
    }

    auto start = std::chrono::steady_clock::now();

    flt.process(dist(rng));

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }

  std::cout << std::format(ROW_FORMAT_STR, flt.name, sumElapsed, nSample, flt.process(0));
}

template<typename Filter, typename Sample, size_t length> void benchParallel()
{
  constexpr Sample sampleRate = Sample(48000);
  std::array<Sample, length> cutoff{};
  std::array<Sample, length> Q{};
  for (size_t idx = 0; idx < length; ++idx) {
    Sample ratio = Sample(idx) / Sample(length);
    cutoff[idx] = std::exp(-2 * ratio) * Sample(1000) / sampleRate;
    Q[idx] = std::exp2(ratio) / std::numbers::sqrt2_v<Sample>;
  }

  constexpr size_t nSample = 20 * size_t(sampleRate);
  constexpr size_t controlInterval = 512; // Must be power of 2.
  auto scalar = std::exp2(Sample(1) / Sample(1 + nSample / controlInterval));

  Filter flt;
  flt.reset(cutoff, Q);
  auto rng = std::minstd_rand(3216748);
  std::uniform_real_distribution<Sample> dist(Sample(-0.25), Sample(0.25));
  double sumElapsed = 0.0;
  std::array<Sample, length> sig{};
  for (size_t i = 0; i < nSample; ++i) {
    if ((i & (controlInterval - 1)) == 0) {
      for (size_t j = 0; j < length; ++j) {
        cutoff[j] *= scalar;
        Q[j] *= scalar;
      }
      flt.prepare(cutoff, Q);
    }

    sig.fill(dist(rng));

    auto start = std::chrono::steady_clock::now();

    flt.process(sig);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }

  std::cout << std::format(
    ROW_FORMAT_STR, flt.name, sumElapsed, nSample,
    std::accumulate(sig.begin(), sig.end(), Sample(0)));
}

int main()
{
  std::cout << "--- Warm up ";
  std::mt19937 rng(std::random_device{}());
  rng.discard(700000);
  std::cout << std::uniform_int_distribution<int>{0, 2147483647}(rng) << "\n";

  std::cout << std::format(
    ROW_FORMAT_STR, "Filter Type", "Time Elapsed [ms]", "Sample Count", "Last Output");

  using Float = float;
  constexpr size_t length = 16;

  benchSerial<RbjBiquad<Float>, Float>();
  benchSerial<TptBiquad<Float>, Float>();
  benchSerial<TptSvf<Float>, Float>();
  benchParallel<RbjBiquadParallel<Float, length>, Float, length>();
  benchParallel<TptBiquadParallel<Float, length>, Float, length>();

  return 0;
}
