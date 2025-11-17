/*
There's a naming convension to add a postfix number like `windowedSincBiquad1` and
`windowedSincBiquad2`. Following is the meaning of postfix:

- 1: Naive implementation.
- 2: Trigonometric functions are expanded for better numerical precision.

For example, `sin(a - b)` is changed to `sin(a) * cos(b) - cos(a) * sin(b)` on 2. This
slightly improves the precision when `a << b`.

Compiler flags like `/fp:fast` or `-ffast-math` breaks the numerical precision of 2.
*/

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <typeinfo>
#include <vector>

#if FP_FAST == 1
  #define FP_OPTION "fast"
#else
  #define FP_OPTION "default"
#endif

template<typename T>
std::vector<T> windowedSincReference(int length, T cutoff, T fraction)
{
  const T mid = T(length / 2 + length % 2) - fraction;

  constexpr T pi = std::numbers::pi_v<T>;
  std::vector<T> fir(length, T(0));
  for (int i = 0; i < length; ++i) {
    const T x = T(i) - mid;
    fir[i] = x == 0 ? T(2) * cutoff : std::sin(T(2) * cutoff * pi * x) / (pi * x);
  }
  return fir;
}

template<typename T> std::vector<T> windowedSincBiquad1(int length, T cutoff, T fraction)
{
  constexpr T pi = std::numbers::pi_v<T>;
  // constexpr T margin = T(1 << 0) * pi * std::numeric_limits<T>::epsilon();
  // if (fraction <= margin) {
  //   fraction = 0;
  // } else if (fraction >= T(1) - margin) {
  //   fraction = T(1);
  // }
  const T mid = fraction - T(length / 2 + length % 2);

  const T omega = T(2) * pi * cutoff;
  const T phi = mid * omega;
  const T k = T(2) * std::cos(omega);
  T u1 = std::sin(phi - omega);
  T u2 = std::sin(phi - T(2) * omega);

  std::vector<T> fir(length, T(0));
  for (int i = 0; i < length; ++i) {
    const T u0 = k * u1 - u2;
    u2 = u1;
    u1 = u0;

    const T x = T(i) + mid;
    fir[i] = x == 0 ? T(2) * cutoff : u0 / (pi * x);
  }
  return fir;
}

template<typename T> std::vector<T> windowedSincBiquad2(int length, T cutoff, T fraction)
{
  constexpr T pi = std::numbers::pi_v<T>;
  const T mid = fraction - T(length / 2 + length % 2);

  const T omega = T(2) * pi * cutoff;
  const T phi = mid * omega;
  const T cs_o = std::cos(omega);
  const T sn_o = std::sin(omega);
  const T cs_p = std::cos(phi);
  const T sn_p = std::sin(phi);
  const T k = T(2) * cs_o;
  T u1 = sn_p * cs_o - cs_p * sn_o;
  const T cs_o2 = T(2) * cs_o * cs_o - T(1);
  const T sn_o2 = T(2) * sn_o * cs_o;
  T u2 = sn_p * cs_o2 - cs_p * sn_o2;

  std::vector<T> fir(length, T(0));
  for (int i = 0; i < length; ++i) {
    const T u0 = k * u1 - u2;
    u2 = u1;
    u1 = u0;

    const T x = T(i) + mid;
    fir[i] = x == 0 ? T(2) * cutoff : u0 / (pi * x);
  }
  return fir;
}

template<typename T> std::vector<T> windowedSincReinsch1(int length, T cutoff, T fraction)
{
  constexpr T pi = std::numbers::pi_v<T>;
  const T mid = fraction - T(length / 2 + length % 2);

  const T omega = T(2) * pi * cutoff;
  const T phi = mid * omega;
  const T A = T(2) * std::sin(omega / T(2));
  const T k = A * A;
  T u = std::sin(phi - omega);
  T v = A * std::cos(phi - omega / T(2));

  std::vector<T> fir(length, T(0));
  for (int i = 0; i < length; ++i) {
    u = u + v;
    v = v - k * u;

    const T x = T(i) + mid;
    fir[i] = x == 0 ? T(2) * cutoff : u / (pi * x);
  }
  return fir;
}

template<typename T> std::vector<T> windowedSincReinsch2(int length, T cutoff, T fraction)
{
  constexpr T pi = std::numbers::pi_v<T>;
  const T mid = fraction - T(length / 2 + length % 2);

  const T omega = T(2) * pi * cutoff;
  const T phi = mid * omega;
  const T A = T(2) * std::sin(omega / T(2));
  const T k = A * A;
  const T cs_o = std::cos(omega);
  const T sn_o = std::sin(omega);
  const T cs_p = std::cos(phi);
  const T sn_p = std::sin(phi);
  T u = sn_p * cs_o - cs_p * sn_o;

  const T cs_o2 = std::cos(omega / T(2));
  const T sn_o2 = std::sin(omega / T(2));
  T v = A * (cs_p * cs_o2 + sn_p * sn_o2);

  std::vector<T> fir(length, T(0));
  for (int i = 0; i < length; ++i) {
    u = u + v;
    v = v - k * u;

    const T x = T(i) + mid;
    fir[i] = x == 0 ? T(2) * cutoff : u / (pi * x);
  }
  return fir;
}

template<typename T> std::vector<T> windowedSincMagic1(int length, T cutoff, T fraction)
{
  constexpr T pi = std::numbers::pi_v<T>;
  const T mid = fraction - T(length / 2 + length % 2);

  const T omega = T(2) * pi * cutoff;
  const T phi = mid * omega;
  const T k = T(2) * std::sin(omega / T(2));
  T u = std::cos(phi - omega * T(3) / T(2));
  T v = std::sin(phi - omega);

  std::vector<T> fir(length, T(0));
  for (int i = 0; i < length; ++i) {
    u -= k * v;
    v += k * u;

    const T x = T(i) + mid;
    fir[i] = x == 0 ? T(2) * cutoff : v / (pi * x);
  }
  return fir;
}

template<typename T> std::vector<T> windowedSincMagic2(int length, T cutoff, T fraction)
{
  constexpr T pi = std::numbers::pi_v<T>;
  const T mid = fraction - T(length / 2 + length % 2);

  const T omega = T(2) * pi * cutoff;
  const T phi = mid * omega;
  const T k = T(2) * std::sin(omega / T(2));
  const T cs_o = std::cos(omega);
  const T sn_o = std::sin(omega);
  const T cs_p = std::cos(phi);
  const T sn_p = std::sin(phi);
  T v = sn_p * cs_o - cs_p * sn_o;
  const T cs_o2 = std::cos(omega * T(3) / T(2));
  const T sn_o2 = std::sin(omega * T(3) / T(2));
  T u = cs_p * cs_o2 + sn_p * sn_o2;

  std::vector<T> fir(length, T(0));
  for (int i = 0; i < length; ++i) {
    u -= k * v;
    v += k * u;

    const T x = T(i) + mid;
    fir[i] = x == 0 ? T(2) * cutoff : v / (pi * x);
  }
  return fir;
}

template<typename T>
std::vector<T> windowedSincCoupledForm1(int length, T cutoff, T fraction)
{
  constexpr T pi = std::numbers::pi_v<T>;
  const T mid = fraction - T(length / 2 + length % 2);

  const T omega = T(2) * pi * cutoff;
  const T phi = mid * omega;
  const T k1 = std::cos(omega);
  const T k2 = std::sin(omega);
  T u = std::cos(phi - omega);
  T v = std::sin(phi - omega);

  std::vector<T> fir(length, T(0));
  for (int i = 0; i < length; ++i) {
    const T u0 = u;
    const T v0 = v;
    u = k1 * u0 - k2 * v0;
    v = k2 * u0 + k1 * v0;

    const T x = T(i) + mid;
    fir[i] = x == 0 ? T(2) * cutoff : v / (pi * x);
  }
  return fir;
}

template<typename T>
std::vector<T> windowedSincCoupledForm2(int length, T cutoff, T fraction)
{
  constexpr T pi = std::numbers::pi_v<T>;
  const T mid = fraction - T(length / 2 + length % 2);

  const T omega = T(2) * pi * cutoff;
  const T phi = mid * omega;
  const T cs_o = std::cos(omega);
  const T sn_o = std::sin(omega);
  const T cs_p = std::cos(phi);
  const T sn_p = std::sin(phi);
  const T k1 = cs_o;
  const T k2 = sn_o;
  T u = cs_p * cs_o + sn_p * sn_o;
  T v = sn_p * cs_o - cs_p * sn_o;

  std::vector<T> fir(length, T(0));
  for (int i = 0; i < length; ++i) {
    const T u0 = u;
    const T v0 = v;
    u = k1 * u0 - k2 * v0;
    v = k2 * u0 + k1 * v0;

    const T x = T(i) + mid;
    fir[i] = x == 0 ? T(2) * cutoff : v / (pi * x);
  }
  return fir;
}

template<typename T>
std::vector<T> windowedSincStableQuad1(int length, T cutoff, T fraction)
{
  constexpr T pi = std::numbers::pi_v<T>;
  const T mid = fraction - T(length / 2 + length % 2);

  const T omega = T(2) * pi * cutoff;
  const T phi = mid * omega;
  const T k1 = std::tan(omega / T(2));
  const T k2 = std::sin(omega);
  T u = std::cos(phi - omega);
  T v = std::sin(phi - omega);

  std::vector<T> fir(length, T(0));
  for (int i = 0; i < length; ++i) {
    const T w = u - k1 * v;
    v += k2 * w;
    u = w - k1 * v;

    const T x = T(i) + mid;
    fir[i] = x == 0 ? T(2) * cutoff : v / (pi * x);
  }
  return fir;
}

template<typename T>
std::vector<T> windowedSincStableQuad2(int length, T cutoff, T fraction)
{
  constexpr T pi = std::numbers::pi_v<T>;
  const T mid = fraction - T(length / 2 + length % 2);

  //
  // `k1` could use following identities but they are less accurate.
  //
  // - `const T k1 = sn_o / (T(1) + cs_o);`
  // - `const T k1 = (T(1) - sn_o) / cs_o;`
  //
  const T omega = T(2) * pi * cutoff;
  const T phi = mid * omega;
  const T k1 = std::tan(omega / T(2));
  const T cs_o = std::cos(omega);
  const T sn_o = std::sin(omega);
  const T cs_p = std::cos(phi);
  const T sn_p = std::sin(phi);
  const T k2 = sn_o;
  T u = cs_p * cs_o + sn_p * sn_o;
  T v = sn_p * cs_o - cs_p * sn_o;

  std::vector<T> fir(length, T(0));
  for (int i = 0; i < length; ++i) {
    const T w = u - k1 * v;
    v += k2 * w;
    u = w - k1 * v;

    const T x = T(i) + mid;
    fir[i] = x == 0 ? T(2) * cutoff : v / (pi * x);
  }
  return fir;
}

template<typename T>
std::string vectorToJsonList(const std::string &name, const std::vector<T> &data)
{
  std::string text = std::format("\"{}\":[", name);
  for (const auto &value : data) text += std::format("{},", value);
  text.pop_back();
  text += "],";
  return text;
}

template<typename T>
std::vector<T> maxErrorAtFractionalDelay(std::function<std::vector<T>(int, T, T)> sincFn)
{
  constexpr int length = 15;
  constexpr T cutoff = T(0.5);
  constexpr T nudge = T(0);

  std::vector<T> data;
  for (T n = -std::numeric_limits<T>::digits; n <= 0; ++n) {
    const auto fraction = std::exp2(T(n)) * (T(1) + nudge) * 2 * std::numbers::pi_v<T>;

    const auto ref = windowedSincReference(length, cutoff, fraction);
    const auto fast = sincFn(length, cutoff, fraction);

    T maxError = 0;
    for (size_t i = 0; i < length; ++i) {
      maxError = std::max(maxError, std::abs(ref[i] - fast[i]));
    }

    data.push_back(maxError);
  }
  return data;
}

template<typename T> void writeJsonMaxError()
{
  std::string text = "{";

#define ADD_DATA(NAME)                                                                   \
  text += vectorToJsonList(#NAME, maxErrorAtFractionalDelay<T>(windowedSinc##NAME<T>));

  ADD_DATA(Biquad1)
  ADD_DATA(Biquad2)
  ADD_DATA(Reinsch1)
  ADD_DATA(Reinsch2)
  ADD_DATA(Magic1)
  ADD_DATA(Magic2)
  ADD_DATA(CoupledForm1)
  ADD_DATA(CoupledForm2)
  ADD_DATA(StableQuad1)
  ADD_DATA(StableQuad2)
#undef ADD_DATA

  text.pop_back();
  text += "}";
  std::ofstream ofs(
    std::format("sinc_max_error_{}_{}.json", FP_OPTION, typeid(T).name()));
  ofs << text;
}

void writeJsonFirCoefficient()
{
  using T = double;
  constexpr int length = 15;
  const T fraction = std::exp2(T(-50)) * T(1 + 0.0) * 2 * std::numbers::pi_v<T>;

  std::vector<T> cutoff(16, T(0));
  for (size_t i = 0; i < cutoff.size(); ++i) cutoff[i] = T(0.5) / T(i + 1);

  std::string text = "{\"ref\":{";
  for (const auto &cut : cutoff) {
    text += vectorToJsonList(
      std::format("{}", cut), windowedSincReference(length, cut, fraction));
  }

#define ADD_DATA(NAME)                                                                   \
  text.pop_back();                                                                       \
  text += "},\"" #NAME "\":{";                                                           \
  for (const auto &cut : cutoff) {                                                       \
    text += vectorToJsonList(                                                            \
      std::format("{}", cut), windowedSinc##NAME(length, cut, fraction));                \
  }

  ADD_DATA(Biquad1)
  ADD_DATA(Biquad2)
  ADD_DATA(Reinsch1)
  ADD_DATA(Reinsch2)
  ADD_DATA(Magic1)
  ADD_DATA(Magic2)
  ADD_DATA(CoupledForm1)
  ADD_DATA(CoupledForm2)
  ADD_DATA(StableQuad1)
  ADD_DATA(StableQuad2)
#undef ADD_DATA

  text.pop_back();
  text += "}}";
  std::ofstream ofs("test_windowed_sinc.json");
  ofs << text;
}

template<typename T> std::vector<T> fillExp(T start, T end, size_t length)
{
  std::vector<T> y(length, T(0));
  const T logStart = std::log2(start);
  const T logEnd = std::log2(end);
  for (size_t i = 0; i < length; ++i) {
    y[i] = std::exp2(std::lerp(logStart, logEnd, T(i) / (length - 1)));
  }
  return y;
}

template<typename T> void findBranchingThreshold()
{
  constexpr T pi = std::numbers::pi_v<T>;
  constexpr int length = 256;
  constexpr int nCutoff = 16;

  auto isPadeBetter = [](int length, T fc, T frac) -> std::vector<T> {
    const T x = frac <= T(0.5) ? frac : frac - T(1);
    const T target = x == 0 ? T(2) * fc : std::sin(T(2) * fc * pi * x) / (pi * x);

    auto fast = windowedSincBiquad1(length, fc, frac);
    const int mid = length / 2 + length % 2 + (frac <= T(0.5) ? 0 : -1);
    const T errorFast = std::abs(target - fast[mid]);

    T q = pi * fc * x;
    q *= q;
    const T pade = T(2) / T(3) * fc * (T(15) - T(7) * q) / (T(5) + q);
    const T errorPade = std::abs(target - pade);

    // // debug
    // std::cout << std::format(
    //   "{:24} | {:24} | {:24} | {:24} | {:24} | {:24} \n", x, mid, target, fast[mid],
    //   pade, errorFast > errorPade);

    return {T(errorFast > errorPade), errorFast, errorPade};
  };

  auto printResult = [](int length, T fc, T branchingPoint, T e1, T e2) {
    std::cout << std::format(
      "{:>6} | {:>24} | {:>24} | {:>24} | {:>24}\n", length, fc, branchingPoint, e1, e2);
  };
  std::cout << std::format(
    "{:^6} | {:^24} | {:^24} | {:^24} | {:^24}\n", "Length", "Cutoff", "Branching Point",
    "Error (fast)", "Error (Pade)");

  auto cutoff = fillExp(T(20) / T(48000), T(0.5), nCutoff);

  // Bisection in [0, 0.5].
  for (const auto &fc : cutoff) {
    T fracL = T(0);
    T fracH = T(0.5);

    auto errorH = isPadeBetter(length, fc, fracH);
    if (errorH[0]) {
      printResult(length, fc, fracH, errorH[1], errorH[2]);
      continue;
    }

    T fracM = 0;
    std::vector<T> errorM;
    while (std::nextafter(fracL, +std::numeric_limits<T>::infinity()) < fracH) {
      fracM = (fracL + fracH) / T(2);
      errorM = isPadeBetter(length, fc, fracM);
      if (errorM[0]) {
        fracL = fracM;
      } else {
        fracH = fracM;
      }
    }
    printResult(length, fc, fracL, errorM[1], errorM[2]);
  }

  std::cout << "---\n";

  // Bisection in [0.5, 1].
  for (const auto &fc : cutoff) {
    T fracL = T(0.5);
    T fracH = T(1);

    auto errorL = isPadeBetter(length, fc, fracL);
    if (errorL[0]) {
      printResult(length, fc, fracL, errorL[1], errorL[2]);
      continue;
    }

    T fracM = 0;
    std::vector<T> errorM;
    while (std::nextafter(fracL, +std::numeric_limits<T>::infinity()) < fracH) {
      fracM = (fracL + fracH) / T(2);
      errorM = isPadeBetter(length, fc, fracM);
      if (errorM[0]) {
        fracH = fracM;
      } else {
        fracL = fracM;
      }
    }
    printResult(length, fc, 1 - fracH, errorM[1], errorM[2]);
  }
}

int main()
{
  // writeJsonFirCoefficient();
  // writeJsonMaxError<float>();
  // writeJsonMaxError<double>();
  findBranchingThreshold<double>();
  return 0;
}
