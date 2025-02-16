#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

#include <sndfile.h>

#include "./lib/cephes.hpp"
#include "./lib/polylogarithm/Li.hpp"

void writeWave(const std::string &filename, std::vector<float> &buffer, int samplerate)
{
  SF_INFO sfinfo;
  memset(&sfinfo, 0, sizeof(sfinfo));
  sfinfo.samplerate = samplerate;
  sfinfo.frames = buffer.size();
  sfinfo.channels = 1;
  sfinfo.format = (SF_FORMAT_WAV | SF_FORMAT_FLOAT);

  SNDFILE *file = sf_open(filename.c_str(), SFM_WRITE, &sfinfo);
  if (!file) {
    std::cout << "Error: sf_open failed." << std::endl;
    exit(EXIT_FAILURE);
  }

  size_t length = sfinfo.channels * buffer.size();
  if (sf_write_float(file, &buffer[0], length) != length)
    std::cout << sf_strerror(file) << std::endl;

  if (sf_close(file) != 0) {
    std::cout << "Error: sf_close failed." << std::endl;
    exit(EXIT_FAILURE);
  }
}

template<typename Sample> class SinOsc {
public:
  Sample phase = 0;
  Sample process(Sample freqNormalized)
  {
    phase += freqNormalized;
    phase -= std::floor(phase);
    return std::sin(2 * std::numbers::pi_v<Sample> * phase);
  }
};

namespace Antiderivatives {

template<typename T> inline T hardclipJ0(T x) { return std::clamp(x, T(-1), T(1)); }
template<typename T> inline T hardclipJ1(T x)
{
  T z = std::abs(x);
  return z < 1 ? x * x / T(2) : z - T(0.5);
}
template<typename T> inline T hardclipJ2(T x)
{
  return std::abs(x) < T(1) ? x * x * x / T(6)
                            : std::copysign(x * x / T(2) + T(1) / T(6), x) - x / T(2);
}

template<typename T> inline T halfrectJ0(T x) { return x < 0 ? 0 : x; }
template<typename T> inline T halfrectJ1(T x) { return x < 0 ? 0 : x * x / T(2); }
template<typename T> inline T halfrectJ2(T x) { return x < 0 ? 0 : x * x * x / T(6); }

template<typename T> inline T powerJ0(T x, T beta = T(2.345))
{
  return std::copysign(std::pow(std::abs(x), beta), x);
}
template<typename T> inline T powerJ1(T x, T beta = T(2.345))
{
  // beta > 0.
  const T b1 = beta + T(1);
  return std::pow(std::abs(x), b1) / b1;
}
template<typename T> inline T powerJ2(T x, T beta = T(2.345))
{
  // beta > 0.
  const T b1 = beta + T(1);
  const T b2 = beta + T(2);
  const T output = std::pow(std::abs(x), b2) / (b1 * b2);
  return std::copysign(output, x);
}

// softclip2 works mostly the same as hard-cliping. However, there's a 2nd order
// (quadratic) region near threshold.
//
// - `h`: threshold of clipping.
// - `ratio`: the ratio between linear region and quadratic curve region. 0 is fully
//   linear, 1 is fully quadratic.
template<typename T> inline T softclip2J0(T x, T h = T(1), T ratio = T(0.5))
{
  const T z = std::abs(x);

  const T a1 = h * ratio;
  if (z <= a1) return x;

  const T a2 = T(2) * h - a1;
  if (z >= a2) return std::copysign(h, x);

  const T C1 = a2 - z;
  return std::copysign(h + C1 * C1 / (a1 - h) / T(4), x);
}
template<typename T> inline T softclip2J1(T x, T h = T(1), T ratio = T(0.5))
{
  const T z = std::abs(x);

  const T a1 = h * ratio;
  if (z <= a1) return z * z / T(2);

  const T a2 = T(2) * h - a1;
  const T C0 = a1 - a2;
  if (z >= a2) return a1 * (a1 / T(2) - h) + h * z + C0 * C0 * C0 / (h - a1) / T(12);

  const T C1 = z - a2;
  return a1 * (a1 / T(2) - h) + h * z + (C0 * C0 * C0 - C1 * C1 * C1) / (h - a1) / T(12);
}
template<typename T> inline T softclip2J2(T x, T h = T(1), T ratio = T(0.5))
{
  const T z = std::abs(x);

  const T a1 = h * ratio;
  if (z <= a1) return x * x * x / 6;

  const T a2 = 2 * h - a1;
  const T C0 = a1 - a2;
  const T C1 = z - a1;
  if (z >= a2) {
    const T output
      = (a1 * a1 * (T(3) * z - T(2) * a1) / T(6) + C1 * C1 * h / T(2) + (C0 * C0 * C0 * (T(4) * z - T(3) * a1 - a2)) / (h - a1) / T(48));
    return std::copysign(output, x);
  }

  const T C2 = z + a1;
  const T output
    = (a1 * a1 * (T(3) * z - T(2) * a1) / T(6) + C1 * C1 * (h / T(2) - (C2 * C2 + T(2) * C0 * C0 - T(4) * a2 * (C0 + z)) / (h - a1) / T(48)));
  return std::copysign(output, x);
}

// softclipN draws a curve that can be used as an input/output response of a compressor.
// `C` is the threshold, `R` and `beta` controls knee, `1 / S` is compression ratio.
//
// - `C`: Clipping threshold.
// - `R`: Ratio of polynomial region.
// - `beta`: Order of polynomial.
// - `S`: Slope of secondary linear region.
template<typename T>
T softclipNJ0(T x0, T C = T(1), T R = T(0.5), T beta = T(2), T S = T(0.1))
{
  const T z = std::abs(x0);

  const T rc = C * R;
  if (z <= rc) return x0;

  const T xc = rc + beta * (C - rc);
  const T A = (rc - C) / std::pow(xc - rc, beta);
  const T xs = xc - std::pow(-S / (A * beta), T(1) / (beta - T(1)));
  return z < xs ? std::copysign(A * std::pow(xc - z, beta) + C, x0)
                : std::copysign(A * std::pow(xc - xs, beta) + C + S * (z - xs), x0);
}
template<typename T>
T softclipNJ1(T x0, T C = T(1), T R = T(0.5), T beta = T(2), T S = T(0.1))
{
  const T z = std::abs(x0);

  const T rc = C * R;
  if (z <= rc) return x0 * x0 / T(2);

  const T xc = rc + beta * (C - rc);
  const T Q0 = xc - rc;
  const T A = (rc - C) / std::pow(Q0, beta);
  const T xs = xc - std::pow(-S / (A * beta), T(1) / (beta - T(1)));
  const T b1 = T(1) + beta;
  const T Q0_pow_b1 = std::pow(Q0, b1);
  return z < xs
    ? A * (Q0_pow_b1 - std::pow(xc - z, b1)) / b1 + rc * rc / T(2) + C * (z - rc)
    : (A * Q0_pow_b1 / b1 + C * Q0 + S * (z * z - xc * xc) / T(2) + rc * rc / T(2)
       + (z - xc) * (A * std::pow(xc - xs, beta) + C - S * xs));
}
template<typename T>
T softclipNJ2(T x0, T C = T(1), T R = T(0.5), T beta = T(2), T S = T(0.1))
{
  const T z = std::abs(x0);

  const T rc = C * R;
  if (z <= rc) return x0 * x0 * x0 / T(6);

  const T xc = rc + beta * (C - rc);
  const T Q0 = xc - rc;
  const T A = (rc - C) / std::pow(Q0, beta);
  const T xs = xc - std::pow(-S / (A * beta), 1 / (beta - 1));
  const T b1 = T(1) + beta;
  const T b2 = T(2) + beta;
  if (z < xs) {
    const T Q1 = z - rc;
    const T output
      = (A * ((std::pow(xc - z, b2) - std::pow(Q0, b2)) / (b1 * b2) + std::pow(Q0, b1) * Q1 / b1) + rc * rc * (z / T(2) - rc / T(3)) + C * Q1 * Q1 / T(2));
    return std::copysign(output, x0);
  }

  const T Q2 = xc - xs;
  const T Q2_pow_beta = std::pow(Q2, beta);
  const T output
    = (A * std::pow(Q0, b2) * (T(1) - T(1) / b2) / b1 + C * Q0 * Q0 / T(2) + S * (z * z * z - xc * xc * xc) / T(6) + rc * rc * (xc / T(2) - rc / T(3)) + (z - xc) * (A * (std::pow(Q0, b1) / b1 - xc * Q2_pow_beta) - C * rc + S * xc * (xs - xc / T(2)) + rc * rc / T(2) + (z + xc) / T(2) * (A * Q2_pow_beta + C - S * xs)));
  return std::copysign(output, x0);
}

template<typename T> inline T tanhJ0(T x) { return std::tanh(x); }
template<typename T> inline T tanhJ1(T x) { return std::log(std::cosh(x)); }
template<typename T> inline T tanhJ2(T x)
{
  const T e2x = std::exp(T(2) * x);
  return x * std::log(std::cosh(x) / (e2x + T(1))) + (x * x - cephes::spence(e2x)) / T(2);
}

template<typename T> inline T atanJ0(T x)
{
  return T(2) / std::numbers::pi_v<T> * std::atan(x);
}
template<typename T> inline T atanJ1(T x)
{
  return T(2) / std::numbers::pi_v<T> * (x * std::atan(x) - std::log1p(x * x) / T(2));
}
template<typename T> inline T atanJ2(T x)
{
  return (x - x * std::log1p(x * x) + (x * x - T(1)) * std::atan(x))
    / std::numbers::pi_v<T>;
}

template<typename T> inline T algebraicJ0(T x) { return x / (std::abs(x) + T(1)); }
template<typename T> inline T algebraicJ1(T x)
{
  const T z = std::abs(x);
  return z - std::log1p(z);
}
template<typename T> inline T algebraicJ2(T x)
{
  const T z = std::abs(x);
  const T w = std::log1p(z);
  return std::copysign(z * (T(1) + z / T(2) - w) - w, x);
}

template<typename T> inline T softplusJ0(T x) { return std::log(std::exp(x) + T(1)); }
template<typename T> inline T softplusJ1(T x)
{
  return T(1.6449340668482264) - cephes::spence(std::exp(x));
}
template<typename T> inline T softplusJ2(T x)
{
  return -polylogarithm::Li3(-std::exp(x));
}

template<typename T> inline T swishJ0(T x, T beta = T(2))
{
  return x == 0 ? T(0.5) : x / (std::exp(-x * beta) + T(1));
}
template<typename T> inline T swishJ1(T x, T beta = T(2))
{
  const T exb = std::exp(x * beta);
  return (x * beta * std::log1p(exb) + cephes::spence(exb) - T(1.6449340668482264))
    / (beta * beta);
}
template<typename T> inline T swishJ2(T x, T beta = T(2))
{
  const T exb = -std::exp(x * beta);
  return (T(2) * polylogarithm::Li3(exb) - x * beta * polylogarithm::Li2(exb))
    / (beta * beta * beta);
}

template<typename T> inline T exppolyJ0(T x, T beta = T(2))
{
  const T z = std::abs(x);
  return std::copysign(std::pow(z, beta) * std::exp(-z), x);
}
template<typename T> inline T exppolyJ1(T x, T beta = T(2))
{
  // Return with implicit conversion from double to T. The ported cephes functions are
  // double precision only. Letting implicit conversion to happen so that the compiler
  // warns.
  const T b1 = beta + T(1);
  return std::tgamma(b1) * (cephes::igamc(b1, 0) - cephes::igamc(b1, std::abs(x)));
}
template<typename T> inline T exppolyJ2(T x, T beta = T(2))
{
  const T z = std::abs(x);
  const T b1 = beta + T(1);
  const T b2 = beta + T(2);
  // implicit conversion from double to T on the line below.
  const T output = std::tgamma(b2) * (cephes::igamc(b2, z) - cephes::igamc(b2, 0))
    + z * std::tgamma(b1) * (cephes::igamc(b1, 0) - cephes::igamc(b1, z));
  return std::copysign(output, x);
}

template<typename T> inline T cosdecayJ0(T x)
{
  return x == 0 ? 0 : (T(1) - std::cos(x)) / x;
}
template<typename T> inline T cosdecayJ1(T x)
{
  if (x == 0) return 0;
  const T z = std::abs(x);
  const T ci = cephes::Ci(z); // implicit conversion from double to T.
  return std::log(z) - ci;
}
template<typename T> inline T cosdecayJ2(T x)
{
  if (x == 0) return 0;
  const T z = std::abs(x);
  const T ci = cephes::Ci(z); // implicit conversion from double to T.
  return std::copysign(std::sin(z) + z * (std::log(z) - ci - T(1)), x);
}

template<typename T> inline T log1pJ0(T x)
{
  return std::copysign(std::log1p(std::abs(x)), x);
}
template<typename T> inline T log1pJ1(T x)
{
  const T z = std::abs(x);
  return (z + T(1)) * std::log1p(z) - z;
}
template<typename T> inline T log1pJ2(T x)
{
  const T z = std::abs(x);
  const T z1 = z + T(1);
  return std::copysign(
    (T(2) * z1 * z1 * std::log1p(z) - (T(3) * z + T(2)) * z) / T(4), x);
}

} // namespace Antiderivatives

template<typename Sample> class SaturatorAdaa1 {
private:
  static constexpr Sample eps = Sample(std::numeric_limits<float>::epsilon());
  Sample x1 = 0;
  Sample s1 = 0;

public:
  void reset()
  {
    x1 = 0;
    s1 = 0;
  }

  template<typename Sample, typename J0, typename J1>
  Sample process(Sample input, J0 f0, J1 f1)
  {
    const auto d0 = input - x1;
    const auto s0 = f1(input);
    const auto output = (x1 == 0 && s1 == 0) || (std::abs(d0) < eps)
      ? f0(Sample(0.5) * (input + x1))
      : (s0 - s1) / d0;
    s1 = s0;
    x1 = input;
    return output;
  }

#define UHHYOU_ADAA1_ARG_0(NAME)                                                         \
  Sample NAME(Sample input)                                                              \
  {                                                                                      \
    using namespace Antiderivatives;                                                     \
    return process<Sample>(input, NAME##J0<Sample>, NAME##J1<Sample>);                   \
  }

  UHHYOU_ADAA1_ARG_0(hardclip);
  UHHYOU_ADAA1_ARG_0(halfrect);
  UHHYOU_ADAA1_ARG_0(tanh);
  UHHYOU_ADAA1_ARG_0(atan);
  UHHYOU_ADAA1_ARG_0(algebraic);
  UHHYOU_ADAA1_ARG_0(softplus);
  UHHYOU_ADAA1_ARG_0(cosdecay);
  UHHYOU_ADAA1_ARG_0(log1p);
#undef UHHYOU_ADAA1_ARG_0

#define UHHYOU_ADAA1_ARG_BETA(NAME)                                                      \
  Sample NAME(Sample input, Sample beta = Sample(2))                                     \
  {                                                                                      \
    return process<Sample>(                                                              \
      input, [&](Sample x) { return Antiderivatives::NAME##J0<Sample>(x, beta); },       \
      [&](Sample x) { return Antiderivatives::NAME##J1<Sample>(x, beta); });             \
  }

  UHHYOU_ADAA1_ARG_BETA(power);
  UHHYOU_ADAA1_ARG_BETA(swish);
  UHHYOU_ADAA1_ARG_BETA(exppoly);
#undef UHHYOU_ADAA1_ARG_BETA

  Sample softclip2(Sample input, Sample threshold = Sample(1), Sample ratio = Sample(0.5))
  {
    using namespace Antiderivatives;
    return process<Sample>(
      input, [&](Sample x) { return softclip2J0<Sample>(x, threshold, ratio); },
      [&](Sample x) { return softclip2J1<Sample>(x, threshold, ratio); });
  }

  Sample softclipN(
    Sample input,
    Sample threshold = Sample(1),
    Sample ratio = Sample(0.5),
    Sample beta = Sample(2),
    Sample slope = Sample(0.1))
  {
    using namespace Antiderivatives;
    return process<Sample>(
      input,
      [&](Sample x) { return softclipNJ0<Sample>(x, threshold, ratio, beta, slope); },
      [&](Sample x) { return softclipNJ1<Sample>(x, threshold, ratio, beta, slope); });
  }
};

template<typename Sample> class SaturatorAdaa2 {
private:
  static constexpr Sample eps = Sample(0.0000152587890625);
  Sample x1 = 0;
  Sample x2 = 0;
  Sample s1 = 0;

public:
  void reset()
  {
    x1 = 0;
    x2 = 0;
    s1 = 0;
  }

  template<typename Sample, typename J0, typename J1, typename J2>
  Sample process(Sample input, J0 f0, J1 f1, J2 f2)
  {
    const Sample f2_x1 = f2(x1);
    const Sample d0 = input - x1;
    const Sample s0
      = std::abs(d0) < eps ? f1(Sample(0.5) * (input + x1)) : (f2(input) - f2_x1) / d0;

    const Sample d1 = input - x2;
    Sample output;
    if (x1 == 0 && x2 == 0) {
      output = f0((input + Sample(2) * x1 + x2) / Sample(4));
    } else if (std::abs(d1) < eps) {
      const Sample x_bar = Sample(0.5) * (input + x2);
      const Sample delta = x_bar - x1;
      output = delta < eps
        ? f0((x_bar + x1) / Sample(2))
        : (Sample(2) / delta) * (f1(x_bar) + (f2_x1 - f2(x_bar)) / delta);
    } else {
      output = Sample(2) * (s0 - s1) / d1;
    }
    s1 = s0;
    x2 = x1;
    x1 = input;
    return output;
  }

#define UHHYOU_ADAA2_ARG_0(NAME)                                                         \
  Sample NAME(Sample input)                                                              \
  {                                                                                      \
    using namespace Antiderivatives;                                                     \
    return process<Sample>(input, NAME##J0<Sample>, NAME##J1<Sample>, NAME##J2<Sample>); \
  }

  UHHYOU_ADAA2_ARG_0(hardclip);
  UHHYOU_ADAA2_ARG_0(halfrect);
  UHHYOU_ADAA2_ARG_0(tanh);
  UHHYOU_ADAA2_ARG_0(atan);
  UHHYOU_ADAA2_ARG_0(algebraic);
  UHHYOU_ADAA2_ARG_0(softplus);
  UHHYOU_ADAA2_ARG_0(cosdecay);
  UHHYOU_ADAA2_ARG_0(log1p);
#undef UHHYOU_ADAA2_ARG_0

#define UHHYOU_ADAA2_ARG_BETA(NAME)                                                      \
  Sample NAME(Sample input, Sample beta = Sample(2))                                     \
  {                                                                                      \
    return process<Sample>(                                                              \
      input, [&](Sample x) { return Antiderivatives::NAME##J0<Sample>(x, beta); },       \
      [&](Sample x) { return Antiderivatives::NAME##J1<Sample>(x, beta); },              \
      [&](Sample x) { return Antiderivatives::NAME##J2<Sample>(x, beta); });             \
  }

  UHHYOU_ADAA2_ARG_BETA(power);
  UHHYOU_ADAA2_ARG_BETA(swish);
  UHHYOU_ADAA2_ARG_BETA(exppoly);
#undef UHHYOU_ADAA2_ARG_BETA

  Sample softclip2(Sample input, Sample threshold = Sample(1), Sample ratio = Sample(0.5))
  {
    using namespace Antiderivatives;
    return process<Sample>(
      input, [&](Sample x) { return softclip2J0<Sample>(x, threshold, ratio); },
      [&](Sample x) { return softclip2J1<Sample>(x, threshold, ratio); },
      [&](Sample x) { return softclip2J2<Sample>(x, threshold, ratio); });
  }

  Sample softclipN(
    Sample input,
    Sample threshold = Sample(1),
    Sample ratio = Sample(0.5),
    Sample beta = Sample(2),
    Sample slope = Sample(0.1))
  {
    using namespace Antiderivatives;
    return process<Sample>(
      input,
      [&](Sample x) { return softclipNJ0<Sample>(x, threshold, ratio, beta, slope); },
      [&](Sample x) { return softclipNJ1<Sample>(x, threshold, ratio, beta, slope); },
      [&](Sample x) { return softclipNJ2<Sample>(x, threshold, ratio, beta, slope); });
  }
};

template<typename SaturatorType, typename ClipBinder>
void render(std::string name, int order, ClipBinder fn)
{
  constexpr double sampleRate = double(48000);
  constexpr double freqNormalized = double(1661) / sampleRate;
  constexpr double boost = 10.0;
  constexpr double gain = 0.1;

  SinOsc<double> osc;
  SaturatorType saturator;

  std::vector<float> output(sampleRate);
  for (size_t i = 0; i < output.size(); ++i) {
    double sig = osc.process(freqNormalized);
    sig = fn(saturator, boost * sig);
    output[i] = float(gain * sig);
  }
  writeWave(std::format("{}_{}.wav", name, order), output, int(sampleRate));
}

int main()
{
#define RENDER_LEAF(NAME, ORDER)                                                         \
  render<SaturatorAdaa##ORDER<double>>(                                                  \
    #NAME, ORDER,                                                                        \
    [](SaturatorAdaa##ORDER<double> &sat, double sig) { return sat.NAME(sig); });

#define RENDER(NAME)                                                                     \
  RENDER_LEAF(NAME, 1)                                                                   \
  RENDER_LEAF(NAME, 2)

  RENDER(hardclip);
  RENDER(halfrect);
  RENDER(tanh);
  RENDER(atan);
  RENDER(algebraic);
  RENDER(softplus);
  RENDER(cosdecay);
  RENDER(log1p);
  RENDER(power);
  RENDER(swish);
  RENDER(exppoly);
  RENDER(softclip2);
  RENDER(softclipN);

  return 0;
}
