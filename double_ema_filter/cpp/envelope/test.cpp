#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sndfile.h>
#include <string>
#include <utility>
#include <vector>

void writeWave(std::string filename, std::vector<float> &buffer, const size_t &samplerate)
{
  SF_INFO sfinfo;
  memset(&sfinfo, 0, sizeof(sfinfo));
  sfinfo.samplerate = int(samplerate);
  sfinfo.frames = buffer.size();
  sfinfo.channels = 1;
  sfinfo.format = (SF_FORMAT_WAV | SF_FORMAT_FLOAT);

  SNDFILE *file = sf_open(filename.c_str(), SFM_WRITE, &sfinfo);
  if (!file) {
    std::cout << "Error: sf_open failed." << std::endl;
    exit(EXIT_FAILURE);
  }

  size_t length = sfinfo.channels * buffer.size();
  if (sf_write_float(file, &buffer[0], length) != sf_count_t(length))
    std::cout << sf_strerror(file) << std::endl;

  if (sf_close(file) != 0) {
    std::cout << "Error: sf_close failed." << std::endl;
    exit(EXIT_FAILURE);
  }
}

constexpr double halfpi = 1.57079632679489661923;
constexpr double pi = 3.14159265358979323846;
constexpr double twopi = 6.28318530717958647692;

/**
Brent's method to find local minimum of scalar function. Translated from
`scipy.optimize.minimize_scalar`.
*/
template<typename T, typename F> std::pair<T, T> minimizeScalarBrent(F func)
{
  //////////////////////////////////////////////////////////////////
  // Bracket
  //////////////////////////////////////////////////////////////////
  constexpr T grow_limit = T(110);
  constexpr size_t maxiterBracket = 1000;

  constexpr auto _gold = T(1.618034); // golden ratio: (1.0+sqrt(5.0))/2.0
  constexpr auto _verysmall_num = std::numeric_limits<T>::epsilon();

  auto xa = T(0);
  auto xb = T(1);
  auto fa = func(xa);
  auto fb = func(xb);
  if (fa < fb) {
    std::swap(xa, xb);
    std::swap(fa, fb);
  }
  auto xc = xb + _gold * (xb - xa);
  auto fc = func(xc);
  size_t iter = 0;

  while (fc < fb && iter <= maxiterBracket) {
    auto tmp1 = (xb - xa) * (fb - fc);
    auto tmp2 = (xb - xc) * (fb - fa);
    auto val = tmp2 - tmp1;

    auto denom = std::abs(val) < _verysmall_num ? T(2) * _verysmall_num : T(2) * val;

    auto w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom;
    auto wlim = xb + grow_limit * (xc - xb);

    T fw;
    if ((w - xc) * (xb - w) > 0) {
      fw = func(w);
      if (fw < fc) {
        xa = xb;
        xb = w;
        fa = fb;
        fb = fw;
        break;
      } else if (fw > fb) {
        xc = w;
        fc = fw;
        break;
      }
      w = xc + _gold * (xc - xb);
      fw = func(w);

    } else if ((w - wlim) * (wlim - xc) >= 0) {
      w = wlim;
      fw = func(w);
    } else if ((w - wlim) * (xc - w) > 0) {
      fw = func(w);
      if (fw < fc) {
        xb = xc;
        xc = w;
        w = xc + _gold * (xc - xb);
        fb = fc;
        fc = fw;
        fw = func(w);
      }
    } else {
      w = xc + _gold * (xc - xb);
      fw = func(w);
    }
    xa = xb;
    xb = xc;
    xc = w;
    fa = fb;
    fb = fc;
    fc = fw;
  }

  //////////////////////////////////////////////////////////////////
  // Brent's algorithm
  //////////////////////////////////////////////////////////////////
  constexpr auto _mintol = T(1.0e-11);
  constexpr auto _cg = T(0.3819660);
  constexpr size_t maxiter = 64;

  auto x = xb;
  auto w = xb;
  auto v = xb;

  auto fw = func(x);
  auto fv = fw;
  auto fx = fw;

  auto a = xa;
  auto b = xc;
  if (a >= b) std::swap(a, b);

  auto deltax = T(0);
  auto rat = T(0);
  iter = 0;

  while (iter < maxiter) {
    auto tol1 = T(1.48e-8) * std::abs(x) + _mintol;
    auto tol2 = T(2) * tol1;
    auto xmid = T(0.5) * (a + b);

    // check for convergence.
    if (std::abs(x - xmid) < (tol2 - T(0.5) * (b - a))) break;

    if (std::abs(deltax) <= tol1) {
      deltax = x >= xmid ? a - x : b - x;
      rat = _cg * deltax;
    } else {
      // do a parabolic step.
      auto tmp1 = (x - w) * (fx - fv);
      auto tmp2 = (x - v) * (fx - fw);
      auto p = (x - v) * tmp2 - (x - w) * tmp1;
      tmp2 = T(2) * (tmp2 - tmp1);
      if (tmp2 > 0) p = -p;
      tmp2 = std::abs(tmp2);
      auto dx_temp = deltax;
      deltax = rat;

      // check parabolic fit.
      if (
        p > tmp2 * (a - x) && p < tmp2 * (b - x)
        && std::abs(p) < std::abs(T(0.5) * tmp2 * dx_temp))
      {
        rat = p * T(1) / tmp2; // if parabolic step is useful.
        auto u = x + rat;
        if ((u - a) < tol2 || (b - u) < tol2) rat = xmid - x >= 0 ? tol1 : -tol1;
      } else {
        deltax = x >= xmid ? a - x : b - x;
        rat = _cg * deltax;
      }
    }

    auto u = std::abs(rat) < tol1 //
      ? (rat >= 0 ? x + tol1 : x - tol1)
      : x + rat;

    auto fu = func(u); // calculate new output value

    if (fu > fx) { // if it's bigger than current
      if (u < x)
        a = u;
      else
        b = u;

      if (fu <= fw || w == x) {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    } else {
      if (u >= x)
        a = x;
      else
        b = x;

      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    }

    ++iter;
  }
  return std::pair<T, T>(x, fx);
}

// Analitic solution of DoubleEMAEnvelope output. Negated because `minimizeScalarBrent`
// finds minimum.
template<typename T> T doubleEmaEnvelopeD0Negative(T n, T k_A, T k_D)
{
  auto A = std::pow(T(1) - k_A, n + T(1)) * (k_A * n + k_A + T(1));
  auto D = std::pow(T(1) - k_D, n + T(1)) * (k_D * n + k_D + T(1));
  return (A - T(1)) * D;
}

template<typename T> T samplesToKp(T timeInSamples)
{
  if (timeInSamples < std::numeric_limits<T>::epsilon()) return T(1);
  auto y = T(1) - std::cos(T(twopi) / timeInSamples);
  return -y + std::sqrt(y * (y + T(2)));
}

template<typename Sample> class DoubleEmaADEnvelope {
public: // Make public for test.
  Sample v1_A = 0;
  Sample v2_A = 0;

  Sample v1_D = 0;
  Sample v2_D = 0;

  Sample k_A = Sample(1);
  Sample k_D = Sample(1);

  Sample invPeak = Sample(1); // Gain to normalize peak to 1.

public:
  void reset()
  {
    v1_A = 0;
    v2_A = 0;
    v1_D = 1;
    v2_D = 1;
  }

  void noteOn(Sample attackTimeSamples, Sample decayTimeSamples)
  {
    // Using `double` for minimization. `float` is inaccurate over 10^4 samples.
    auto kA = samplesToKp<double>(attackTimeSamples);
    auto kD = samplesToKp<double>(decayTimeSamples);

    if (kA == 1.0 || kD == 1.0) {
      invPeak = Sample(1);
      k_A = Sample(kA);
      k_D = Sample(kD);
    }

    auto result = minimizeScalarBrent<double>(
      [&](double n) { return doubleEmaEnvelopeD0Negative<double>(n, kA, kD); });

    std::cout << "Result(x,y) : " << result.first << ", " << result.second << "\n";

    auto peak = -result.second;
    invPeak
      = peak < std::numeric_limits<Sample>::epsilon() ? Sample(1) : Sample(1.0 / peak);
    k_A = Sample(kA);
    k_D = Sample(kD);

    reset();
  }

  Sample process()
  {
    v1_A += k_A * (Sample(1) - v1_A);
    v2_A += k_A * (v1_A - v2_A);

    v1_D += k_D * (Sample(0) - v1_D);
    v2_D += k_D * (v1_D - v2_D);

    return invPeak * v2_A * v2_D;
  }
};

int main()
{
  std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  float sampleRate = 48000;

  DoubleEmaADEnvelope<double> env;
  std::vector<float> output(sampleRate);

  env.noteOn(10 * sampleRate, 1 * sampleRate);
  for (size_t i = 0; i < output.size(); ++i) output[i] = env.process();

  auto peak = *std::max_element(output.begin(), output.end());
  std::cout << "Real Peak  : " << peak << "\n"        //
            << "invPeak    : " << env.invPeak << "\n" //
            << "Normalzied : " << peak * env.invPeak << "\n";

  writeWave("snd/output.wav", output, size_t(sampleRate));

  return 0;
}
