#include <algorithm>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>

template<typename T> struct BracketResult {
  T xa;
  T xb;
  T xc;
  T fa;
  T fb;
  T fc;
  size_t funcalls;

  BracketResult(T xa, T xb, T xc, T fa, T fb, T fc, size_t funcalls)
    : xa(xa), xb(xb), xc(xc), fa(fa), fb(fb), fc(fc), funcalls(funcalls)
  {
  }
};

template<typename T, typename F>
BracketResult<T> bracket(F func, T xa = T(0), T xb = T(1))
{
  constexpr T grow_limit = T(110);
  constexpr size_t maxiter = 1000;

  auto _gold = T(1.618034); // golden ratio: (1.0+sqrt(5.0))/2.0
  auto _verysmall_num = std::numeric_limits<T>::epsilon();
  auto fa = func(xa);
  auto fb = func(xb);
  if (fa < fb) {
    std::swap(xa, xb);
    std::swap(fa, fb);
  }
  auto xc = xb + _gold * (xb - xa);
  auto fc = func(xc);
  size_t funcalls = 3;
  size_t iter = 0;

  while (fc < fb) {
    auto tmp1 = (xb - xa) * (fb - fc);
    auto tmp2 = (xb - xc) * (fb - fa);
    auto val = tmp2 - tmp1;

    auto denom = std::abs(val) < _verysmall_num ? T(2) * _verysmall_num : T(2) * val;

    auto w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom;
    auto wlim = xb + grow_limit * (xc - xb);

    if (iter > maxiter) break;
    ++iter;

    T fw;
    if ((w - xc) * (xb - w) > 0) {
      fw = func(w);
      funcalls += 1;
      if (fw < fc) {
        xa = xb;
        xb = w;
        fa = fb;
        fb = fw;
        return BracketResult<T>(xa, xb, xc, fa, fb, fc, funcalls);
      } else if (fw > fb) {
        xc = w;
        fc = fw;
        return BracketResult<T>(xa, xb, xc, fa, fb, fc, funcalls);
      }
      w = xc + _gold * (xc - xb);
      fw = func(w);
      funcalls += 1;
    } else if ((w - wlim) * (wlim - xc) >= 0) {
      w = wlim;
      fw = func(w);
      funcalls += 1;
    } else if ((w - wlim) * (xc - w) > 0) {
      fw = func(w);
      funcalls += 1;
      if (fw < fc) {
        xb = xc;
        xc = w;
        w = xc + _gold * (xc - xb);
        fb = fc;
        fc = fw;
        fw = func(w);
        funcalls += 1;
      }
    } else {
      w = xc + _gold * (xc - xb);
      fw = func(w);
      funcalls += 1;
    }
    xa = xb;
    xb = xc;
    xc = w;
    fa = fb;
    fb = fc;
    fc = fw;
  }
  return BracketResult<T>(xa, xb, xc, fa, fb, fc, funcalls);
}

template<typename T> struct BrentSciPyResult {
  T xmin;
  T fval;
  size_t iter;
  size_t funcalls;

  BrentSciPyResult() {}

  BrentSciPyResult(T xmin, T fval, size_t iter, size_t funcalls)
    : xmin(xmin), fval(fval), iter(iter), funcalls(funcalls)
  {
  }
};

template<typename T, typename F> BrentSciPyResult<T> brentSciPy(F func)
{
  constexpr auto _mintol = T(1.0e-11);
  constexpr auto _cg = T(0.3819660);
  constexpr size_t maxiter = 500;

  // set up for optimization
  auto bracketInfo = bracket<T>(func);
  auto &xa = bracketInfo.xa;
  auto &xb = bracketInfo.xb;
  auto &xc = bracketInfo.xc;
  auto &funcalls = bracketInfo.funcalls;

  //////////////////////////////////////////////////////////////////
  // BEGIN CORE ALGORITHM
  //////////////////////////////////////////////////////////////////
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
  funcalls += 1;
  size_t iter = 0;

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
    funcalls += 1;

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

    iter += 1;
  }
  /////////////////////////////////
  // END CORE ALGORITHM
  /////////////////////////////////

  return BrentSciPyResult<T>(x, fx, iter, funcalls);
}

template<typename T, typename F> std::pair<T, T> brentSciPySimplified(F func)
{
  //////////////////////////////////////////////////////////////////
  // Bracket
  //////////////////////////////////////////////////////////////////
  // set up for optimization
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

template<typename T> T doubleEmaEnvelopeD0Negative(T n, T k_A, T k_D)
{
  auto A = std::pow(1 - k_A, n + 1) * (k_A * n + k_A + 1);
  auto D = std::pow(1 - k_D, n + 1) * (k_D * n + k_D + 1);
  return (A - 1) * D;
}

template<typename T> T samplesToKp(T timeInSamples)
{
  constexpr T eps = std::numeric_limits<T>::epsilon();
  if (timeInSamples < eps) return T(1);
  constexpr T twopi = T(2) * std::numbers::pi_v<T>;
  auto y = 1 - std::cos(twopi / timeInSamples);
  return -y + std::sqrt(y * (y + 2));
}

template<typename T> T doubleEmaEnvelopeD1(T n, T k_A, T k_D)
{
  auto D_0 = std::pow(T(1) - k_D, n + T(1));
  auto A_0 = std::pow(T(1) - k_A, n + T(1));
  auto D_1 = k_D * n + k_D + T(1);
  auto A_1 = k_A * n + k_A + T(1);
  auto D_2 = std::log(T(1) - k_D);
  auto A_2 = std::log(T(1) - k_A);
  return D_0 * ((k_D + D_1 * D_2) * (T(1) - A_0 * A_1) - D_1 * A_0 * (k_A + A_1 * A_2));
}

int main()
{
  using boost::math::tools::bracket_and_solve_root;
  using boost::math::tools::brent_find_minima;
  using boost::math::tools::eps_tolerance;
  using boost::math::tools::toms748_solve;
  using Sample = double;

  const size_t nLoop = 1024;

  double sampleA = 0 * 48000;
  double sampleD = 0.01 * 48000;
  Sample k_A = samplesToKp(sampleA);
  Sample k_D = samplesToKp(sampleD);

  // if (sampleD == 0) peak = 0;
  // else if (sampleA == 0) peak = 1;
  // else peak = findPeak();
  std::cout << k_A << ", " << k_D << "\n";

  std::pair<Sample, Sample> result;
  double sumElapsed = 0.0;
  std::cout.precision(std::numeric_limits<Sample>::digits10);

  std::cout << "--- Warm Up\n";
  sumElapsed = 0.0;
  for (size_t n = 0; n < nLoop; ++n) {
    auto start = std::chrono::steady_clock::now();

    result = brentSciPySimplified<Sample>(
      [&](Sample n) { return doubleEmaEnvelopeD0Negative(n, k_A, k_D); });

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }
  std::cout << sumElapsed / double(nLoop) << "[ms]\n";
  std::cout << "time  : " << result.first << "\nvalue : " << result.second << std::endl;

  std::cout << "--- Boost\n";
  const int bits = std::numeric_limits<Sample>::digits;
  for (size_t n = 0; n < nLoop; ++n) {
    auto start = std::chrono::steady_clock::now();

    // // `bracket` degrades accuracy in some cases.
    // // One example case is when `sampleA = 1100000` and `sampleD = 1000000`.
    // auto func = [&](Sample n) { return doubleEmaEnvelopeD0Negative(n, k_A, k_D); };
    // auto br = bracket<Sample>(func);
    // result = brent_find_minima(func, br.xa, br.xb, bits);

    result = brent_find_minima(
      [&](Sample n) { return doubleEmaEnvelopeD0Negative(n, k_A, k_D); }, Sample(0),
      std::max(sampleA, sampleD), bits);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }
  std::cout << sumElapsed / double(nLoop) << "[ms]\n";
  std::cout << "time  : " << result.first << "\nvalue : " << result.second << std::endl;

  std::cout << "--- SciPy\n";
  BrentSciPyResult<Sample> resultSciPy;
  sumElapsed = 0.0;
  for (size_t n = 0; n < nLoop; ++n) {
    auto start = std::chrono::steady_clock::now();

    resultSciPy = brentSciPy<Sample>(
      [&](Sample n) { return doubleEmaEnvelopeD0Negative(n, k_A, k_D); });

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }
  std::cout << sumElapsed / double(nLoop) << "[ms]\n";
  std::cout << "time  : " << resultSciPy.xmin << "\nvalue : " << resultSciPy.fval
            << std::endl;

  std::cout << "--- SciPy Simplified\n";
  sumElapsed = 0.0;
  for (size_t n = 0; n < nLoop; ++n) {
    auto start = std::chrono::steady_clock::now();

    result = brentSciPySimplified<Sample>(
      [&](Sample n) { return doubleEmaEnvelopeD0Negative(n, k_A, k_D); });

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }
  std::cout << sumElapsed / double(nLoop) << "[ms]\n";
  std::cout << "time  : " << result.first << "\nvalue : " << result.second << std::endl;

  std::cout << "--- Boost TOMS 748\n";
  sumElapsed = 0.0;
  eps_tolerance<Sample> tol;
  std::uintmax_t maxiter = 64;
  Sample outputValue = 0;
  for (size_t n = 0; n < nLoop; ++n) {
    auto start = std::chrono::steady_clock::now();

    // result = toms748_solve(
    //   [&](Sample n) { return doubleEmaEnvelopeD1(n, k_A, k_D); }, Sample(0),
    //   std::max(sampleA, sampleD), tol, maxiter);
    result = bracket_and_solve_root(
      [&](Sample n) { return doubleEmaEnvelopeD1(n, k_A, k_D); },
      std::max(sampleA, sampleD) / Sample(2), Sample(2), false, tol, maxiter);
    outputValue
      = doubleEmaEnvelopeD0Negative((result.first + result.second) / Sample(2), k_A, k_D);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }
  std::cout << sumElapsed / double(nLoop) << "[ms]\n";
  std::cout << "time  : " << result.first << "\nvalue : " << outputValue << std::endl;

  return 0;
}
