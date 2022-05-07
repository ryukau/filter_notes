#include <cmath>
#include <limits>

// In C++20, use `std::lerp` in <cmath>.
template<typename T> T lerp(T y0, T y1, T t) { return y0 + t * (y1 - y0); }

// Catmull-Rom cubic interpolation.
template<typename T> T cubicInterp(T y0, T y1, T y2, T y3, T t)
{
  T t2 = t * t;
  T c0 = y1 - y2;
  T c1 = (y2 - y0) / T(2);
  T c2 = c0 + c1;
  T c3 = c0 + c2 + (y3 - y1) / T(2);
  return c3 * t * t2 - (c2 + c3) * t2 + c1 * t + y1;
}

float cinterp(float y0, float y1, float y2, float y3, float t)
{
  return cubicInterp<float>(y0, y1, y2, y3, t);
}

// 3rd order Lagrange interpolation.
template<typename T> T lagrange3Interp(T y0, T y1, T y2, T y3, T t)
{
  T u = T(1) + t;
  T d0 = y0 - y1;
  T d1 = d0 - (y1 - y2);
  T d2 = d1 - ((y1 - y2) - (y2 - y3));
  return y0 - u * (d0 + (T(1) - u) / T(2) * (d1 + (T(2) - u) / T(3) * d2));
}

float l3interp(float y0, float y1, float y2, float y3, float t)
{
  return lagrange3Interp<float>(y0, y1, y2, y3, t);
}

// PCHIP interpolation. Monotonic between control points.
template<typename T> T pchipInterp(T y0, T y1, T y2, T y3, T t)
{
  T m0 = y1 - y0;
  T m1 = y2 - y1;
  T m2 = y3 - y2;

  T dk0 = m0 * m1 <= 0 ? 0 : T(2) * (m0 * m1) / (m0 + m1);
  T dk1 = m1 * m2 <= 0 ? 0 : T(2) * (m1 * m2) / (m1 + m2);

  T t2 = t * t;
  T c0 = y1 - y2;
  T c1 = dk0;
  T c2 = c0 + c1;
  T c3 = c0 + c2 + dk1;
  return c3 * t * t2 - (c2 + c3) * t2 + c1 * t + y1;
}

template<typename T> T akimaInterp(T y0, T y1, T y2, T y3, T y4, T y5, T t)
{
  T m0 = y1 - y0;
  T m1 = y2 - y1;
  T m2 = y3 - y2;
  T m3 = y4 - y3;
  T m4 = y5 - y4;

  T w2 = std::fabs(m1 - m0);
  T w3 = std::fabs(m2 - m1);
  T w4 = std::fabs(m3 - m2);
  T w5 = std::fabs(m4 - m3);

  T b2 = w2 + w4 == 0 ? m1 : (w4 * m1 + w2 * m2) / (w2 + w4);
  T b3 = w3 + w5 == 0 ? m2 : (w5 * m2 + w3 * m3) / (w3 + w5);

  T c2 = T(3) * m2 - T(2) * b2 - b3;

  T d2 = b2 + b3 - T(2) * m2;

  return ((t * d2 + c2) * t + b2) * t + y2;
}

template<typename T> T uniformBsplineInterp(T y0, T y1, T y2, T y3, T t)
{
  T t2 = t * t;
  T t3 = t2 * t;

  T b0 = T(1) / T(6) - t / T(2) + t2 / T(2) - t3 / T(6);
  T b1 = T(2) / T(3) - t2 + t3 / T(2);
  T b2 = T(1) / T(6) + t / T(2) + t2 / T(2) - t3 / T(2);
  T b3 = t3 / T(6);
  return b0 * y0 + b1 * y1 + b2 * y2 + b3 * y3;
}

template<typename T> T centripetalCatmullRomInterp(T y0, T y1, T y2, T y3, T t, T alpha)
{
  if (t <= 0) return y1;
  if (t >= T(1)) return y2;

  //
  // abs is only applicable for 1D case. In 2D or higher, use dot product as following.
  // ```
  // d1 = y0 - y1;
  // t1 = t0 + pow(dot(d1, d1), alpha / dimension);
  // # Same for t2 and t3.
  // ```
  //
  T t0 = 0;

  T t1 = t0 + std::pow(std::fabs(y1 - y0), alpha);
  if (t1 == t0) t1 += std::numeric_limits<T>::epsilon();

  T t2 = t1 + std::pow(std::fabs(y2 - y1), alpha);
  if (t2 == t1) t2 += t1 * std::numeric_limits<T>::epsilon();

  T t3 = t2 + std::pow(std::fabs(y3 - y2), alpha);
  if (t3 == t2) t3 += t2 * std::numeric_limits<T>::epsilon();

  // It must be t0 < t1 < t2 < t3 when reaching here. Otherwise, they cause 0 division.

  t = t1 + t * (t2 - t1);

  T A1 = (t1 - t) / (t1 - t0) * y0 + (t - t0) / (t1 - t0) * y1;
  T A2 = (t2 - t) / (t2 - t1) * y1 + (t - t1) / (t2 - t1) * y2;
  T A3 = (t3 - t) / (t3 - t2) * y2 + (t - t2) / (t3 - t2) * y3;
  T B1 = (t2 - t) / (t2 - t0) * A1 + (t - t0) / (t2 - t0) * A2;
  T B2 = (t3 - t) / (t3 - t1) * A2 + (t - t1) / (t3 - t1) * A3;
  return (t2 - t) / (t2 - t1) * B1 + (t - t1) / (t2 - t1) * B2;
}

template<typename T> inline T sinc(T n)
{
  constexpr T pi = T(3.141592653589793);
  return n == 0 ? T(1) : std::sin(pi * n) / (pi * n);
}

template<typename T> T sinc4Interp(T y0, T y1, T y2, T y3, T t)
{
  constexpr T e = T(2.718281828459045);
  T s0 = sinc(t + T(1)) / e;
  T s1 = sinc(t);
  T s2 = sinc(t - T(1));
  T s3 = sinc(t - T(2)) / e;
  return (y0 * s0 + y1 * s1 + y2 * s2 + y3 * s3) / (s0 + s1 + s2 + s3);
}
