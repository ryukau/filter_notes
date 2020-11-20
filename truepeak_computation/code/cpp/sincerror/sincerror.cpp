#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

template<typename FLOAT> void testDigits()
{
  FLOAT a_val = 1.0f;
  int a_exp;
  std::frexp(a_val, &a_exp);

  FLOAT b_val = 1.0f / std::pow(2.0f, 24.0f);
  int b_exp;
  std::frexp(b_val, &b_exp);

  std::cout << "digits: " << std::numeric_limits<FLOAT>::digits << "\n"
            << "a_exp: " << a_exp << "\n"
            << "b_exp: " << b_exp << "\n"
            << a_exp - b_exp << " > " << std::numeric_limits<FLOAT>::digits << "="
            << (a_exp - b_exp > std::numeric_limits<FLOAT>::digits) << "\n"
            << std::hexfloat << a_val << " + " << b_val << " = " << a_val + b_val;
}

constexpr double pi = 3.141592653589793238462643383279;
template<typename T> T sinc(T x) { return x != 0 ? std::sin(pi * x) / (pi * x) : T(1); }
template<typename T> T toDecibel(T x) { return T(20) * std::log10(x); }

template<typename FLOAT> void testSincWidth(uint64_t step = 1, uint64_t start = 0)
{
  uint64_t i = start;
  int exponent = 0;
  while (exponent > -std::numeric_limits<FLOAT>::digits) {
    FLOAT value = sinc(FLOAT(i));
    std::frexp(value, &exponent);
    i += step;
  }
  std::cout << i << ": " << sinc(FLOAT(i)) << ", " << exponent << "\n";
}

template<typename FLOAT, int64_t length> void testSincSum(FLOAT fraction = FLOAT(0.5))
{
  FLOAT sum = 0;
  int64_t half = length / 2;
  for (int64_t i = -half; i < half + 1; ++i) sum += std::fabs(sinc(FLOAT(i) + fraction));
  std::cout << std::setprecision(std::numeric_limits<FLOAT>::digits10 + 1) //
            << sum << " (" << toDecibel(sum) << " [dB])\n";
}

int main()
{
  // testDigits<float>();

  // testSincWidth<float>();               // 61394
  // testSincWidth<double>(1, 2223497000); // 2223497001

  constexpr int64_t fs = 48000;
  constexpr int64_t day = 24 * 60 * 60;
  testSincSum<double, fs>();           // 1 sec
  testSincSum<double, 60 * fs>();      // 1 min
  testSincSum<double, 60 * 60 * fs>(); // 1 hour
  testSincSum<double, day * fs>();     // 1 day

  return 0;
}
