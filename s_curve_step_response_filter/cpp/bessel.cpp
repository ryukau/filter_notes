/*
# Test Cases
## delay = 64
Degree 4, double:
gain: 3.6187205595758476302e-07
1, 2, 1, -1.9336785441024335608, 0.9363920816500437283,
1, 2, 1, -1.9113242790337949817, 0.91345800771590612843,

Degree 4, float:
gain: 3.6187492469252902083e-07
1, 2, 1, -1.9336786270141601563, 0.93639218807220458984,
1, 2, 1, -1.9113242626190185547, 0.91345798969268798828,

## delay = 960
Degree 4, float:
gain: 7.7036155232690362027e-12
1, 2, 1, -1.995614171028137207, 0.99562662839889526367,
1, 2, 1, -1.9939744472503662109, 0.99398434162139892578,

## delay = 8192
Degree 4, float:
gain: 1.3322676295501878485e-15
1, 2, 1, -1.9994865655899047852, 0.99948674440383911133,
1, 2, 1, -1.9992929697036743164, 0.99929308891296386719,

## delay = 48000
Degree 4, double:
gain: 1.2361157060797950178e-18
1, 2, 1, -1.9999123409645387373, 0.99991234595034472754,
1, 2, 1, -1.9998793278723459022, 0.99987933183917721003,

Degree 4, float (diverged):
gain: 0
1, 2, 1, -1.999912261962890625, 0.999912261962890625,
1, 2, 1, -1.9998793601989746094, 0.99987936019897460938,

Reference:
- [Design IIR Filters Using Cascaded Biquads - Neil
Robertson](https://www.dsprelated.com/showarticle/1137.php)
- [Design IIR Butterworth Filters Using 12 Lines of Code - Neil
Robertson](https://www.dsprelated.com/showarticle/1119.php)
*/

#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

constexpr double pi = 3.14159265358979323846264338;
constexpr double twopi = 2.0 * pi;

template<typename Sample> struct Bessel4 {
  /*
  The poles of Bessel filter becomes conjugate pair when degree is even.

  For example if degree is 6, then poles are: [p0, p1, p2, p2*, p1*, p0*].
  * is complex conjugate operator.

  ```python
  import numpy
  import scipy.signal as signal

  # z is zero, p is pole, k is gain.
  z, p, k = signal.besselap(order, norm="delay")
  ```
  */
  constexpr static uint8_t degree = 4;
  constexpr static Sample gainAp = Sample(105);
  constexpr static std::array<std::complex<Sample>, degree / 2> analogPole{{
    {Sample(-2.1037893971796273), Sample(+2.6574180418567526)},
    {Sample(-2.8962106028203722), Sample(+0.8672341289345038)},
  }};

  std::array<Sample, 2> x0{};
  std::array<Sample, 2> x1{};
  std::array<Sample, 2> x2{};
  std::array<Sample, 2> y0{};
  std::array<Sample, 2> y1{};
  std::array<Sample, 2> y2{};
  std::array<std::array<Sample, 5>, degree / 2> co;
  Sample gain = 1;

  Bessel4()
  {
    for (auto &coef : co) {
      coef[0] = 1;
      coef[1] = 2;
      coef[2] = 1;
      coef[3] = 0;
      coef[4] = 0;
    };
  }

  void reset(Sample value)
  {
    x0.fill(value);
    x1.fill(value);
    x2.fill(value);
    y0.fill(value);
    y1.fill(value);
    y2.fill(value);
  }

  // Set delay in samples.
  void setDelay(Sample delay)
  {
    auto wo = Sample(1) / delay; // Convert delay to frequency.

    constexpr auto twofs = Sample(2);
    gain = 1;
    for (uint8_t i = 0; i < co.size(); ++i) {
      std::complex<Sample> pole = wo * analogPole[i]; // Apply cutoff.
      pole = (twofs + pole) / (twofs - pole);         // Bilinear transform.
      co[i][3] = -2 * pole.real();
      co[i][4] = std::norm(pole);
      gain *= (Sample(1) + co[i][3] + co[i][4]) / Sample(4);
    }

    // debug print.
    std::cout << std::setprecision(20);
    std::cout << "gain: " << gain << "\n";
    for (uint8_t i = 0; i < co.size(); ++i) {
      for (const auto &coef : co[i]) std::cout << coef << ", ";
      std::cout << "\n";
    }
  }

  Sample process(Sample input)
  {
    x0[0] = input;
    x0[1] = y0[0];

    y0[0] = co[0][0] * x0[0] + co[0][1] * x1[0] + co[0][2] * x2[0] - co[0][3] * y1[0]
      - co[0][4] * y2[0];
    y0[1] = co[1][0] * x0[1] + co[1][1] * x1[1] + co[1][2] * x2[1] - co[1][3] * y1[1]
      - co[1][4] * y2[1];

    x2[0] = x1[0];
    x2[1] = x1[1];

    x1[0] = x0[0];
    x1[1] = x0[1];

    y2[0] = y1[0];
    y2[1] = y1[1];

    y1[0] = y0[0];
    y1[1] = y0[1];

    return y0[1];
  }
};

int main()
{
  Bessel4<double> bessel;

  bessel.setDelay(48000);

  return 0;
}
