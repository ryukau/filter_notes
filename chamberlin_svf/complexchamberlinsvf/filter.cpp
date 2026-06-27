#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>

template<typename T> class ComplexChamberlinSvf {
public:
  using C = std::complex<T>;

private:
  C lp = {0, 0};
  C bp = {0, 0};

public:
  void reset() {
    lp = {0, 0};
    bp = {0, 0};
  }

  C process(C input, C cutoffNorm, T Q) {
    C f = std::complex<T>(2.0, 0.0) * std::sin(std::numbers::pi_v<T> * cutoffNorm);
    T f_abs = std::abs(f);
    if (f_abs > 1.0) {
      f = (f / f_abs) * T(1.0); // Normalize and scale back to 1.0
      f_abs = 1.0;
    }
    T damp_q = std::max(Q, T(2) * f_abs / ((T(2) - f_abs) * (T(2) + f_abs)));

    C hp = input - lp - (bp / damp_q);
    bp += f * hp;
    lp += f * bp;

    return lp;
  }
};
