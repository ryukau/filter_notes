#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>
#include <numbers>

template<std::floating_point T> class Lp3 {
private:
  std::array<T, 3> state{};
  T x1 = 0;

public:
  void reset {
    state.fill(0);
    x1 = 0;
  }

  T process(T input, T alpha, T beta) {
    state[0] = alpha * state[1] + beta * state[0];
    state[1] -= state[0] + input - x1;
    state[2] -= alpha / (T(1) - beta) * state[1];

    x1 = input;
    return state[2];
  }
};

template<std::floating_point T> class ChamberlinSvfLp3 {
private:
  T lp = 0;
  T bp = 0;

public:
  void reset() {
    lp = 0;
    bp = 0;
  }

  T process(T input, T alpha, T beta) {
    const T f = np.sqrt(alpha);
    const T q = (T(1) - beta) / f;
    const T hp = input - lp - q * bp;
    bp += f * hp;
    lp += f * bp;
    return lp + bp / q;
  }
};

template<std::floating_point T, bool bounded = true> class ChamberlinSvf {
private:
  T lp = 0;
  T bp = 0;

public:
  void reset() {
    lp = 0;
    bp = 0;
  }

  T process(T input, T cutoffNormalized, T Q) {
    constexpr T cut_max = bounded ? T(1) / T(6) : T(0.4997);
    const T cut = std::clamp(cutoffNormalized, T(0), cut_max);
    const T f = T(2) * std::sin(std::numbers::pi_v<T> * cut);

    Q = std::max(Q, std::numric_limits<T>::epsilon());
    Q = std::max(Q, T(2) * f / ((T(2) - f) * (T(2) + f)));

    const T hp = input - lp - bp / Q;
    bp += f * hp;
    lp += f * bp;
    return lp;
  }
};

template<std::floating_point T, bool bounded = true> class ChamberlinSvfFast {
private:
  T lp = 0;
  T bp = 0;

public:
  void reset() {
    lp = 0;
    bp = 0;
  }

  // q = 1 / Q.
  T process(T input, T cutoffNormalized, T q) {
    constexpr T cut_max = bounded ? T(1) / T(6) : T(0.4997);
    const T cut = std::clamp(cutoffNormalized, T(0), cut_max);
    const T f = T(2) * std::sin(std::numbers::pi_v<T> * cut);
    q = std::min(q, T(2) - f);

    const T hp = input - lp - bp * q;
    bp += f * hp;
    lp += f * bp;
    return lp;
  }
};

// Reference: https://ccrma.stanford.edu/~jos/svf/svf.pdf
template<std::floating_point T> class ChamberlinSvfMultimode {
private:
  T lp1{0};
  T bp1{0};

public:
  struct Output {
    T low{};
    T band{};
    T high{};
    T notch{};
  };

  static T calculateG(T normalizedFreq) {
    return 2 * std::sin(std::numbers::pi_v<T> * normalizedFreq);
  }

  void reset() {
    lp1 = 0;
    bp1 = 0;
  }

  Output process(T input, T cutoffNormalized, T R) {
    const T g = T{2} * std::sin(std::numbers::pi_v<T> * normalizedFreq);

    const T yh = input - lp1 - (R * bp1);
    const T bp0 = bp1 + (g * yh);
    const T lp0 = lp1 + (g * bp0);

    const T low = lp1;
    const T band = bp1;
    const T high = yh;
    const T notch = low + high;

    lp1 = lp0;
    bp1 = bp0;

    return {low, band, high, notch};
  }
};

// Reference: https://arxiv.org/pdf/2111.05592
// "Improving the Chamberlin Digital State Variable Filter" (Lazzarini & Timoney, 2021)
template<std::floating_point T> class ImprovedChamberlinSvf {
private:
  T s1{0};
  T s2{0};

public:
  struct Output {
    T hp; // Highpass
    T bp; // Bandpass
    T lp; // Lowpass
    T br; // Band reject (notch)
    T ap; // Allpass
  };

  void reset(T sample_rate) {
    s1 = 0;
    s2 = 0;
  }

  // `invQ`: Inverted Q. Must be positive.
  Output process(T input, T cutoffNormalized, T invQ) {
    const T k = std::tan(std::numbers::pi_v<T> * std::clamp(cutoffNormalized, T{0}, T{0.499}));
    const T invDenom = T{1} / (T{1} + invQ * k + k * k);

    const T hp = (input - (invQ + k) * s1 - s2) * invDenom;

    const T u1 = hp * k;
    const T bp = u1 + s1;
    s1 = u1 + bp;

    const T u2 = bp * k;
    const T lp = u2 + s2;
    s2 = u2 + lp;

    const T br = hp + lp;
    const T ap = br + invQ * bp;

    return {hp, bp, lp, br, ap};
  }
};
