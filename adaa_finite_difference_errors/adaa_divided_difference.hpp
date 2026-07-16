#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>
#include <utility>

// Likely unusable for `order >= 3` due to the numerical instability.
template<int order, typename T, typename Funcs> class Processor {
  static_assert(order >= 1, "order must be at least 1.");

private:
  static constexpr T tolerance = T(1e-8);
  static constexpr T inv24 = T(1.0 / 24.0);

  std::array<T, order> x{};
  std::array<T, order> f{};

  template<int N> struct always_false : std::false_type {};

  template<int N> static T call_f(T val) {
    if constexpr (N == 0) {
      return Funcs::f0(val);
    } else if constexpr (N == 1) {
      return Funcs::f1(val);
    } else if constexpr (N == 2) {
      return Funcs::f2(val);
    } else if constexpr (N == 3) {
      return Funcs::f3(val);
    } else if constexpr (N == 4) {
      return Funcs::f4(val); // Extend as needed
    } else {
      static_assert(always_false<N>::value, "Funcs::fN not implemented for this order.");
      return T(0);
    }
  }

  template<int N> static T call_correction_term(T val) {
    if constexpr (N >= 0) {
      return call_f<N>(val);
    } else {
      constexpr int g_idx = -N;
      if constexpr (g_idx == 1) {
        if constexpr (requires { Funcs::g1(val); }) { return Funcs::g1(val); }
      } else if constexpr (g_idx == 2) {
        if constexpr (requires { Funcs::g2(val); }) { return Funcs::g2(val); }
      }
      static_assert(always_false<N>::value,
                    "Must not reach here, unless the order of correction term is manually extended "
                    "in `process_stage`.");
      return T(0);
    }
  }

  template<int step> T process_stage(T input, T currentVal, std::array<T, order>& next_f) {
    constexpr int TargetOrder = order - 1 - step;
    constexpr T Multiplier = T(step + 1);

    next_f[step] = currentVal;

    T d = input - x[step];
    T result;

    if (std::abs(d) < tolerance) { // Fallback branch.
      T mid = (input + x[step]) * T(0.5);
      result = call_f<TargetOrder>(mid);

      T correction = call_correction_term<TargetOrder - 2>(mid);
      result += (d * d * inv24) * correction;

    } else {
      result = Multiplier * (currentVal - f[step]) / d;
    }

    if constexpr (step + 1 < order) { return process_stage<step + 1>(input, result, next_f); }
    return result;
  }

public:
  void reset(T input = T(0)) {
    x.fill(input);

    auto fill_init
      = [&]<size_t... I>(std::index_sequence<I...>) { ((f[I] = call_f<order - I>(input)), ...); };

    fill_init(std::make_index_sequence<order>{});
  }

  T process(T input) {
    T highest_F = call_f<order>(input);

    std::array<T, order> next_f;
    T F0 = process_stage<0>(input, highest_F, next_f);

    f = next_f;

    for (int i = order - 1; i > 0; --i) { x[i] = x[i - 1]; }
    x[0] = input;

    return F0;
  }
};
