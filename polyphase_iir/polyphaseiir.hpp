#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <tuple>
#include <utility>

// Transposed direct form II.
namespace TDF2 {

template<std::floating_point T, typename Sos> class SosFilter {
  struct State {
    T s1 = 0;
    T s2 = 0;
  };
  alignas(64) std::array<State, Sos::co.size()> states{};
  T last = 0;

public:
  void reset() {
    states.fill({});
    last = 0;
  }

  inline void push(T x) {
    for (std::size_t i = 0; i < Sos::co.size(); ++i) {
      const auto& c = Sos::co[i];
      T y = c[0] * x + states[i].s1;
      states[i].s1 = states[i].s2 + c[1] * x - c[3] * y;
      states[i].s2 = c[2] * x - c[4] * y;
      x = y;
    }
    last = x;
  }

  inline T output() const { return last; }

  template<size_t nPhase> inline T process(const std::array<T, nPhase>& inputs) {
    for (const auto& x : inputs | std::views::reverse) { push(x); }
    return output();
  }
};

template<std::floating_point T, typename Sos> class SosFilter_fma {
  struct State {
    T s1 = 0;
    T s2 = 0;
  };
  alignas(64) std::array<State, Sos::co.size()> states{};
  T last = 0;

  template<typename U> static constexpr U constexpr_abs(U val) { return val < 0 ? -val : val; }

  template<std::size_t Is> inline void push_step(T& x) {
    if constexpr (Is < Sos::co.size()) {
      constexpr auto& c = Sos::co[Is];

      T y = std::fma(c[0], x, states[Is].s1);

      /*
      To minimize the rounding error, the product term (multiply part) in the inner `fma` must be
      larger than the product term in the outer `fma`. The heuristic here is to assume that the
      larger coefficent `max(c1, c3)` likely results in a larger term.
      */
      if constexpr (constexpr_abs(c[1]) > constexpr_abs(c[3])) {
        states[Is].s1 = std::fma(-c[3], y, std::fma(c[1], x, states[Is].s2));
      } else {
        states[Is].s1 = std::fma(c[1], x, std::fma(-c[3], y, states[Is].s2));
      }

      if constexpr (constexpr_abs(c[2]) > constexpr_abs(c[4])) {
        states[Is].s2 = std::fma(c[2], x, -c[4] * y);
      } else {
        states[Is].s2 = std::fma(-c[4], y, c[2] * x);
      }

      x = y;
      push_step<Is + 1>(x);
    }
  }

public:
  void reset() {
    states.fill({});
    last = 0;
  }

  inline void push(T x) {
    push_step<0>(x);
    last = x;
  }

  inline T output() const { return last; }

  template<size_t nPhase> inline T process(const std::array<T, nPhase>& inputs) {
    for (const auto& x : inputs | std::views::reverse) { push(x); }
    return output();
  }
};

} // namespace TDF2

namespace Simple {

template<std::floating_point T, auto Coefficients> class FirSection {
private:
  static constexpr std::size_t N = Coefficients.size();
  std::array<T, (N > 0 ? N - 1 : 0)> s{};

  template<std::size_t... I> inline void update_state(T input, std::index_sequence<I...>) {
    ((update_state_element<I>(input)), ...);
  }

  template<std::size_t Index> inline void update_state_element(T input) {
    constexpr auto b = Coefficients[Index + 1];
    if constexpr (Index < N - 2) {
      s[Index] = s[Index + 1] + b * input;
    } else {
      s[Index] = b * input;
    }
  }

public:
  void reset() { s.fill(T(0)); }
  inline T process(T input) {
    if constexpr (N == 0) { return T(0); }
    constexpr auto b0 = Coefficients[0];
    T y = b0 * input;
    if constexpr (N > 1) {
      y += s[0];
      update_state(input, std::make_index_sequence<N - 1>{});
    }
    return y;
  }
};

template<std::floating_point T, auto Coefficients> class IirSection {
private:
  static constexpr std::size_t nSections = Coefficients.size();
  std::array<std::array<T, 2>, nSections> s{};

  template<std::size_t... I> inline void process_sections(T& val, std::index_sequence<I...>) {
    ((process_one_section<I>(val)), ...);
  }

  template<std::size_t Index> inline void process_one_section(T& val) {
    constexpr auto a1 = Coefficients[Index][0];
    constexpr auto a2 = Coefficients[Index][1];

    T y_out = val + s[Index][0];
    s[Index][0] = s[Index][1] - a1 * y_out;
    s[Index][1] = -a2 * y_out;

    val = y_out;
  }

public:
  void reset() {
    for (auto& sec : s) { sec.fill(T(0)); }
  }

  inline T process(T input) {
    if constexpr (nSections == 0) { return input; }
    T y = input;
    process_sections(y, std::make_index_sequence<nSections>{});
    return y;
  }
};

template<std::floating_point T, typename Coefficients> class PolyphaseIir {
public:
  static constexpr int nPhase = Coefficients::nPhase;

private:
  template<std::size_t... I> static auto make_fir_tuple(std::index_sequence<I...>) {
    return std::tuple<FirSection<T, Coefficients::branches[I]>...>{};
  }
  using FirTuple = decltype(make_fir_tuple(std::make_index_sequence<nPhase>{}));

  FirTuple fir_branches;
  IirSection<T, Coefficients::denom> iir_filter;

  template<std::size_t Index> inline T process_branch(T input) {
    return std::get<Index>(fir_branches).process(input);
  }

  template<std::size_t... I>
  inline T sum_branches(const std::array<T, nPhase>& inputs, std::index_sequence<I...>) {
    return (... + process_branch<I>(inputs[I]));
  }

  template<std::size_t... I> inline void reset_branches(std::index_sequence<I...>) {
    ((std::get<I>(fir_branches).reset()), ...);
  }

public:
  void reset() {
    reset_branches(std::make_index_sequence<nPhase>{});
    iir_filter.reset();
  }

  inline T process(const std::array<T, nPhase>& inputs) {
    T fir_sum = sum_branches(inputs, std::make_index_sequence<nPhase>{});
    return iir_filter.process(fir_sum);
  }
};

} // namespace Simple

namespace KahanSum {

template<typename T, auto Coefficients> class FirSection {
private:
  static constexpr std::size_t N = Coefficients.size();
  std::array<T, (N > 0 ? N - 1 : 0)> s{};
  std::array<T, (N > 0 ? N - 1 : 0)> c{}; // Kahan compensation

  template<std::size_t... I> inline void update_state(T input, std::index_sequence<I...>) {
    ((update_state_element<I>(input)), ...);
  }

  template<std::size_t Index> inline void update_state_element(T input) {
    constexpr auto b = Coefficients[Index + 1];
    T term = b * input;

    if constexpr (Index < N - 2) {
      T next_s = s[Index + 1];
      T next_c = c[Index + 1];

      T y = term - next_c;
      T t = next_s + y;

      c[Index] = (t - next_s) - y;
      s[Index] = t;
    } else {
      s[Index] = term;
      c[Index] = T(0);
    }
  }

public:
  void reset() {
    s.fill(T(0));
    c.fill(T(0));
  }

  inline T process(T input) {
    if constexpr (N == 0) { return T(0); }

    constexpr auto b0 = Coefficients[0];
    T term = b0 * input;
    T y = term;

    if constexpr (N > 1) {
      T corrected_term = term - c[0];
      y = s[0] + corrected_term;

      update_state(input, std::make_index_sequence<N - 1>{});
    }
    return y;
  }
};

template<typename T, auto Coefficients> class IirSection {
private:
  static constexpr std::size_t nSections = Coefficients.size();
  std::array<std::array<T, 2>, nSections> s{};

  template<std::size_t... I> inline void process_sections(T& val, std::index_sequence<I...>) {
    ((process_one_section<I>(val)), ...);
  }

  template<std::size_t Index> inline void process_one_section(T& val) {
    constexpr auto a1 = Coefficients[Index][0];
    constexpr auto a2 = Coefficients[Index][1];

    T y_out = val + s[Index][0];
    s[Index][0] = s[Index][1] - a1 * y_out;
    s[Index][1] = -a2 * y_out;

    val = y_out;
  }

public:
  void reset() {
    for (auto& sec : s) { sec.fill(T(0)); }
  }

  inline T process(T input) {
    if constexpr (nSections == 0) { return input; }
    T y = input;
    process_sections(y, std::make_index_sequence<nSections>{});
    return y;
  }
};

template<typename T, typename Coefficients> class PolyphaseIir {
public:
  static constexpr int nPhase = Coefficients::nPhase;

private:
  template<std::size_t... I> static auto make_fir_tuple(std::index_sequence<I...>) {
    return std::tuple<FirSection<T, Coefficients::branches[I]>...>{};
  }
  using FirTuple = decltype(make_fir_tuple(std::make_index_sequence<nPhase>{}));

  FirTuple fir_branches;
  IirSection<T, Coefficients::denom> iir_filter;

  template<std::size_t Index> inline T process_branch(T input) {
    return std::get<Index>(fir_branches).process(input);
  }

  template<std::size_t... I>
  inline T sum_branches(const std::array<T, nPhase>& inputs, std::index_sequence<I...>) {
    T sum = 0;
    T c = 0; // Kahan compensation

    auto kahan_accumulate = [&](T input) {
      T y = input - c;
      T t = sum + y;
      c = (t - sum) - y;
      sum = t;
    };

    ((kahan_accumulate(process_branch<I>(inputs[I]))), ...);

    return sum;
  }

  template<std::size_t... I> inline void reset_branches(std::index_sequence<I...>) {
    ((std::get<I>(fir_branches).reset()), ...);
  }

public:
  void reset() {
    reset_branches(std::make_index_sequence<nPhase>{});
    iir_filter.reset();
  }

  inline T process(const std::array<T, nPhase>& inputs) {
    T fir_sum = sum_branches(inputs, std::make_index_sequence<nPhase>{});
    return iir_filter.process(fir_sum);
  }
};

} // namespace KahanSum

namespace Simple_fma {

template<std::floating_point T, auto Coefficients> class FirSection {
private:
  static constexpr std::size_t N = Coefficients.size();

  alignas(64) std::array<T, (N > 0 ? N - 1 : 0)> s{};

public:
  void reset() { s.fill(T(0)); }

  inline T process(T input) {
    if constexpr (N == 0) { return T(0); }

    constexpr auto b0 = Coefficients[0];
    if constexpr (N == 1) { return b0 * input; }

    T y = std::fma(b0, input, s[0]);
    for (std::size_t i = 0; i < N - 2; ++i) {
      s[i] = std::fma(Coefficients[i + 1], input, s[i + 1]);
    }
    s[N - 2] = Coefficients[N - 1] * input;
    return y;
  }
};

template<std::floating_point T, auto Coefficients> class IirSection {
private:
  static constexpr std::size_t nSections = Coefficients.size();
  std::array<std::array<T, 2>, nSections> s{};

  template<std::size_t... I> inline void process_sections(T& val, std::index_sequence<I...>) {
    ((process_one_section<I>(val)), ...);
  }

  template<std::size_t Index> inline void process_one_section(T& val) {
    T y_out = val + s[Index][0];
    s[Index][0] = std::fma(-Coefficients[Index][0], y_out, s[Index][1]);
    s[Index][1] = -Coefficients[Index][1] * y_out;
    val = y_out;
  }

public:
  void reset() {
    for (auto& sec : s) { sec.fill(T(0)); }
  }

  inline T process(T input) {
    if constexpr (nSections == 0) { return input; }
    T y = input;
    process_sections(y, std::make_index_sequence<nSections>{});
    return y;
  }
};

template<std::floating_point T, typename Coefficients> class PolyphaseIir {
public:
  static constexpr int nPhase = Coefficients::nPhase;

private:
  template<std::size_t... I> static auto make_fir_tuple(std::index_sequence<I...>) {
    return std::tuple<FirSection<T, Coefficients::branches[I]>...>{};
  }
  using FirTuple = decltype(make_fir_tuple(std::make_index_sequence<nPhase>{}));

  FirTuple fir_branches;
  IirSection<T, Coefficients::denom> iir_filter;

  template<std::size_t Index> inline T process_branch(T input) {
    return std::get<Index>(fir_branches).process(input);
  }

  template<std::size_t... I>
  inline T sum_branches(const std::array<T, nPhase>& inputs, std::index_sequence<I...>) {
    return (... + process_branch<I>(inputs[I]));
  }

  template<std::size_t... I> inline void reset_branches(std::index_sequence<I...>) {
    ((std::get<I>(fir_branches).reset()), ...);
  }

public:
  void reset() {
    reset_branches(std::make_index_sequence<nPhase>{});
    iir_filter.reset();
  }

  inline T process(const std::array<T, nPhase>& inputs) {
    T fir_sum = sum_branches(inputs, std::make_index_sequence<nPhase>{});
    return iir_filter.process(fir_sum);
  }
};

} // namespace Simple_fma

namespace KahanSum_fma {

template<typename T, auto Coefficients> class FirSection {
private:
  static constexpr std::size_t N = Coefficients.size();
  std::array<T, (N > 0 ? N - 1 : 0)> s{};
  std::array<T, (N > 0 ? N - 1 : 0)> c{}; // Kahan compensation

  template<std::size_t... I> inline void update_state(T input, std::index_sequence<I...>) {
    ((update_state_element<I>(input)), ...);
  }

  template<std::size_t Index> inline void update_state_element(T input) {
    constexpr auto b = Coefficients[Index + 1];

    if constexpr (Index < N - 2) {
      T next_s = s[Index + 1];
      T next_c = c[Index + 1];

      // Fuse multiplication and subtraction of previous compensation
      T y = std::fma(b, input, -next_c);
      T t = next_s + y;

      c[Index] = (t - next_s) - y;
      s[Index] = t;
    } else {
      s[Index] = b * input;
      c[Index] = T(0);
    }
  }

public:
  void reset() {
    s.fill(T(0));
    c.fill(T(0));
  }

  inline T process(T input) {
    if constexpr (N == 0) { return T(0); }

    constexpr auto b0 = Coefficients[0];

    if constexpr (N > 1) {
      // Fuse multiplication and subtraction of compensation error
      T y = std::fma(b0, input, -c[0]);
      T out_val = s[0] + y;

      update_state(input, std::make_index_sequence<N - 1>{});
      return out_val;
    } else {
      return b0 * input;
    }
  }
};

template<typename T, auto Coefficients> class IirSection {
private:
  static constexpr std::size_t nSections = Coefficients.size();
  std::array<std::array<T, 2>, nSections> s{};

  template<std::size_t... I> inline void process_sections(T& val, std::index_sequence<I...>) {
    ((process_one_section<I>(val)), ...);
  }

  template<std::size_t Index> inline void process_one_section(T& val) {
    constexpr auto a1 = Coefficients[Index][0];
    constexpr auto a2 = Coefficients[Index][1];

    T y_out = val + s[Index][0];
    s[Index][0] = std::fma(-a1, y_out, s[Index][1]);
    s[Index][1] = -a2 * y_out;

    val = y_out;
  }

public:
  void reset() {
    for (auto& sec : s) { sec.fill(T(0)); }
  }

  inline T process(T input) {
    if constexpr (nSections == 0) { return input; }
    T y = input;
    process_sections(y, std::make_index_sequence<nSections>{});
    return y;
  }
};

template<typename T, typename Coefficients> class PolyphaseIir {
public:
  static constexpr int nPhase = Coefficients::nPhase;

private:
  template<std::size_t... I> static auto make_fir_tuple(std::index_sequence<I...>) {
    return std::tuple<FirSection<T, Coefficients::branches[I]>...>{};
  }
  using FirTuple = decltype(make_fir_tuple(std::make_index_sequence<nPhase>{}));

  FirTuple fir_branches;
  IirSection<T, Coefficients::denom> iir_filter;

  template<std::size_t Index> inline T process_branch(T input) {
    return std::get<Index>(fir_branches).process(input);
  }

  template<std::size_t... I>
  inline T sum_branches(const std::array<T, nPhase>& inputs, std::index_sequence<I...>) {
    T sum = 0;
    T c = 0; // Kahan compensation

    auto kahan_accumulate = [&](T input) {
      T y = input - c;
      T t = sum + y;
      c = (t - sum) - y;
      sum = t;
    };

    ((kahan_accumulate(process_branch<I>(inputs[I]))), ...);

    return sum;
  }

  template<std::size_t... I> inline void reset_branches(std::index_sequence<I...>) {
    ((std::get<I>(fir_branches).reset()), ...);
  }

public:
  void reset() {
    reset_branches(std::make_index_sequence<nPhase>{});
    iir_filter.reset();
  }

  inline T process(const std::array<T, nPhase>& inputs) {
    T fir_sum = sum_branches(inputs, std::make_index_sequence<nPhase>{});
    return iir_filter.process(fir_sum);
  }
};

} // namespace KahanSum_fma

namespace SosSimple {

template<typename T> inline constexpr T constexpr_abs(T val) { return val < 0 ? -val : val; }

template<std::floating_point T, typename Coefficients> class PolyphaseIirSos {
public:
  static constexpr std::size_t nPhase = Coefficients::nPhase;
  static constexpr std::size_t nSections = Coefficients::nSections;

private:
  struct SectionState {
    T s1 = T(0);
    T s2 = T(0);
  };

  alignas(64) std::array<std::array<SectionState, nSections>, nPhase> states{};

  template<std::size_t PhaseIdx, std::size_t SecIdx> inline void process_section(T& x) {
    constexpr auto& co = Coefficients::co[PhaseIdx][SecIdx]; // [b0, b1, b2, a1, a2]
    auto& st = states[PhaseIdx][SecIdx];

    T y = std::fma(co[0], x, st.s1);

    if constexpr (constexpr_abs(co[1]) > constexpr_abs(co[3])) {
      st.s1 = std::fma(-co[3], y, std::fma(co[1], x, st.s2));
    } else {
      st.s1 = std::fma(co[1], x, std::fma(-co[3], y, st.s2));
    }

    if constexpr (constexpr_abs(co[2]) > constexpr_abs(co[4])) {
      st.s2 = std::fma(co[2], x, -co[4] * y);
    } else {
      st.s2 = std::fma(-co[4], y, co[2] * x);
    }

    x = y;
  }

  template<std::size_t PhaseIdx, std::size_t SecIdx> inline void process_phase_sections(T& x) {
    if constexpr (SecIdx < nSections) {
      process_section<PhaseIdx, SecIdx>(x);
      process_phase_sections<PhaseIdx, SecIdx + 1>(x);
    }
  }

  template<std::size_t PhaseIdx>
  inline void sum_phases(const std::array<T, nPhase>& inputs, T& sum) {
    if constexpr (PhaseIdx < nPhase) {
      T x = inputs[PhaseIdx];
      process_phase_sections<PhaseIdx, 0>(x);
      sum += x;
      sum_phases<PhaseIdx + 1>(inputs, sum);
    }
  }

public:
  void reset() {
    for (auto& phase_state : states) { phase_state.fill({}); }
  }

  inline T process(const std::array<T, nPhase>& inputs) {
    T sum = T(0);
    sum_phases<0>(inputs, sum);
    return sum;
  }
};

} // namespace SosSimple

namespace SosKahanSum {

template<typename T> inline constexpr T constexpr_abs(T val) { return val < 0 ? -val : val; }

template<std::floating_point T, typename Coefficients> class PolyphaseIirSos {
public:
  static constexpr std::size_t nPhase = Coefficients::nPhase;
  static constexpr std::size_t nSections = Coefficients::nSections;

private:
  struct SectionState {
    T s1 = T(0);
    T s2 = T(0);
  };

  alignas(64) std::array<std::array<SectionState, nSections>, nPhase> states{};

  template<std::size_t PhaseIdx, std::size_t SecIdx> inline void process_section(T& x) {
    constexpr auto& co = Coefficients::co[PhaseIdx][SecIdx]; // [b0, b1, b2, a1, a2]

    T y = std::fma(co[0], x, states[PhaseIdx][SecIdx].s1);

    if constexpr (constexpr_abs(co[1]) > constexpr_abs(co[3])) {
      states[PhaseIdx][SecIdx].s1
        = std::fma(-co[3], y, std::fma(co[1], x, states[PhaseIdx][SecIdx].s2));
    } else {
      states[PhaseIdx][SecIdx].s1
        = std::fma(co[1], x, std::fma(-co[3], y, states[PhaseIdx][SecIdx].s2));
    }

    if constexpr (constexpr_abs(co[2]) > constexpr_abs(co[4])) {
      states[PhaseIdx][SecIdx].s2 = std::fma(co[2], x, -co[4] * y);
    } else {
      states[PhaseIdx][SecIdx].s2 = std::fma(-co[4], y, co[2] * x);
    }

    x = y;
  }

  template<std::size_t PhaseIdx, std::size_t SecIdx> inline void process_phase_sections(T& x) {
    if constexpr (SecIdx < nSections) {
      process_section<PhaseIdx, SecIdx>(x);
      process_phase_sections<PhaseIdx, SecIdx + 1>(x);
    }
  }

  template<std::size_t PhaseIdx>
  inline void sum_phases(const std::array<T, nPhase>& inputs, T& sum, T& c) {
    if constexpr (PhaseIdx < nPhase) {
      T x = inputs[PhaseIdx];
      process_phase_sections<PhaseIdx, 0>(x);

      T y = x - c;
      T t = sum + y;
      c = (t - sum) - y;
      sum = t;

      sum_phases<PhaseIdx + 1>(inputs, sum, c);
    }
  }

public:
  void reset() {
    for (auto& phase_state : states) { phase_state.fill({}); }
  }

  inline T process(const std::array<T, nPhase>& inputs) {
    T sum = T(0);
    T c = T(0);
    sum_phases<0>(inputs, sum, c);
    return sum;
  }
};

} // namespace SosKahanSum
