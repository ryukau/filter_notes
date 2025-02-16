// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#pragma once
#include <cmath>
#include <cstdint>

namespace polylogarithm {

constexpr double PI = 3.14159265358979324;
constexpr double PI2 = 2 * PI;

inline constexpr bool is_even(int64_t n) noexcept { return n % 2 == 0; }

} // namespace polylogarithm
