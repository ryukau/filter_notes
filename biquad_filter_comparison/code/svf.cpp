#include <algorithm>
#include <array>
#include <cmath>

/**
Translation of SVF in Faust filter.lib.
https://faustlibraries.grame.fr/libs/filters/#svf-filters

Arguments of `process`.
- fs             : Sampling rate.
- v0             : Input.
- normalizedFreq : Normalized frequency. In theory, range is in [0, 1].
- Q              : Quality factor, or resonance.
- shelvingGaindB : Gain for bell and shelving filters.

List of `type`.
- 0: LP
- 1: BP
- 2: HP
- 3: Notch
- 4: Peak
- 5: AP
- 6: Bell
- 7: Low-shelf
- 8: High-shelf
*/
template<typename Sample, size_t type> class SVF {
private:
  Sample ic1eq = Sample(0);
  Sample ic2eq = Sample(0);

public:
  SVF() { static_assert(type <= 8, "SVF type must be less than or equal to 8."); }

  void reset() { ic1eq = ic2eq = Sample(0); }

  Sample
  process(Sample v0, Sample normalizedFreq, Sample Q, Sample shelvingGaindB = Sample(0))
  {
    auto A = Sample(1);
    if constexpr (type >= 6) A = std::pow(Sample(10), shelvingGaindB / Sample(40));

    auto g = std::tan(normalizedFreq * Sample(pi));
    if constexpr (type == 7) {
      g /= std::sqrt(A);
    } else if (type == 8) {
      g *= std::sqrt(A);
    }

    auto k = Sample(1) / Q;
    if constexpr (type == 6) k /= A;

    // tick.
    auto v1 = (ic1eq + g * (v0 - ic2eq)) / (Sample(1) + g * (g + k));
    auto v2 = ic2eq + g * v1;

    ic1eq = Sample(2) * v1 - ic1eq;
    ic2eq = Sample(2) * v2 - ic2eq;

    // Mix.
    if constexpr (type == 0) {
      return v2;
    } else if (type == 1) {
      return v1;
    } else if (type == 2) {
      return v0 - k * v1 - v2;
    } else if (type == 3) {
      return v0 - k * v1;
    } else if (type == 4) {
      return v0 - k * v1 - Sample(2) * v2;
    } else if (type == 5) {
      return v0 - Sample(2) * k * v1;
    } else if (type == 6) {
      return v0 + k * (A * A - Sample(1)) * v1;
    } else if (type == 7) {
      return v0 + (A - Sample(1)) * k * v1 + (A * A - Sample(1)) * v2;
    } else if (type == 8) {
      return A * A * (v0 - k * v1 - v2) + A * k * v1 + v2;
    }
    return Sample(0); // Shouldn't reach here.
  }
};

int main() { return 0; }
