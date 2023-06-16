#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

enum class MatchedBiquadType {
  orfanidisPeaking,
  massbergLowpass,
  matchedLowpass,
  matchedHighpass,
  matchedBandpass,
  matchedPeaking,
  simpleMatchedLowpass,
  simpleMatchedHighpass,
  simpleMatchedBandpass,
  matchedShelvingOnePole,
};

std::vector<std::string> matchedBiquadNames{
  "orfanidisPeaking",       "massbergLowpass",       "matchedLowpass",
  "matchedHighpass",        "matchedBandpass",       "matchedPeaking",
  "simpleMatchedLowpass",   "simpleMatchedHighpass", "simpleMatchedBandpass",
  "matchedShelvingOnePole",
};

template<typename T> class MatchedBiquad {
private:
  static constexpr T pi = std::numbers::pi_v<T>;

  T x1 = 0;
  T x2 = 0;
  T y1 = 0;
  T y2 = 0;
  std::array<T, 5> co{}; // {b0, b1, b2, -a1, -a2}.

public:
  // cutoffRadian: Cutoff frequency in radian.
  // Q: Quality factor, or resonance.
  // G: Peak gain.
  void setOrfanidisPeaking(T cutoffRadian, T Q, T G)
  {
    // w <- small omega, W <- capital omega, in referenced paper.

    const auto &w0 = cutoffRadian;
    auto G0 = T(1);   // DC gain.
    auto G1 = T(1);   // Nyquist gain.
    auto dw = w0 / Q; // Bandwidth.
    auto GB = std::sqrt(G);

    auto W0 = std::tan(T(0.5) * w0);

    auto G_2 = G * G;
    auto GB_2 = GB * GB;
    auto G0_2 = G0 * G0;
    auto G1_2 = G1 * G1;

    auto W_2 = std::sqrt((G_2 - G1_2) / (G_2 - G0_2)) * W0 * W0;
    auto dW
      = (T(1) + std::sqrt((GB_2 - G0_2) / (GB_2 - G1_2)) * W_2) * std::tan(T(0.5) * dw);

    auto C = (dW * dW) * std::abs(GB_2 - G1_2)
      - T(2) * W_2
        * (std::abs(GB_2 - G0 * G1) - std::sqrt((GB_2 - G0_2) * (GB_2 - G1_2)));
    auto D
      = T(2) * W_2 * (std::abs(G_2 - G0 * G1) - std::sqrt((G_2 - G0_2) * (G_2 - G1_2)));

    auto A = std::sqrt((C + D) / std::abs(G_2 - GB_2));
    auto B = std::sqrt((G_2 * C + GB_2 * D) / std::abs(G_2 - GB_2));

    auto a0_inv = T(1) / (T(1) + W_2 + A);
    co[0] = a0_inv * (G1 + G0 * W_2 + B);     // b0
    co[1] = a0_inv * T(-2) * (G1 - G0 * W_2); // b1
    co[2] = a0_inv * (G1 + G0 * W_2 - B);     // b2
    co[3] = -a0_inv * T(-2) * (T(1) - W_2);   // -a1
    co[4] = -a0_inv * (T(1) + W_2 - A);       // -a2
  }

  void setMassbergLowpass(T cutoffRadian, T Q)
  {
    const auto &w0 = cutoffRadian;

    auto Q2 = Q * Q;
    auto t1 = pi * pi / (w0 * w0);
    auto t2 = T(1) - t1;
    auto t3 = t1 / Q2;
    auto g1 = T(1) / (std::sqrt(t2 * t2 + t3 * t3));

    T Ws = T(0);
    if (Q > std::numbers::sqrt2_v<T> / T(2)) {
      auto gr = T(2) * Q2 / std::sqrt(T(4) * Q2 - T(1));
      auto w_r = w0 * std::sqrt(T(1) - T(1) / (T(2) * Q2));
      Ws = std::tan(w_r / T(2))
        * std::pow((gr * gr - g1 * g1) / (gr * gr - T(1)), T(0.25));
    } else {
      auto w_m = w0
        * std::sqrt(
                   T(1) - T(0.5) / Q2
                   + std::sqrt((T(1) - T(4) * Q2) / (T(4) * Q2 * Q2) + T(1) / g1));
      Ws
        = std::min(T(0.5) * w0 * std::pow(T(1) - g1 * g1, T(0.25)), std::tan(w_m / T(2)));
    }

    auto w_z = T(2) * std::atan(Ws / std::sqrt(g1));
    auto z_tmp1 = w_z * w_z / (w0 * w0);
    auto z_tmp2 = T(1) - z_tmp1;
    auto gz = T(1) / (z_tmp2 * z_tmp2 + z_tmp1 / Q2);

    auto w_p = T(2) * std::atan(Ws);
    auto p_tmp1 = w_p * w_p / (w0 * w0);
    auto p_tmp2 = T(1) - p_tmp1;
    auto gp = T(1) / (p_tmp2 * p_tmp2 + p_tmp1 / Q2);

    auto gz_2 = gz * gz;
    auto gp_2 = gp * gp;
    auto beta = g1 - T(1);
    auto Qz = std::sqrt(g1 * g1 * (gp_2 - gz_2) / (gz_2 * (g1 + gp_2) * beta * beta));
    auto Qp = std::sqrt(g1 * (gp_2 - gz_2) / ((g1 + gz_2) * beta * beta));

    auto Ws_2 = Ws * Ws;
    auto sqrt_g1 = std::sqrt(g1);
    auto beta0 = Ws_2 + sqrt_g1 * Ws / Qz + g1;
    auto beta1 = T(2) * (Ws_2 - g1);
    auto beta2 = Ws_2 - sqrt_g1 * Ws / Qz + g1;
    auto gamma = Ws_2 + Ws / Qp + T(1);
    auto alpha1 = T(2) * (Ws_2 - T(1));
    auto alpha2 = Ws_2 - Ws / Qp + T(1);

    co[0] = beta0 / gamma;   // b0
    co[1] = beta1 / gamma;   // b1
    co[2] = beta2 / gamma;   // b2
    co[3] = -alpha1 / gamma; // -a1
    co[4] = -alpha2 / gamma; // -a2
  }

#define SOLVE_DENOM_SIMPLE                                                               \
  const auto &w0 = cutoffRadian;                                                         \
                                                                                         \
  auto q = T(0.5) / Q;                                                                   \
  auto a1 = T(-2) * std::exp(-q * w0);                                                   \
  if (q <= T(1)) {                                                                       \
    a1 *= std::cos(std::sqrt(T(1) - q * q) * w0);                                        \
  } else {                                                                               \
    a1 *= std::cosh(std::sqrt(q * q - T(1)) * w0);                                       \
  }                                                                                      \
  auto a2 = std::exp(T(-2) * q * w0);

#define SOLVE_DENOM_FULL                                                                 \
  SOLVE_DENOM_SIMPLE;                                                                    \
                                                                                         \
  auto sn = std::sin(w0 / T(2));                                                         \
  auto phi0 = T(1) - sn * sn;                                                            \
  auto phi1 = sn * sn;                                                                   \
  auto phi2 = T(4) * phi0 * phi1;                                                        \
                                                                                         \
  auto A0 = (T(1) + a1 + a2);                                                            \
  A0 *= A0;                                                                              \
  auto A1 = (T(1) - a1 + a2);                                                            \
  A1 *= A1;                                                                              \
  auto A2 = T(-4) * a2;

  void setMatchedLowpass(T cutoffRadian, T Q)
  {
    SOLVE_DENOM_FULL;

    auto sqrt_B0 = T(1) + a1 + a2;
    auto B0 = A0;

    auto R1 = Q * Q * (A0 * phi0 + A1 * phi1 + A2 * phi2);
    auto B1 = (R1 - B0 * phi0) / phi1;

    co[0] = T(0.5) * (sqrt_B0 + std::sqrt(B1)); // b0
    co[1] = sqrt_B0 - co[0];                    // b1
    co[2] = 0;                                  // b2
    co[3] = -a1;
    co[4] = -a2;
  }

  void setMatchedHighpass(T cutoffRadian, T Q)
  {
    SOLVE_DENOM_FULL;

    co[0] = Q * std::sqrt(A0 * phi0 + A1 * phi1 + A2 * phi2) / (T(4) * phi1); // b0
    co[1] = T(-2) * co[0];                                                    // b1
    co[2] = co[0];                                                            // b2
    co[3] = -a1;
    co[4] = -a2;
  }

  void setMatchedBandpass(T cutoffRadian, T Q)
  {
    SOLVE_DENOM_FULL;

    auto R1 = A0 * phi0 + A1 * phi1 + A2 * phi2;
    auto R2 = -A0 + A1 + T(4) * (phi0 - phi1) * A2;

    auto B2 = (R1 - R2 * phi1) / (T(4) * phi1 * phi1);
    auto B1 = R2 - T(4) * (phi0 - phi1) * B2;

    co[1] = T(-0.5) * std::sqrt(B1);                          // b1
    co[0] = T(0.5) * (std::sqrt(B2 + co[1] * co[1]) - co[1]); // b0
    co[2] = -co[0] - co[1];                                   // b2
    co[3] = -a1;
    co[4] = -a2;
  }

  // G: Gain.
  void setMatchedPeaking(T cutoffRadian, T Q, T G)
  {
    SOLVE_DENOM_FULL;

    auto R1 = G * G * (A0 * phi0 + A1 * phi1 + A2 * phi2);
    auto R2 = G * G * (-A0 + A1 + T(4) * (phi0 - phi1) * A2);

    const auto &B0 = A0;
    auto B2 = (R1 - R2 * phi1 - B0) / (T(4) * phi1 * phi1);
    auto B1 = R2 + B0 - T(4) * (phi0 - phi1) * B2;

    auto sqrt_B0 = T(1) + a1 + a2;
    auto sqrt_B1 = std::sqrt(B1);

    auto W = T(0.5) * (sqrt_B0 + sqrt_B1);

    co[0] = T(0.5) * (W + std::sqrt(W * W + B2)); // b0
    co[1] = T(0.5) * (sqrt_B0 - sqrt_B1);         // b1
    co[2] = -B2 / (T(4) * co[0]);                 // b2
    co[3] = -a1;
    co[4] = -a2;
  }

#undef SOLVE_DENOM_FULL

  void setSimpleMatchedLowpass(T cutoffRadian, T Q)
  {
    SOLVE_DENOM_SIMPLE;
    auto w0_2 = w0 * w0;

    auto r0 = T(1) + a1 + a2;
    auto r1 = (T(1) - a1 + a2) * w0_2
      / std::sqrt((T(1) - w0_2) * (T(1) - w0_2) + w0_2 / (Q * Q));

    co[0] = T(0.5) * (r0 + r1); // b0
    co[1] = r0 - co[0];         // b1
    co[2] = 0;                  // b2
    co[3] = -a1;
    co[4] = -a2;
  }

  void setSimpleMatchedHighpass(T cutoffRadian, T Q)
  {
    SOLVE_DENOM_SIMPLE;
    auto w0_2 = w0 * w0;

    auto r1
      = (T(1) - a1 + a2) / std::sqrt((T(1) - w0_2) * (T(1) - w0_2) + w0_2 / (Q * Q));

    co[0] = T(0.25) * r1;  // b0
    co[1] = T(-2) * co[0]; // b1
    co[2] = co[0];         // b2
    co[3] = -a1;
    co[4] = -a2;
  }

  void setSimpleMatchedBandpass(T cutoffRadian, T Q)
  {
    SOLVE_DENOM_SIMPLE;
    auto w0_2 = w0 * w0;

    auto r0 = (T(1) + a1 + a2) / (w0 * Q);
    auto r1 = (T(1) - a1 + a2) * w0
      / (Q * std::sqrt((T(1) - w0_2) * (T(1) - w0_2) + w0_2 / (Q * Q)));

    co[0] = T(0.5) * r0 + T(0.25) * r1; // b0
    co[1] = T(-0.5) * r1;               // b1
    co[2] = -co[0] - co[1];             // b2
    co[3] = -a1;
    co[4] = -a2;
  }

#undef SOLVE_DENOM_SIMPLE

  void setMatchedShelvingOnePole(T cutoffNormalized, T G)
  {
    const auto &fc = cutoffNormalized;
    // constexpr T fm = T(0.9);
    constexpr T phim = T(1.9510565162951536); // 1 - cos(pi * fm).

    constexpr T pp = T(2) / (pi * pi);
    constexpr T xi = pp / (phim * phim) - T(1) / phim;

    auto alpha = xi + pp / (G * fc * fc);
    auto beta = xi + pp * G / (fc * fc);

    auto neg_a1 = alpha / (T(1) + alpha + std::sqrt(T(1) + T(2) * alpha));
    auto b = -beta / (T(1) + beta + std::sqrt(T(1) + T(2) * beta));
    co[0] = (T(1) - neg_a1) / (T(1) + b); // b0
    co[1] = b * co[0];                    // b1
    co[2] = 0;                            // b2
    co[3] = neg_a1;                       // -a1
    co[4] = 0;                            // -a2
  }

  void setCutoff(MatchedBiquadType type, T normalizedFreq, T Q, T gainDecibel)
  {
    auto cutoffRadian = T(2) * pi * normalizedFreq;
    auto gainAmp = std::pow(T(10), gainDecibel / T(20));

    if (type == MatchedBiquadType::orfanidisPeaking) {
      setOrfanidisPeaking(cutoffRadian, Q, gainAmp);
    } else if (type == MatchedBiquadType::massbergLowpass) {
      setMassbergLowpass(cutoffRadian, Q);
    } else if (type == MatchedBiquadType::matchedLowpass) {
      setMatchedLowpass(cutoffRadian, Q);
    } else if (type == MatchedBiquadType::matchedHighpass) {
      setMatchedHighpass(cutoffRadian, Q);
    } else if (type == MatchedBiquadType::matchedBandpass) {
      setMatchedBandpass(cutoffRadian, Q);
    } else if (type == MatchedBiquadType::matchedPeaking) {
      setMatchedPeaking(cutoffRadian, Q, gainAmp);
    } else if (type == MatchedBiquadType::simpleMatchedLowpass) {
      setSimpleMatchedLowpass(cutoffRadian, Q);
    } else if (type == MatchedBiquadType::simpleMatchedHighpass) {
      setSimpleMatchedHighpass(cutoffRadian, Q);
    } else if (type == MatchedBiquadType::simpleMatchedBandpass) {
      setSimpleMatchedBandpass(cutoffRadian, Q);
    } else if (type == MatchedBiquadType::matchedShelvingOnePole) {
      setMatchedShelvingOnePole(T(2) * normalizedFreq, gainAmp);
    }
  }

  void reset(T value = 0)
  {
    x1 = value;
    x2 = value;
    y1 = value;
    y2 = value;
  }

  T process(T x0)
  {
    auto y0 = co[0] * x0 + co[1] * x1 + co[2] * x2 + co[3] * y1 + co[4] * y2;

    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;

    return y0;
  }
};

template<typename Float>
std::string vectorToJsonList(const std::string &name, const std::vector<Float> &data)
{
  std::string text = std::format("\"{}\":[", name);
  for (const auto &value : data) text += std::format("{},", value);
  text.pop_back();
  text += "],";
  return text;
}

template<typename Float>
auto getImpulseResponse(
  MatchedBiquadType filterType,
  Float sampleRate,
  Float cutoffNormalized,
  Float Q,
  Float gainDecibel)
{
  MatchedBiquad<Float> flt;
  flt.setCutoff(filterType, cutoffNormalized, Q, gainDecibel);

  std::vector<Float> sig;
  sig.resize(size_t(sampleRate));
  sig[0] = flt.process(Float(1));
  for (size_t i = 1; i < sig.size(); ++i) sig[i] = flt.process(0);
  for (auto &v : sig) v = std::isfinite(v) ? v : 0;
  return vectorToJsonList<Float>(
    matchedBiquadNames[static_cast<size_t>(filterType)], sig);
}

int main()
{
  using Float = double;
  constexpr Float sampleRate = Float(48000);
  Float cutoffNormalized = Float(1000) / sampleRate; // [rad/2pi].
  Float Q = std::numbers::sqrt2_v<Float> / Float(2);
  Float gainDecibel = Float(20);

  auto text = std::format(
    "{{\"cutoffNormalized\":{},\"Q\":{},\"gainDecibel\":{},\"impulseResponse\":{{",
    cutoffNormalized, Q, gainDecibel);

  auto getIR = std::bind(
    getImpulseResponse<Float>, std::placeholders::_1, sampleRate, cutoffNormalized, Q,
    gainDecibel);

  text += getIR(MatchedBiquadType::orfanidisPeaking);
  text += getIR(MatchedBiquadType::massbergLowpass);
  text += getIR(MatchedBiquadType::matchedLowpass);
  text += getIR(MatchedBiquadType::matchedHighpass);
  text += getIR(MatchedBiquadType::matchedBandpass);
  text += getIR(MatchedBiquadType::matchedPeaking);
  text += getIR(MatchedBiquadType::simpleMatchedLowpass);
  text += getIR(MatchedBiquadType::simpleMatchedHighpass);
  text += getIR(MatchedBiquadType::simpleMatchedBandpass);
  text += getIR(MatchedBiquadType::matchedShelvingOnePole);

  text.pop_back();
  text += "}}";
  std::ofstream ofs("impulse_response.json");
  ofs << text;

  return 0;
}
