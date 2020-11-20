#include <algorithm>
#include <array>
#include <iostream>

template<typename Sample> struct SOCPFIR {
  constexpr static size_t bufferSize = 7;
  constexpr static size_t intDelay = 3;

  constexpr static std::array<std::array<Sample, bufferSize>, 4> coefficient{{
    {Sample(0.03396642725330925), Sample(-0.12673821137646601),
     Sample(0.5759982312324312), Sample(0.6592123095604063), Sample(-0.19435321143573606),
     Sample(0.0782612693103079), Sample(-0.025807862651826587)},
    {Sample(0.021616078095824397), Sample(-0.07539816970638001),
     Sample(0.2653441329619578), Sample(0.9081714824861011), Sample(-0.16017585860369898),
     Sample(0.059489586593950955), Sample(-0.018863293456169244)},
    {Sample(-0.018863293456169286), Sample(0.05948958659395098),
     Sample(-0.16017585860369907), Sample(0.908171482486101), Sample(0.2653441329619578),
     Sample(-0.07539816970638011), Sample(0.02161607809582444)},
    {Sample(-0.02580786265182662), Sample(0.07826126931030812),
     Sample(-0.1943532114357363), Sample(0.6592123095604064), Sample(0.5759982312324308),
     Sample(-0.12673821137646582), Sample(0.033966427253309124)},
  }};
};

template<typename Sample, typename FractionalDelayFIR> class TruePeakMeterFIR {
  std::array<Sample, FractionalDelayFIR::bufferSize> buf{};

public:
  void reset() { buf.fill(Sample(0)); }

  Sample process(Sample input)
  {
    for (size_t i = 0; i < buf.size() - 1; ++i) buf[i] = buf[i + 1];
    buf.back() = input;

    Sample max = std::fabs(buf[FractionalDelayFIR::intDelay]);
    for (const auto &phase : FractionalDelayFIR::coefficient) {
      Sample sum = 0;
      for (size_t i = 0; i < phase.size(); ++i) sum += buf[i] * phase[i];
      max = std::max(max, std::fabs(sum));
    }
    return max;
  }
};

int main()
{
  // 1 サンプルだけ計算する例。
  TruePeakMeterFIR<float, SOCPFIR<float>> bs1770;
  std::cout << bs1770.process(1.0f) << "\n";

  return 0;
}
