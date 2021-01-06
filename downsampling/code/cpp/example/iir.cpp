#include <array>
#include <iostream>
#include <numeric>
#include <vector>

/*
---python
import numpy
from scipy import signal
sos = signal.ellip(12, 0.01, 100, 0.4, "low", output="sos", fs=8)
---
*/
template<typename Sample> struct SosEllipticDecimation8 {
  constexpr static size_t nSection = 6;

  constexpr static std::array<std::array<Sample, 5>, nSection> co{{
    {Sample(2.0478175697515062e-05), Sample(5.313899743869536e-06),
     Sample(2.047817569751506e-05), Sample(-1.7517115385026274),
     Sample(0.7700392502488301)},
    {Sample(1.0), Sample(-1.438980096034124), Sample(1.0), Sample(-1.774966834164871),
     Sample(0.8105614413580351)},
    {Sample(1.0), Sample(-1.7274028784750513), Sample(1.0), Sample(-1.808053369266683),
     Sample(0.867640704095363)},
    {Sample(1.0), Sample(-1.8120823742998813), Sample(1.0000000000000004),
     Sample(-1.8388054096377313), Sample(0.9192163691033324)},
    {Sample(1.0), Sample(-1.8441397779182482), Sample(1.0000000000000002),
     Sample(-1.8637392843186766), Sample(0.9580528037069385)},
    {Sample(1.0), Sample(-1.855844538281174), Sample(1.0), Sample(-1.8856454058448944),
     Sample(0.9871065355804314)},
  }};
};

// SOS: Second order sections.
template<typename Sample, typename IIR, size_t oversample = 8> class SosFilter {
public:
  void reset()
  {
    x0.fill(0);
    x1.fill(0);
    x2.fill(0);
    y0.fill(0);
    y1.fill(0);
    y2.fill(0);
  }

  void push(Sample input)
  {
    x0[0] = input;
    for (size_t i = 1; i < IIR::nSection; ++i) x0[i] = y0[i - 1];

    for (size_t i = 0; i < IIR::nSection; ++i) {
      y0[i] = +IIR::co[i][0] * x0[i] + IIR::co[i][1] * x1[i] + IIR::co[i][2] * x2[i]
        - IIR::co[i][3] * y1[i] - IIR::co[i][4] * y2[i];
    }

    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }

  inline Sample output() { return y0[IIR::nSection - 1]; }

  Sample process(const std::array<Sample, oversample> &input)
  {
    for (const auto &value : input) push(value);
    return output();
  }

  std::array<Sample, IIR::nSection> x0{};
  std::array<Sample, IIR::nSection> x1{};
  std::array<Sample, IIR::nSection> x2{};
  std::array<Sample, IIR::nSection> y0{};
  std::array<Sample, IIR::nSection> y1{};
  std::array<Sample, IIR::nSection> y2{};
};

int main()
{
  SosFilter<float, SosEllipticDecimation8<float>> filter;

  std::vector<float> data(64);
  std::iota(data.begin(), data.end(), 0.0f);

  std::cout << "\n--- Use `process` to process each n samples at once.\n";
  std::array<float, 8> input;
  for (size_t i = 0; i < data.size(); i += 8) {
    for (size_t j = 0; j < input.size(); ++j) input[j] = data[i + j];
    std::cout << filter.process(input) << "\n";
  }

  std::cout << "\n--- Use `push` and `output` to get output for each sample.\n";
  filter.reset();
  for (size_t i = 0; i < data.size(); ++i) {
    filter.push(data[i]);
    if (i % 8 == 7) std::cout << filter.output() << "\n";
  }

  return 0;
}
