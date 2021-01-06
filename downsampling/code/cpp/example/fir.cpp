#include <array>
#include <iostream>
#include <numeric>
#include <vector>

template<typename Sample> struct FirRemez64Decimation2 {
  constexpr static size_t bufferSize = 32;
  constexpr static size_t nPhase = 2;

  constexpr static std::array<std::array<Sample, bufferSize>, nPhase> coefficient{{
    {Sample(-0.00012586358900962356), Sample(0.0001828145651312807),
     Sample(8.205416319885108e-05),   Sample(-0.0007239553549487794),
     Sample(0.0018923639640423555),   Sample(-0.00353192867622875),
     Sample(0.005331459582347676),    Sample(-0.006667992332855854),
     Sample(0.006636334343333397),    Sample(-0.004162169091786725),
     Sample(-0.0018313041529962014),  Sample(0.012291623826628002),
     Sample(-0.02812500771752753),    Sample(0.05102267479538475),
     Sample(-0.08758235952225661),    Sample(0.18605559261058696),
     Sample(0.4035431721330982),      Sample(-0.03627960335535102),
     Sample(-0.006631516329342658),   Sample(0.02064755543805059),
     Sample(-0.023817278019366648),   Sample(0.021451149873294335),
     Sample(-0.0164261830615856),     Sample(0.010706859341764241),
     Sample(-0.005602179907252967),   Sample(0.001817817136086864),
     Sample(0.00046580645271121403),  Sample(-0.0014510531292568874),
     Sample(0.0015411334474003102),   Sample(-0.0011736418939169283),
     Sample(0.0006849141024872853),   Sample(-0.00039292037784520776)},
    {Sample(-0.00039292037784520776), Sample(0.0006849141024872853),
     Sample(-0.0011736418939169283),  Sample(0.0015411334474003102),
     Sample(-0.0014510531292568874),  Sample(0.00046580645271121403),
     Sample(0.001817817136086864),    Sample(-0.005602179907252967),
     Sample(0.010706859341764241),    Sample(-0.0164261830615856),
     Sample(0.021451149873294335),    Sample(-0.023817278019366648),
     Sample(0.02064755543805059),     Sample(-0.006631516329342658),
     Sample(-0.03627960335535102),    Sample(0.4035431721330982),
     Sample(0.18605559261058696),     Sample(-0.08758235952225661),
     Sample(0.05102267479538475),     Sample(-0.02812500771752753),
     Sample(0.012291623826628002),    Sample(-0.0018313041529962014),
     Sample(-0.004162169091786725),   Sample(0.006636334343333397),
     Sample(-0.006667992332855854),   Sample(0.005331459582347676),
     Sample(-0.00353192867622875),    Sample(0.0018923639640423555),
     Sample(-0.0007239553549487794),  Sample(8.205416319885108e-05),
     Sample(0.0001828145651312807),   Sample(-0.00012586358900962356)},
  }};
};

template<typename Sample, typename FIR> class FirFilter {
  std::array<std::array<Sample, FIR::bufferSize>, FIR::nPhase> buf;

public:
  FirFilter() { reset(); }

  void reset()
  {
    for (auto &bf : buf) bf.fill(Sample(0));
  }

  Sample process(std::array<Sample, FIR::nPhase> &input)
  {
    Sample sum = 0;

    for (size_t ph = 0; ph < FIR::nPhase; ++ph) {
      auto &bf = buf[ph];
      for (size_t i = 0; i < bf.size() - 1; ++i) bf[i] = bf[i + 1];
      bf.back() = input[ph];

      const auto &fir = FIR::coefficient[ph];
      for (size_t i = 0; i < fir.size(); ++i) sum += bf[i] * fir[i];
    }

    return sum;
  }
};

int main()
{
  FirFilter<float, FirRemez64Decimation2<float>> filter;

  std::vector<float> data(64);
  std::iota(data.begin(), data.end(), 0.0f);

  std::array<float, 2> input;
  for (size_t i = 0; i < data.size(); i += 2) {
    input[0] = data[i];
    input[1] = data[i + 1];
    std::cout << filter.process(input) << "\n";
  }

  return 0;
}
