#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <sndfile.h>
#include <vector>

constexpr double halfpi = 1.57079632679489661923;
constexpr double pi = 3.14159265358979323846;
constexpr double twopi = 6.28318530717958647692;

void writeWave(std::string filename, std::vector<float> &buffer, const size_t &samplerate)
{
  SF_INFO sfinfo;
  memset(&sfinfo, 0, sizeof(sfinfo));
  sfinfo.samplerate = int(samplerate);
  sfinfo.frames = buffer.size();
  sfinfo.channels = 1;
  sfinfo.format = (SF_FORMAT_WAV | SF_FORMAT_FLOAT);

  SNDFILE *file = sf_open(filename.c_str(), SFM_WRITE, &sfinfo);
  if (!file) {
    std::cout << "Error: sf_open failed." << std::endl;
    exit(EXIT_FAILURE);
  }

  size_t length = sfinfo.channels * buffer.size();
  if (sf_write_float(file, &buffer[0], length) != length)
    std::cout << sf_strerror(file) << std::endl;

  if (sf_close(file) != 0) {
    std::cout << "Error: sf_close failed." << std::endl;
    exit(EXIT_FAILURE);
  }
}

template<typename Sample> inline Sample frequencyToMidinote(Sample freq)
{
  return Sample(12) * std::log2(freq / Sample(440)) + Sample(69);
}

template<typename Sample> inline Sample midinoteToFrequency(Sample note)
{
  return Sample(440) * std::exp2((note - Sample(69)) / Sample(12));
}

template<typename T> struct HalfBandCoefficient {
  static constexpr std::array<T, 9> h0_a{
    T(0.0765690656031399), T(0.264282270318935),  T(0.47939467893641907),
    T(0.661681722389424),  T(0.7924031566294969), T(0.8776927911111817),
    T(0.9308500986629166), T(0.9640156636878193), T(0.9862978287283355),
  };
  static constexpr std::array<T, 10> h1_a{
    T(0.019911761024506557), T(0.16170648261075027), T(0.37320978687920564),
    T(0.5766558985008232),   T(0.7334355636406803),  T(0.8399227128761151),
    T(0.9074601780285125),   T(0.9492937701934973),  T(0.9760539731706528),
    T(0.9955323321150525),
  };
};

template<typename T> struct HalfBandCoefficientHiir {
  static constexpr std::array<T, 9> h0_a{
    T(0.0765690656031399), T(0.264282270318935),  T(0.4793946789364191),
    T(0.661681722389424),  T(0.792403156629497),  T(0.8776927911111816),
    T(0.9308500986629166), T(0.9640156636878193), T(0.9862978287283355),
  };
  static constexpr std::array<T, 10> h1_a{
    T(0.019911761024506557), T(0.16170648261075027), T(0.37320978687920564),
    T(0.5766558985008232),   T(0.7334355636406803),  T(0.8399227128761151),
    T(0.9074601780285125),   T(0.9492937701934973),  T(0.9760539731706528),
    T(0.9955323321150525),
  };
};

template<typename Sample> class FirstOrderAllpass {
private:
  Sample x1 = 0;
  Sample y1 = 0;

public:
  void reset()
  {
    x1 = 0;
    y1 = 0;
  }

  Sample process(Sample x0, Sample a)
  {
    y1 = a * (x0 - y1) + x1;
    x1 = x0;
    return y1;
  }
};

template<typename Sample, typename Coefficient> class HalfBandIIRSplit {
private:
  std::array<FirstOrderAllpass<Sample>, Coefficient::h0_a.size()> ap0;
  std::array<FirstOrderAllpass<Sample>, Coefficient::h1_a.size()> ap1;

public:
  void reset()
  {
    for (auto &ap : ap0) ap.reset();
    for (auto &ap : ap1) ap.reset();
  }

  // input[0] must be earlier sample.
  Sample process(std::array<Sample, 2> &input)
  {
    auto s0 = input[0];
    for (size_t i = 0; i < ap0.size(); ++i) s0 = ap0[i].process(s0, Coefficient::h0_a[i]);
    auto s1 = input[1];
    for (size_t i = 0; i < ap1.size(); ++i) s1 = ap1[i].process(s1, Coefficient::h1_a[i]);
    return Sample(0.5) * (s0 + s1);
  }
};

template<typename Sample, size_t nSection> class FirstOrderAllpassSections {
private:
  std::array<Sample, nSection> x{};
  std::array<Sample, nSection> y{};

public:
  void reset()
  {
    x.fill(0);
    y.fill(0);
  }

  Sample process(Sample input, const std::array<Sample, nSection> &a)
  {
    for (size_t i = 0; i < nSection; ++i) {
      y[i] = a[i] * (input - y[i]) + x[i];
      x[i] = input;
      input = y[i];
    }
    return y.back();
  }
};

template<typename Sample, typename Coefficient> class HalfBandIIRArray {
private:
  FirstOrderAllpassSections<Sample, Coefficient::h0_a.size()> ap0;
  FirstOrderAllpassSections<Sample, Coefficient::h1_a.size()> ap1;

public:
  void reset()
  {
    ap0.reset();
    ap1.reset();
  }

  // input[0] must be earlier sample.
  Sample process(std::array<Sample, 2> &input)
  {
    auto s0 = ap0.process(input[0], Coefficient::h0_a);
    auto s1 = ap1.process(input[1], Coefficient::h1_a);
    return Sample(0.5) * (s0 + s1);
  }
};

template<typename Sample> class TestOsc {
public:
  Sample phase = 0;
  Sample delta = 0;

  void reset() { phase = 0; }

  Sample process(Sample sampleRate, Sample frequency)
  {
    phase += frequency / sampleRate;
    phase -= std::floor(phase);
    return std::sin(Sample(twopi) * phase);
  }
};

constexpr size_t nLoop = 1;

template<typename Filter>
void bench(float sampleRate, std::string name, std::vector<float> &src)
{
  Filter halfbandiir;

  std::vector<float> wav;
  wav.resize(src.size() / 2);

  double sumElapsed = 0.0;
  std::array<float, 2> phases{};
  for (size_t loop = 0; loop < nLoop; ++loop) {
    halfbandiir.reset();
    for (size_t i = 0; i < src.size(); i += 2) {
      phases[0] = src[i];
      phases[1] = src[i + 1];

      auto start = std::chrono::steady_clock::now();
      wav[i / 2] = halfbandiir.process(phases);
      auto end = std::chrono::steady_clock::now();

      std::chrono::duration<double, std::milli> elapsed = end - start;
      sumElapsed += elapsed.count();
    }
  }
  std::cout << name << ": " << wav.size() << "[sample], " << sumElapsed << "[ms], "
            << "\n";

  writeWave("snd/" + name + ".wav", wav, size_t(sampleRate));
}

std::vector<float> generate2xSourceSignal(float sampleRate)
{
  auto upRate = float(2) * sampleRate;

  std::vector<float> wav;
  wav.resize(size_t(upRate));

  TestOsc<float> osc;
  osc.reset();

  auto freqSlope = sampleRate / upRate;
  for (size_t i = 0; i < wav.size(); ++i) {
    wav[i] = osc.process(upRate, float(i) * freqSlope);
  }

  writeWave("snd/Source.wav", wav, size_t(upRate));
  return wav;
}

int main()
{
  constexpr float sampleRate = 48000.0f;
  auto source = generate2xSourceSignal(sampleRate);

  std::cout << "--- Warm up\n";
  bench<HalfBandIIRArray<float, HalfBandCoefficient<float>>>(sampleRate, "Array", source);

  std::cout << "\n--- Benchmark\n";
  bench<HalfBandIIRSplit<float, HalfBandCoefficient<float>>>(sampleRate, "Split", source);
  bench<HalfBandIIRArray<float, HalfBandCoefficient<float>>>(sampleRate, "Array", source);
  bench<HalfBandIIRSplit<float, HalfBandCoefficientHiir<float>>>(
    sampleRate, "SplitHiir", source);
  bench<HalfBandIIRArray<float, HalfBandCoefficientHiir<float>>>(
    sampleRate, "ArrayHiir", source);
  return 0;
}
