/*
Reference:
- [Design IIR Filters Using Cascaded Biquads - Neil
Robertson](https://www.dsprelated.com/showarticle/1137.php)
- [Design IIR Butterworth Filters Using 12 Lines of Code - Neil
Robertson](https://www.dsprelated.com/showarticle/1119.php)
*/

#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <sndfile.h>
#include <vector>

constexpr double pi = 3.14159265358979323846264338;
constexpr double twopi = 2.0 * pi;

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

template<typename Sample> struct Bessel4 {
  // The poles of Bessel filter becomes conjugate pair when degree is even.
  //
  // For example if degree is 6, then poles are: [p0, p1, p2, p2*, p1*, p0*].
  // * is complex conjugate operator.
  //
  // Use double. float may be sufficient if delay time is short (delayInSamples <~ 1000).
  //
  // ```python
  // import numpy
  // import scipy.signal as signal
  //
  // # z is zero, p is pole, k is gain.
  // z, p, k = signal.besselap(order, norm="delay")
  // ```
  constexpr static uint8_t halfDegree = 2; // 4 の半分。

  // 複素共役の極は無くても計算できる。
  constexpr static std::array<std::complex<Sample>, halfDegree> analogPole{{
    {Sample(-2.1037893971796273), Sample(+2.6574180418567526)},
    {Sample(-2.8962106028203722), Sample(+0.8672341289345038)},
  }};

  std::array<Sample, 2> x0{};
  std::array<Sample, 2> x1{};
  std::array<Sample, 2> x2{};
  std::array<Sample, 2> y0{};
  std::array<Sample, 2> y1{};
  std::array<Sample, 2> y2{};
  std::array<std::array<Sample, 2>, halfDegree> co; // フィルタ係数 a1 と a2 。
  Sample gain = 1;

  Bessel4()
  {
    for (auto &coef : co) coef.fill(0);
  }

  void reset(Sample value = 0)
  {
    x0.fill(value);
    x1.fill(value);
    x2.fill(value);
    y0.fill(value);
    y1.fill(value);
    y2.fill(value);
  }

  void setDelay(Sample delayInSamples)
  {
    auto wo = Sample(2) / delayInSamples; // 遅延サンプル数を周波数に変換。

    gain = Sample(1);
    for (uint8_t i = 0; i < co.size(); ++i) {
      std::complex<Sample> pole = wo * analogPole[i]; // カットオフ周波数の適用。
      pole = (Sample(2) + pole) / (Sample(2) - pole); // バイリニア変換。
      co[i][0] = Sample(-2) * pole.real();
      co[i][1] = std::norm(pole);
      gain *= (Sample(1) + co[i][0] + co[i][1]) / Sample(4);
    }
  }

  Sample process(Sample input)
  {
    x0[0] = input;

    for (uint8_t i = 0; i < halfDegree; ++i) {
      y0[i] = x0[i] + Sample(2) * x1[i] + x2[i] - co[i][0] * y1[i] - co[i][1] * y2[i];

      x2[i] = x1[i];
      x1[i] = x0[i];
      y2[i] = y1[i];
      y1[i] = y0[i];

      if (i + 1 < halfDegree) x0[i + 1] = y0[i];
    }

    return gain * y0[1];
  }
};

int main()
{
  size_t sampleRate = 48000;
  size_t delay = 1024;

  std::vector<float> wav;
  wav.resize(2 * delay);

  Bessel4<float> bessel;
  bessel.setDelay(float(delay));

  for (size_t i = 0; i < wav.size(); ++i) wav[i] = bessel.process(1.0f);

  writeWave("snd/bessel.wav", wav, sampleRate);

  return 0;
}
