#include "lib/LambertW.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sndfile.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

int32_t
writeWave(const char *filename, std::vector<float> &buffer, const size_t &samplerate)
{
  SF_INFO sfinfo;
  memset(&sfinfo, 0, sizeof(sfinfo));
  sfinfo.samplerate = samplerate;
  sfinfo.frames = buffer.size();
  sfinfo.channels = 1;
  sfinfo.format = (SF_FORMAT_WAV | SF_FORMAT_FLOAT);

  SNDFILE *file = sf_open(filename, SFM_WRITE, &sfinfo);
  if (!file) {
    std::cout << "Error: sf_open failed." << std::endl;
    return 1;
  }

  size_t length = sfinfo.channels * buffer.size();
  if (sf_write_float(file, &buffer[0], length) != length)
    std::cout << sf_strerror(file) << std::endl;

  sf_close(file);

  return 0;
}

template<typename Sample> class ExpAD {
private:
  Sample gain = 0;
  Sample valueA = 0;
  Sample alphaA = 0;
  Sample valueD = 0;
  Sample alphaD = 0;

public:
  bool isTerminated() { return valueD <= std::numeric_limits<Sample>::epsilon(); }

  void reset1(Sample sampleRate, Sample attackSeconds, Sample decaySeconds)
  {
    constexpr Sample epsilon = Sample(1e-5);

    valueA = Sample(1);
    alphaA
      = std::pow(epsilon, Sample(1) / std::max(Sample(1), attackSeconds * sampleRate));

    valueD = Sample(1);
    alphaD
      = std::pow(epsilon, Sample(1) / std::max(Sample(1), decaySeconds * sampleRate));

    const auto log_a = std::log(alphaA);
    const auto log_d = std::log(alphaD);
    const auto t_p = std::log(log_d / (log_a + log_d)) / log_a;
    gain = Sample(1) / ((Sample(1) - std::pow(alphaA, t_p)) * std::pow(alphaD, t_p));
  }

  void reset2(Sample sampleRate, Sample attackSeconds, Sample decaySeconds)
  {
    constexpr Sample epsilon = Sample(1e-5);

    const auto a_ = std::log(epsilon) / attackSeconds;
    const auto d_ = std::log(epsilon) / decaySeconds;

    valueA = Sample(1);
    alphaA = std::exp(a_ / sampleRate);

    valueD = Sample(1);
    alphaD = std::exp(d_ / sampleRate);

    const auto t_p = -std::log1p(a_ / d_) / a_;
    gain = Sample(1) / ((Sample(1) - std::exp(a_ * t_p)) * std::exp(d_ * t_p));
  }

  void reset3(Sample sampleRate, Sample peakSeconds, Sample releaseSeconds)
  {
    constexpr Sample epsilon = std::numeric_limits<Sample>::epsilon();

    const auto decaySeconds = releaseSeconds - std::log(epsilon) * peakSeconds;
    const auto d_ = std::log(epsilon) / decaySeconds;
    const auto x_ = d_ * peakSeconds;
    const auto a_ = Sample(utl::LambertW(-1, x_ * std::exp(x_))) / peakSeconds - d_;

    const auto attackSeconds = -std::log(epsilon) / std::log(-a_);
    valueA = Sample(1);
    alphaA = std::exp(a_ / sampleRate);

    valueD = Sample(1);
    alphaD = std::exp(d_ / sampleRate);

    gain = Sample(1)
      / ((Sample(1) - std::exp(a_ * peakSeconds)) * std::exp(d_ * peakSeconds));
  }

  Sample process()
  {
    valueA *= alphaA;
    valueD *= alphaD;
    return gain * (Sample(1) - valueA) * valueD;
  }
};

int main()
{
  constexpr float sampleRate = 48000;
  std::vector<float> wav(sampleRate);
  ExpAD<double> envelope;

  envelope.reset1(sampleRate, float(1), float(2));
  for (auto &buf : wav) buf = envelope.process();
  writeWave("snd/ExpAD1.wav", wav, sampleRate);

  envelope.reset2(sampleRate, float(1), float(2));
  for (auto &buf : wav) buf = envelope.process();
  writeWave("snd/ExpAD2.wav", wav, sampleRate);

  envelope.reset3(sampleRate, float(0.1), float(2));
  for (auto &buf : wav) buf = envelope.process();
  writeWave("snd/ExpAD3.wav", wav, sampleRate);

  return 0;
}
