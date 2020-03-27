#include <cmath>
#include <cstring>
#include <iostream>
#include <random>
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
    return EXIT_FAILURE;
  }

  size_t length = sfinfo.channels * buffer.size();
  if (sf_write_float(file, &buffer[0], length) != length)
    std::cout << sf_strerror(file) << std::endl;

  sf_close(file);

  return 0;
}

constexpr float sampleRate = 48000.0f;
constexpr size_t nFrame = 512;
constexpr size_t nBuffer = 32;

#include <iostream>

template<typename Sample> class RateLimiter {
public:
  void setup(Sample sampleRate) { dt = Sample(1) / sampleRate; }
  void reset() { y1 = 0; }

  void setRate(Sample risingSlewRate, Sample fallingSlewRate)
  {
    rise = risingSlewRate;
    fall = fallingSlewRate;
  }

  // Time to reach 1 from 0. risingTime and fallingTime should be positive.
  void setTime(Sample sampleRate, Sample risingTime, Sample fallingTime)
  {
    rise = Sample(1) / (sampleRate * risingTime);
    fall = -Sample(1) / (sampleRate * fallingTime);
  }

  Sample process(Sample input)
  {
    auto rate = (input - y1) / dt;
    if (rate > rise)
      y1 += dt * rise;
    else if (rate < fall)
      y1 += dt * fall;
    else
      y1 = input;
    return y1;
  }

private:
  Sample dt = 1.0 / 44100.0;
  Sample y1 = 0;
  Sample rise = 10000;  // 正の値であるべき。
  Sample fall = -10000; // 負の値であるべき。
};

struct DSP {
  float gain = 1.0f;
  RateLimiter<float> slewGain;

  void process(const size_t frame, float *out)
  {
    for (size_t i = 0; i < frame; ++i) out[i] = slewGain.process(gain);
  }
};

int main()
{
  std::vector<float> wav(nFrame * nBuffer);

  float out[nFrame];

  std::minstd_rand rng{0};
  std::uniform_real_distribution<float> gainDist(0.0f, 1.0f);
  DSP dsp;
  dsp.slewGain.setTime(sampleRate, 8e-7f, 4e-7f);

  float gainValue = 1.0f;
  for (size_t idx = 0; idx < nBuffer; ++idx) {
    if (idx % 4 == 0) dsp.gain = gainDist(rng);
    dsp.process(nFrame, out);
    std::memcpy(&wav[idx * nFrame], out, sizeof(float) * nFrame);
  }

  return writeWave("snd/rate_limiter.wav", wav, sampleRate);
}
