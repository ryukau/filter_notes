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

struct DSP {
  float gain = 1.0f;

  void process(const size_t frame, float *out)
  {
    for (size_t i = 0; i < frame; ++i) out[i] = gain;
  }
};

int main()
{
  std::vector<float> wav(nFrame * nBuffer);

  float out[nFrame];

  std::minstd_rand rng{0};
  std::uniform_real_distribution<float> gainDist(0.0f, 1.0f);
  DSP dsp;

  float gainValue = 1.0f;
  for (size_t idx = 0; idx < nBuffer; ++idx) {
    if (idx % 4 == 0) gainValue = gainDist(rng);
    dsp.gain = gainValue;
    dsp.process(nFrame, out);
    std::memcpy(&wav[idx * nFrame], out, sizeof(float) * nFrame);
  }

  return writeWave("snd/naive.wav", wav, sampleRate);
}
