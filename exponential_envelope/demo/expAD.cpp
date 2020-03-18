#include <cmath>
#include <iostream>
#include <sndfile.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

constexpr float sampleRate = 44100;

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

class ExpAD {
public:
  void setup(float sampleRate) { this->sampleRate = sampleRate; }
  bool isTerminated() { return valueD <= threshold; }

  // attack and decay in seconds.
  void reset(float attack, float decay)
  {
    valueA = 1.0f;
    if (attack < 1e-5) attack = 1e-5;
    alphaA = powf(threshold, 1.0f / (attack * sampleRate));

    valueD = 1.0f;
    if (decay < 1e-5) decay = 1e-5;
    alphaD = powf(threshold, 1.0f / (decay * sampleRate));

    if (attack <= 0.0f) {
      gain = 1.0f;
    } else if (decay <= 0.0f) {
      gain = 0.0f;
    } else {
      auto log_a = logf(alphaA);
      auto log_d = logf(alphaD);
      auto t_p = logf(log_d / (log_a + log_d)) / log_a;
      gain = 1.0f / ((1.0f - powf(alphaA, t_p)) * powf(alphaD, t_p));
    }
  }

  float process()
  {
    valueA *= alphaA;
    valueD *= alphaD;
    return gain * (1.0f - threshold - valueA) * (valueD - threshold);
  }

protected:
  const float threshold = 1e-5;
  float sampleRate = 44100;
  float gain = 0;
  float valueA = 0;
  float alphaA = 0;
  float valueD = 0;
  float alphaD = 0;
};

int main()
{
  std::vector<float> wav(sampleRate);

  ExpAD envelope;
  envelope.setup(sampleRate);
  envelope.reset(1.0f, 2.0f);

  for (auto &buf : wav) buf = envelope.process();

  writeWave("snd/ExpAD.wav", wav, sampleRate);

  return 0;
}
