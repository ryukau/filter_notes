#include <cfloat>
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

class ExpPolyEnvelope {
public:
  // attack の単位は秒。
  // curve は任意の値。 β に相当。
  // attack と curve が大きいと計算結果が inf になるときがあるので注意。
  void reset(double sampleRate, double attack, double curve)
  {
    alpha = attack * curve;

    auto betaMin = alpha / getTerminationTime();
    if (curve < betaMin) curve = betaMin;

    peak = pow(alpha / curve, alpha) * exp(-alpha);
    gamma = exp(-curve / sampleRate);
    tick = 1.0 / sampleRate;

    time = 0.0;
    value = 1.0;
  }

  bool isReleasing() { return time >= attack; }
  double getTerminationTime() { return pow(DBL_MAX, 1.0 / alpha); }

  double process()
  {
    auto output = pow(time, alpha) * value / peak;
    if (!std::isfinite(output)) return 0.0; // 念のため。
    time += tick;
    value *= gamma;
    return output;
  }

protected:
  double value = 0;
  double peak = 1;
  double gamma = 0;
  double attack = 0;
  double tick = 0;
  double alpha = 0;
  double time = 0;
};

int main()
{
  std::vector<float> wav(4 * size_t(sampleRate));

  ExpPolyEnvelope envelope;

  envelope.reset(sampleRate, 1.0f, 4.0f);
  for (auto &buf : wav) buf = envelope.process();
  writeWave("snd/ExpPoly.wav", wav, sampleRate);

  return 0;
}
