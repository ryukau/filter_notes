/*
$ g++ -lsndfile -O3 bench_smoother.cpp ; ./a.out ; python plot.py
index # dry run
Total[ms]49.287917
Average[ms]0.048133

ramp
Total[ms]51.467152
Average[ms]0.050261

index
Total[ms]49.488216
Average[ms]0.048328

rampCommon
Total[ms]105.385523
Average[ms]0.102916

indexCommon
Total[ms]29.523655
Average[ms]0.028832
*/

#include <chrono>
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
constexpr size_t nBuffer = 1024;

template<typename Sample> struct SmootherIndexCommon {
public:
  static void setSampleRate(Sample _sampleRate, Sample time = 0.04)
  {
    sampleRate = _sampleRate;
    setTime(time);
  }

  static void setTime(Sample seconds) { timeInSamples = seconds * sampleRate; }
  static void setBufferSize(Sample _bufferSize) { bufferSize = _bufferSize; }

  static Sample sampleRate;
  static Sample timeInSamples;
  static Sample bufferSize;
};

template<typename Sample> Sample SmootherIndexCommon<Sample>::sampleRate = 44100.0;
template<typename Sample> Sample SmootherIndexCommon<Sample>::timeInSamples = 0.0;
template<typename Sample> Sample SmootherIndexCommon<Sample>::bufferSize = 44100.0;

template<typename Sample> class SmootherIndex {
public:
  using Common = SmootherIndexCommon<Sample>;

  virtual inline Sample getValue() { return value; }
  virtual void reset(Sample value) { this->value = value; }
  virtual void refresh() { push(v0); }

  virtual void push(Sample newTarget)
  {
    v1 = v0;
    v0 = (Common::timeInSamples >= Common::bufferSize)
      ? (newTarget - v0) * Common::bufferSize / Common::timeInSamples + v0
      : newTarget;
  }

  Sample process(float index)
  {
    return value = v1 + index / Common::bufferSize * (v0 - v1);
  }

protected:
  Sample value = 1.0;
  Sample v0 = 1;
  Sample v1 = 1;
};

template<typename Sample> class SmootherIndexLocal {
public:
  void setSampleRate(Sample sampleRate, Sample time = 0.04)
  {
    this->sampleRate = sampleRate;
    setTime(time);
  }

  void setTime(Sample seconds) { timeInSamples = seconds * sampleRate; }
  void setBufferSize(Sample bufferSize) { this->bufferSize = bufferSize; }
  void reset(Sample value) { this->value = value; }
  void refresh() { push(v0); }
  inline Sample getValue() { return value; }

  void push(Sample newTarget)
  {
    v1 = v0;
    v0 = (timeInSamples >= bufferSize)
      ? (newTarget - v0) * bufferSize / timeInSamples + v0
      : newTarget;
  }

  Sample process(float index) { return value = v1 + index / bufferSize * (v0 - v1); }

protected:
  Sample sampleRate = 44100;
  Sample timeInSamples = -1;
  Sample bufferSize = 0;
  Sample v0 = 1;
  Sample v1 = 1;
  Sample value = 0;
};

template<typename Sample> class SmootherRampCommon {
public:
  static void setSampleRate(Sample _sampleRate, Sample time = 0.04)
  {
    sampleRate = _sampleRate;
    setTime(time);
  }

  static void setTime(Sample seconds) { timeInSamples = seconds * sampleRate; }
  static void setBufferSize(Sample _bufferSize) { bufferSize = _bufferSize; }

  static Sample sampleRate;
  static Sample timeInSamples;
  static Sample bufferSize;
};

template<typename Sample> Sample SmootherRampCommon<Sample>::sampleRate = 44100.0;
template<typename Sample> Sample SmootherRampCommon<Sample>::timeInSamples = 0.0;
template<typename Sample> Sample SmootherRampCommon<Sample>::bufferSize = 44100.0;

template<typename Sample> class SmootherRamp {
public:
  using Common = SmootherRampCommon<Sample>;

  virtual inline Sample getValue() { return value; }
  virtual void reset(Sample value) { this->value = value; }
  virtual void refresh() { push(target); }

  virtual void push(Sample newTarget)
  {
    target = newTarget;
    if (Common::timeInSamples < Common::bufferSize)
      value = target;
    else
      ramp = (target - value) / Common::timeInSamples;
  }

  virtual Sample process()
  {
    if (value == target) return value;
    value += ramp;

    if (fabsf(value - target) < Sample(1e-5)) value = target;
    return value;
  }

protected:
  Sample value = 1.0;
  Sample target = 1.0;
  Sample ramp = 0.0;
};

template<typename Sample> class SmootherRampLocal {
public:
  void setSampleRate(Sample sampleRate, Sample time = 0.04)
  {
    this->sampleRate = sampleRate;
    setTime(time);
  }

  void setTime(Sample seconds) { timeInSamples = seconds * sampleRate; }
  void setBufferSize(Sample bufferSize) { this->bufferSize = bufferSize; }

  void reset(Sample value)
  {
    this->value = target = value;
    ramp = 0;
  }

  void refresh() { push(target); }

  void push(Sample newTarget)
  {
    target = newTarget;
    if (timeInSamples < bufferSize)
      value = target;
    else
      ramp = (target - value) / timeInSamples;
  }

  inline Sample getValue() { return value; }

  Sample process()
  {
    if (value == target) return value;
    value += ramp;

    auto diff = value - target;
    if (diff < 0) diff = -diff;
    if (diff < 1e-5) value = target;
    return value;
  }

protected:
  Sample sampleRate = 44100;
  Sample timeInSamples = -1;
  Sample bufferSize = 0;
  Sample value = 1.0;
  Sample target = 1.0;
  Sample ramp = 0.0;
};

void benchRamp()
{
  std::vector<float> wav(nFrame * nBuffer);

  float out[nFrame];

  std::minstd_rand rng{0};
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  std::vector<SmootherRampLocal<float>> interp(64);
  for (auto &terp : interp) terp.setSampleRate(sampleRate, 0.04);

  double sumElapsed = 0.0;
  float newValue = 1.0f;
  for (size_t idx = 0; idx < nBuffer; ++idx) {
    if (idx % 16 == 0) newValue = dist(rng);
    for (auto &terp : interp) {
      terp.setBufferSize(nFrame);
      terp.push(newValue);
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < nFrame; ++i) {
      out[i] = 0;
      for (auto &terp : interp) out[i] += terp.process();
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = finish - start;
    sumElapsed += elapsed.count();

    std::memcpy(&wav[idx * nFrame], out, sizeof(float) * nFrame);
  }

  std::cout << "ramp\n"
            << "Total[ms]" << std::to_string(sumElapsed) << "\n"
            << "Average[ms]" << std::to_string(sumElapsed / nBuffer) << "\n\n";

  writeWave("snd/smootherRamp.wav", wav, sampleRate);
}

void benchIndex()
{
  std::vector<float> wav(nFrame * nBuffer);

  float out[nFrame];

  std::minstd_rand rng{0};
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  std::vector<SmootherIndexLocal<float>> interp(64);
  for (auto &terp : interp) terp.setSampleRate(sampleRate, 0.04);

  double sumElapsed = 0.0;
  float newValue = 1.0f;
  for (size_t idx = 0; idx < nBuffer; ++idx) {
    if (idx % 16 == 0) newValue = dist(rng);
    for (auto &terp : interp) {
      terp.setBufferSize(nFrame);
      terp.push(newValue);
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < nFrame; ++i) {
      out[i] = 0;
      for (auto &terp : interp) out[i] += terp.process(i);
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = finish - start;
    sumElapsed += elapsed.count();

    std::memcpy(&wav[idx * nFrame], out, sizeof(float) * nFrame);
  }

  std::cout << "index\n"
            << "Total[ms]" << std::to_string(sumElapsed) << "\n"
            << "Average[ms]" << std::to_string(sumElapsed / nBuffer) << "\n\n";

  writeWave("snd/smootherIndex.wav", wav, sampleRate);
}

void benchRampCommon()
{
  std::vector<float> wav(nFrame * nBuffer);

  float out[nFrame];

  std::minstd_rand rng{0};
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  std::vector<SmootherRamp<float>> interp(64);
  SmootherRampCommon<float>::setSampleRate(sampleRate, 0.04);

  double sumElapsed = 0.0;
  float newValue = 1.0f;
  for (size_t idx = 0; idx < nBuffer; ++idx) {
    if (idx % 16 == 0) newValue = dist(rng);
    SmootherRampCommon<float>::setBufferSize(nFrame);
    for (auto &terp : interp) terp.push(newValue);

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < nFrame; ++i) {
      out[i] = 0;
      for (auto &terp : interp) out[i] += terp.process();
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = finish - start;
    sumElapsed += elapsed.count();

    std::memcpy(&wav[idx * nFrame], out, sizeof(float) * nFrame);
  }

  std::cout << "rampCommon\n"
            << "Total[ms]" << std::to_string(sumElapsed) << "\n"
            << "Average[ms]" << std::to_string(sumElapsed / nBuffer) << "\n\n";

  writeWave("snd/smootherRampCommon.wav", wav, sampleRate);
}

void benchIndexCommon()
{
  std::vector<float> wav(nFrame * nBuffer);

  float out[nFrame];

  std::minstd_rand rng{0};
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  std::vector<SmootherIndex<float>> interp(64);
  SmootherIndexCommon<float>::setSampleRate(sampleRate, 0.04);

  double sumElapsed = 0.0;
  float newValue = 1.0f;
  for (size_t idx = 0; idx < nBuffer; ++idx) {
    if (idx % 16 == 0) newValue = dist(rng);
    SmootherIndexCommon<float>::setBufferSize(nFrame);
    for (auto &terp : interp) terp.push(newValue);

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < nFrame; ++i) {
      out[i] = 0;
      for (auto &terp : interp) out[i] += terp.process(i);
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = finish - start;
    sumElapsed += elapsed.count();

    std::memcpy(&wav[idx * nFrame], out, sizeof(float) * nFrame);
  }

  std::cout << "indexCommon\n"
            << "Total[ms]" << std::to_string(sumElapsed) << "\n"
            << "Average[ms]" << std::to_string(sumElapsed / nBuffer) << "\n\n";

  writeWave("snd/smootherIndexCommon.wav", wav, sampleRate);
}

int main()
{
  benchIndex(); // dry run.

  benchRamp();
  benchIndex();
  benchRampCommon();
  benchIndexCommon();
}
