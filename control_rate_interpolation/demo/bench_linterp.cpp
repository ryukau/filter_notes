/*
$ g++ -lsndfile -O3 bench.cpp; ./a.out
ramp
Total[ms]40.644869
Average[ms]0.039692

index
Total[ms]20.039410
Average[ms]0.019570

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

struct LinearInterpRamp {
  float target = 0;
  float ramp = 0;
  float value = 0;
  size_t frame = 0;

  LinearInterpRamp(float value = 0) : target(value), value(value) {}

  void push(size_t frame, float newTarget)
  {
    target = newTarget;
    ramp = (target - value) / frame;
  }

  float process()
  {
    if (value == target) return value;
    value += ramp;
    if (fabsf(value - target) < 1e-5) value = target;
    return value;
  }
};

struct LinearInterpIndex {
  float v0 = 0;
  float v1 = 0;

  LinearInterpIndex(float value = 0) : v0(value), v1(value) {}

  void push(float newTarget)
  {
    v1 = v0;
    v0 = newTarget;
  }

  float process(float frame, float index) { return v1 + index / frame * (v0 - v1); }
};

void benchRamp()
{
  std::vector<float> wav(nFrame * nBuffer);

  float out[nFrame];

  std::minstd_rand rng{0};
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  std::vector<LinearInterpRamp> interp(64);

  double sumElapsed = 0.0;
  float newValue = 1.0f;
  for (size_t idx = 0; idx < nBuffer; ++idx) {
    if (idx % 16 == 0) newValue = dist(rng);
    for (auto &terp : interp) terp.push(nFrame, newValue);

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

  writeWave("snd/linterpRamp.wav", wav, sampleRate);
}

void benchIndex()
{
  std::vector<float> wav(nFrame * nBuffer);

  float out[nFrame];

  std::minstd_rand rng{0};
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  std::vector<LinearInterpIndex> interp(64);

  double sumElapsed = 0.0;
  float newValue = 1.0f;
  for (size_t idx = 0; idx < nBuffer; ++idx) {
    if (idx % 16 == 0) newValue = dist(rng);
    for (auto &terp : interp) terp.push(newValue);

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < nFrame; ++i) {
      out[i] = 0;
      for (auto &terp : interp) out[i] += terp.process(nFrame, i);
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = finish - start;
    sumElapsed += elapsed.count();

    std::memcpy(&wav[idx * nFrame], out, sizeof(float) * nFrame);
  }

  std::cout << "index\n"
            << "Total[ms]" << std::to_string(sumElapsed) << "\n"
            << "Average[ms]" << std::to_string(sumElapsed / nBuffer) << "\n\n";

  writeWave("snd/linterpIndex.wav", wav, sampleRate);
}

int main()
{
  benchRamp();
  benchIndex();
}
