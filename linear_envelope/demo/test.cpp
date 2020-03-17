#include "linear_envelope.hpp"

#include <cmath>
#include <iostream>
#include <sndfile.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

using Sample = float;

constexpr Sample sampleRate = 44100.0;

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

void testADSR(
  float attackTime = 1.0f,
  float decayTime = 1.0f,
  float sustainLevel = 0.5f,
  float releaseTime = 2.0f,
  float releaseAt = 3.0f,
  float noteFreq = 100.0f)
{
  LinearADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.trigger(attackTime, decayTime, sustainLevel, releaseTime, noteFreq);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(releaseAt * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("snd/ADSR.wav", buffer, sampleRate);
}

void testReleaseWhileAttack()
{
  LinearADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.trigger(4.0f, 1.0f, 0.5f, 2.0f, 100.0f);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(3.0f * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("snd/ReleaseWhileAttack.wav", buffer, sampleRate);
}

void testReleaseWhileDecay()
{
  LinearADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.trigger(1.0f, 1.0f, 0.5f, 0.5f, 100.0f);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(1.5f * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("snd/ReleaseWhileDecay.wav", buffer, sampleRate);
}

void testTriggerWhileRelease()
{
  LinearADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.trigger(1.0f, 1.0f, 0.5f, 2.0f, 100.0f);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(3.0f * sampleRate)) envelope.release();
    if (idx == size_t(4.0f * sampleRate))
      envelope.trigger(0.5f, 0.5f, 0.3f, 2.0f, 100.0f);
    if (idx == size_t(6.0f * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("snd/TriggerWhileRelease.wav", buffer, sampleRate);
}

template<typename Sample> struct PController {
  Sample kp = 1; // Range in [0, 1].
  Sample value = 0;

  void setup(Sample p) { kp = std::clamp<Sample>(p, Sample(0), Sample(1)); }
  void reset() { value = 0; }
  Sample process(Sample input) { return value += kp * (input - value); }
};

void testPController(
  float attackTime = 1.0f,
  float decayTime = 1.0f,
  float sustainLevel = 0.5f,
  float releaseTime = 2.0f,
  float releaseAt = 3.0f,
  float noteFreq = 100.0f)
{
  LinearADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.trigger(attackTime, decayTime, sustainLevel, releaseTime, noteFreq);

  PController<float> pController;
  pController.setup(16 / sampleRate);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(releaseAt * sampleRate)) envelope.release();
    buffer[idx] = pController.process(envelope.process());
  }
  writeWave("snd/PController.wav", buffer, sampleRate);
}

int main()
{
  testADSR();
  testReleaseWhileAttack();
  testReleaseWhileDecay();
  testTriggerWhileRelease();
  testPController();
  return 0;
}
