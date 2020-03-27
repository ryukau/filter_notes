#include "exponential_envelope.hpp"

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
  ExpADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.reset(attackTime, decayTime, sustainLevel, releaseTime, noteFreq, 0);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(releaseAt * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("snd/ADSR.wav", buffer, sampleRate);
}

void testReleaseWhileAttack()
{
  ExpADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.reset(4.0f, 1.0f, 0.5f, 2.0f, 100.0f, 0.5f);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(3.5f * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("snd/ReleaseWhileAttack.wav", buffer, sampleRate);
}

void testReleaseWhileDecay()
{
  ExpADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.reset(1.0f, 1.0f, 0.5f, 0.5f, 100.0f, 0);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(1.5f * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("snd/ReleaseWhileDecay.wav", buffer, sampleRate);
}

void testTriggerWhileRelease()
{
  ExpADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.reset(1.0f, 1.0f, 0.5f, 8.0f, 100.0f, -1.0f);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(3.0f * sampleRate)) envelope.release();
    if (idx == size_t(4.0f * sampleRate))
      envelope.reset(0.5f, 0.5f, 0.3f, 2.0f, 100.0f, 2.0f);
    if (idx == size_t(6.0f * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("snd/TriggerWhileRelease.wav", buffer, sampleRate);
}

void testChangeSustain()
{
  const float atk = 0.1f;
  const float dec = 0.1f;
  const float rel = 0.01f;
  const float noteFreq = 100.0f;

  ExpADSREnvelope<float> envelope;
  envelope.setup(sampleRate);
  envelope.reset(atk, dec, 0.5f, rel, noteFreq, 0);

  std::vector<float> buffer(size_t(8 * sampleRate));
  for (size_t idx = 0; idx < buffer.size(); ++idx) {
    if (idx == size_t(1.0f * sampleRate)) envelope.set(atk, dec, 0.3f, rel, noteFreq);
    if (idx == size_t(2.0f * sampleRate)) envelope.set(atk, dec, 1.3f, rel, noteFreq);
    if (idx == size_t(3.0f * sampleRate)) envelope.set(atk, dec, -0.3f, rel, noteFreq);
    if (idx == size_t(4.0f * sampleRate)) envelope.set(atk, dec, 0.5f, rel, noteFreq);
    if (idx == size_t(5.0f * sampleRate)) envelope.set(atk, dec, 0.1f, rel, noteFreq);
    if (idx == size_t(6.0f * sampleRate)) envelope.release();
    buffer[idx] = envelope.process();
  }
  writeWave("snd/ChangeSustain.wav", buffer, sampleRate);
}

int main()
{
  testADSR();
  testReleaseWhileAttack();
  testReleaseWhileDecay();
  testTriggerWhileRelease();
  testChangeSustain();
  return 0;
}
