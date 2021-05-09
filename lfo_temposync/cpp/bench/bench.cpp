#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <deque>
#include <iostream>
#include <limits>
#include <sndfile.h>
#include <string>
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
  if (sf_write_float(file, &buffer[0], length) != (int)length)
    std::cout << sf_strerror(file) << std::endl;

  if (sf_close(file) != 0) {
    std::cout << "Error: sf_close failed." << std::endl;
    exit(EXIT_FAILURE);
  }
}

struct SoundFile {
  int samplerate = 0;
  int channels = 0;
  size_t frames = 0;
  std::vector<float> data;

  SoundFile(std::string path)
  {
    SF_INFO sfinfo;
    memset(&sfinfo, 0, sizeof(sfinfo));

    SNDFILE *file = sf_open(path.c_str(), SFM_READ, &sfinfo);
    if (!file) {
      std::cout << "Error: sf_open failed." << std::endl;
      // exit(EXIT_FAILURE);
    }

    samplerate = sfinfo.samplerate;
    channels = sfinfo.channels;
    frames = sfinfo.frames;

    size_t items = channels * frames;
    std::vector<float> raw(items);
    sf_read_float(file, &raw[0], items);

    data.resize(frames);
    for (size_t i = 0; i < data.size(); ++i) data[i] = raw[channels * i];

    if (sf_close(file) != 0) {
      std::cout << "Error: sf_close failed." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
};

/**
v* for velocity, p* for position.
- p1 is true phase position.
- v1 is true phase increment in a frame.
- p2 is a phase position which chases p1 while transitioning.
- v2 is a temporary phase increment while braking.
*/
template<typename Sample> class TempoSynchronizer {
private:
  enum class State { steady, decelerating, accelerating };

  State state = State::steady;
  Sample v1 = 0;
  Sample p1 = 0;
  Sample v2 = 0;
  Sample p2 = 0;
  Sample lastTempo = 0;
  Sample lastSync = 0;

  Sample midTime = 0; // In samples.
  Sample midVelocity = 0;
  Sample counter = 0;

public:
  void reset(Sample sampleRate, Sample tempo, Sample sync)
  {
    v1 = tempo / (Sample(60) * sampleRate * sync);
    p1 = 0;
    v2 = v1;
    p2 = 0;
    lastTempo = tempo;
    lastSync = sync;
  }

  // Must call this method at the start of each DSP processing cycle.
  void prepare(Sample sampleRate, Sample tempo, Sample sync, Sample beatsElapsed)
  {
    v1 = tempo / (Sample(60) * sampleRate * sync);
    p1 = beatsElapsed / sync;
    p1 -= std::floor(p1);

    if (lastTempo != tempo || lastSync != sync) {
      v2 = lastTempo / (Sample(60) * sampleRate * lastSync);

      if (state == State::steady) {
        p2 = beatsElapsed / lastSync;
        p2 -= std::floor(p2);
      }

      state = State::decelerating;

      midTime = Sample(0.5) / v1;
      auto distance = p1 - p2;
      auto k = std::ceil((v1 + v2) * Sample(0.5) * midTime - distance);
      midVelocity = (distance + k) / midTime - (v1 + v2) * Sample(0.5);

      counter = 0;
    }
    lastTempo = tempo;
    lastSync = sync;
  }

  Sample process(size_t &st)
  {
    st = static_cast<size_t>(state); // debug

    Sample outPhase;
    switch (state) {
      default:
      case State::steady: {
        outPhase = p1;

        p1 += v1;
        p1 -= std::floor(p1);
      } break;

      case State::decelerating: {
        outPhase = p2;

        p2 += v2 + (midVelocity - v2) * counter / midTime;
        p2 -= std::floor(p2);

        if (++counter >= midTime) {
          state = State::accelerating;
          counter = 0;
        }
      } break;

      case State::accelerating: {
        outPhase = p2;

        p2 += midVelocity + (v1 - midVelocity) * counter / midTime;
        p2 -= std::floor(p2);

        if (++counter >= midTime) {
          state = State::steady;
          p1 = p2;
        }
      } break;
    }
    return outPhase;
  }
};

template<typename Sample> class LineGenerator {
private:
  struct Point {
    size_t frame = 0;
    Sample value = 0;
  };

  inline void addPoint(size_t frame, Sample value)
  {
    Point pt;
    pt.frame = frame;
    pt.value = value;
    points.push_back(pt);
  }

public:
  enum class FillType { step, slope };

  // For each point:
  // - Index 0 is time in samples.
  // - Index 1 is output value.
  std::deque<Point> points;

  LineGenerator(Sample sampleRate, Sample interval, std::vector<Sample> &values)
  {
    for (size_t idx = 0; idx < values.size(); ++idx)
      addPoint(size_t(Sample(idx) * sampleRate * interval), values[idx]);
  }

  Sample process(size_t timeElapsedInSamples)
  {
    if (points.size() <= 0) return Sample(0);
    if (points.size() == 1) return points[0].value;
    if (points[1].frame <= timeElapsedInSamples) points.pop_front();

    if (points.size() <= 0) return Sample(0);
    auto fraction = Sample(timeElapsedInSamples - points[0].frame)
      / Sample(points[1].frame - points[0].frame);
    return points[0].value + fraction * (points[1].value - points[0].value);
  }
};

template<typename Sample> class StepGenerator {
private:
  struct Point {
    size_t frame = 0;
    Sample value = 0;
  };

  inline void addPoint(size_t frame, Sample value)
  {
    Point pt;
    pt.frame = frame;
    pt.value = value;
    points.push_back(pt);
  }

public:
  // For each point:
  // - Index 0 is time in samples.
  // - Index 1 is output value.
  std::deque<Point> points;

  StepGenerator(Sample sampleRate, Sample interval, std::vector<Sample> &values)
  {
    for (size_t idx = 0; idx < values.size(); ++idx)
      addPoint(size_t(Sample(idx) * sampleRate * interval), values[idx]);
  }

  Sample process(size_t timeElapsedInSamples)
  {
    if (points.size() <= 0) return Sample(0);
    if (points.size() == 1) return points[0].value;
    if (points[1].frame <= timeElapsedInSamples) points.pop_front();
    if (points.size() <= 0) return Sample(0);
    return points[0].value;
  }
};

template<typename Synchronizer> void bench(std::string name)
{
  constexpr float sampleRate = 48000.0f;
  constexpr size_t nLoop = 1;
  constexpr size_t bufferSize = 512;

  constexpr float carrierFreq = 1000.0f;
  constexpr float carrierDelta = float(twopi) * carrierFreq / sampleRate;

  constexpr float lineInterval = 1.0f;

  std::vector<float> tempoValues
    = {120.0f, 120.0f, 120.0f, 120.0f, 120.0f, 120.0f, 120.0f, 120.0f, 120.0f};
  StepGenerator<float> tempoGen(sampleRate, lineInterval, tempoValues);

  std::vector<float> syncValues
    = {0.25f, 0.25f, 8.0f, 8.0f, 8.0f, 1.0f, 1.0f, 1.0f, 1.0f};
  StepGenerator<float> syncGen(sampleRate, lineInterval, syncValues);

  if (tempoValues.size() != syncValues.size()) {
    std::cerr << "Size mismatch of tempoValues and syncValues.\n";
    return;
  }

  Synchronizer syncer;

  double sumElapsed = 0.0;

  size_t wavSize = (tempoValues.size() - 1) * size_t(sampleRate);
  wavSize = bufferSize * (wavSize / bufferSize);

  std::vector<float> wav;
  wav.reserve(wavSize);

  std::vector<float> sinwave;
  sinwave.reserve(wavSize);

  std::vector<float> state;
  state.reserve(wavSize);

  std::vector<float> syncout;
  syncout.reserve(wavSize);

  std::vector<float> amout;
  amout.reserve(wavSize);

  float carrierPhase = -carrierDelta;
  for (size_t loop = 0; loop < nLoop; ++loop) {
    size_t frame = 0;
    double beatsElapsed = 0.0f;
    syncer.reset(sampleRate, tempoGen.process(frame), syncGen.process(frame));
    while (frame < wavSize) {
      auto sync = syncGen.process(frame);
      auto tempo = tempoGen.process(frame);
      syncer.prepare(sampleRate, tempo, sync, float(beatsElapsed));
      for (size_t i = 0; i < bufferSize; ++i) {
        carrierPhase += carrierDelta;
        if (carrierPhase >= float(twopi)) carrierPhase -= float(twopi);
        auto carrierSig = std::sin(carrierPhase);

        auto start = std::chrono::steady_clock::now();

        size_t st = 0;

        float phase = syncer.process(st);
        wav.push_back(phase);
        sinwave.push_back(std::sin(float(twopi) * phase));
        state.push_back(float(st));
        syncout.push_back(sync);
        amout.push_back(carrierSig * 0.5f * (std::sin(float(twopi) * phase) + 1.0f));

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        sumElapsed += elapsed.count();

        ++frame;
      }
      beatsElapsed += tempo * bufferSize / (60.0 * sampleRate);
    }
  }
  std::cout << name << ": " << wav.size() << "[sample], " << sumElapsed << "[ms], "
            << "\n";

  std::string filename = "snd/step_" + name + ".wav";
  writeWave(filename, wav, size_t(sampleRate));

  filename = "snd/sinwave_" + name + ".wav";
  writeWave(filename, sinwave, size_t(sampleRate));

  filename = "snd/state_" + name + ".wav";
  writeWave(filename, state, size_t(sampleRate));

  filename = "snd/syncout_" + name + ".wav";
  writeWave(filename, syncout, size_t(sampleRate));

  filename = "snd/amout_" + name + ".wav";
  writeWave(filename, amout, size_t(sampleRate));
}

int main()
{
  std::cout << "--- Warm up\n";
  bench<TempoSynchronizer<float>>("Sync");

  return 0;
}
