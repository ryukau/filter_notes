#include <algorithm>
#include <array>
#include <chrono>
#include <deque>
#include <iostream>
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
  if (sf_write_float(file, &buffer[0], length) != length)
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

template<typename Sample> class Delay {
public:
  std::vector<Sample> buf;
  size_t wptr = 0;
  size_t rptr = 0;

  Delay(size_t size = 65536) : buf(size) {}

  void resize(size_t size)
  {
    buf.resize(size);
    wptr = 0;
    rptr = 0;
  }

  void reset() { std::fill(buf.begin(), buf.end(), Sample(0)); }

  void setFrames(size_t delayFrames)
  {
    if (delayFrames >= buf.size()) delayFrames = buf.size();
    rptr = wptr - delayFrames;
    if (rptr >= buf.size()) rptr += buf.size(); // Unsigned overflow case.
  }

  Sample process(Sample input)
  {
    if (++wptr >= buf.size()) wptr -= buf.size();
    buf[wptr] = input;

    if (++rptr >= buf.size()) rptr -= buf.size();
    return buf[rptr];
  }
};

template<typename T> struct RingQueue {
  std::vector<T> buf;
  size_t wptr = 0;
  size_t rptr = 0;

  void resize(size_t size) { buf.resize(size); }

  void reset(T value = 0) { std::fill(buf.begin(), buf.end(), value); }

  inline size_t size()
  {
    auto sz = wptr - rptr;
    if (sz >= buf.size()) sz += buf.size(); // Unsigned overflow case.
    return sz;
  }

  inline bool empty() { return wptr == rptr; }

  T &front() { return buf[increment(rptr)]; }
  T &back() { return buf[wptr]; }

  inline size_t increment(size_t idx)
  {
    if (++idx >= buf.size()) idx -= buf.size();
    return idx;
  }

  inline size_t decrement(size_t idx)
  {
    if (--idx >= buf.size()) idx += buf.size(); // Unsigned overflow case.
    return idx;
  }

  void push_back(T value)
  {
    wptr = increment(wptr);
    buf[wptr] = value;
  }

  T pop_front()
  {
    rptr = increment(rptr);
    return buf[rptr];
  }

  T pop_back()
  {
    wptr = decrement(wptr);
    return buf[wptr];
  }
};

/*
Ideal peak hold.

Example:
```
PeakHold<float> hold;
hold.setFrames(64);

// In processing loop.
output = hold.process(input);
```
*/
template<typename Sample> struct PeakHold {
  Sample neutral = 0;
  Delay<Sample> delay;
  RingQueue<Sample> hold;

  PeakHold(size_t size = 65536)
  {
    resize(size);
    setFrames(1); // Same as bypass.
  }

  void resize(size_t size)
  {
    delay.resize(size);
    hold.resize(size);
  }

  void reset()
  {
    delay.reset();
    hold.reset(neutral);
  }

  void setFrames(size_t frames) { delay.setFrames(frames); }

  Sample process(Sample x0)
  {
    if (!hold.empty()) {
      for (size_t idx = hold.size(); idx > 0; --idx) {
        if (hold.back() < x0)
          hold.pop_back();
        else
          break;
      }
    }

    hold.push_back(x0);

    auto delayOut = delay.process(x0);
    if (!hold.empty() && delayOut == hold.front()) hold.pop_front();

    return !hold.empty() ? hold.front() : neutral;
  }
};

template<typename Sample> struct PeakHoldDeque {
  Sample neutral = 0;
  Delay<Sample> delay;
  std::deque<Sample> hold;

  PeakHoldDeque(size_t size = 65536) { resize(size); }

  void resize(size_t size) { delay.resize(size); }

  void reset()
  {
    delay.reset();
    hold.resize(0);
  }

  void setFrames(uint32_t frames) { delay.setFrames(frames); }

  Sample process(Sample x0)
  {
    if (!hold.empty()) {
      for (size_t idx = hold.size(); idx > 0; --idx) {
        if (hold.back() < x0)
          hold.pop_back();
        else
          break;
      }
    }

    hold.push_back(x0);

    auto delayOut = delay.process(x0);
    if (!hold.empty() && delayOut == hold.front()) hold.pop_front();

    return !hold.empty() ? hold.front() : neutral;
  }
};

template<typename Sample> struct PeakHoldNaive {
  Delay<Sample> delay;

  PeakHoldNaive(size_t size = 65536)
  {
    delay.buf.reserve(size);
    setFrames(1);
  }

  void reset() { delay.reset(); }
  void setFrames(uint32_t frames) { delay.resize(frames); }

  Sample process(Sample x0)
  {
    delay.process(x0);

    // std::max_element is slow.
    Sample max = 0;
    for (const auto &value : delay.buf) {
      if (max < value) max = value;
    }
    return max;
  }
};

template<typename Sample> struct ForwardHold {
  Sample hold = 0;
  Sample neutral = 0;
  size_t holdFrames = 1;
  size_t counter = 0;

  void setFrames(size_t frames) { holdFrames = frames; }

  void reset()
  {
    counter = 0;
    hold = neutral;
  }

  Sample process(Sample input)
  {
    if (counter > 0)
      --counter;
    else
      hold = neutral;

    if (hold <= input) {
      hold = input;
      counter = holdFrames;
    }

    return hold;
  }
};

constexpr size_t nLoop = 4; // Set 2 or greater to test `reset` method.

template<typename Hold> void bench(std::string name, SoundFile &snd)
{
  const auto &signal = snd.data;
  std::vector<float> wav(signal.size());

  Hold hold;
  hold.setFrames(64);

  double sumElapsed = 0.0;
  for (size_t loop = 0; loop < nLoop; ++loop) {
    hold.reset();

    for (size_t idx = 0; idx < wav.size(); ++idx) {
      auto start = std::chrono::steady_clock::now();

      wav[idx] = hold.process(signal[idx]);

      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end - start;
      sumElapsed += elapsed.count();
    }
  }
  std::cout << name << ": " << signal.size() << "[sample], " << sumElapsed << "[ms], "
            << "\n";

  name = "snd/" + name + ".wav";
  writeWave(name, wav, snd.samplerate);
}

int main()
{
  SoundFile snd("snd/input.wav");

  std::cout << "--- Warm up\n";
  bench<ForwardHold<float>>("Forward", snd);

  std::cout << "\n--- Benchmark\n";
  bench<PeakHold<float>>("Ideal", snd);
  bench<PeakHoldDeque<float>>("IdealDeque", snd);
  // bench<PeakHoldNaive<float>>("Naive", snd);

  return 0;
}
