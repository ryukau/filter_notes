/*
Reference: https://signalsmith-audio.co.uk/writing/2022/constant-time-peak-hold/

Below is the license of `SignalsmithAudio::PeakHold`.

```
MIT License

Copyright (c) 2021 Geraint Luff / Signalsmith Audio Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
*/

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <format>
#include <iostream>
#include <limits>
#include <random>
#include <sndfile.h>
#include <string>
#include <vector>

void writeWave(
  std::string filename, const std::vector<float> &buffer, const int samplerate)
{
  SF_INFO sfinfo;
  memset(&sfinfo, 0, sizeof(sfinfo));
  sfinfo.samplerate = samplerate;
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

namespace Uhhyou {

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
Peak hold for limiter.
- Only outputs 0 when `setFrames(0)`.
- Bypass input when `setFrames(1)`.
*/
template<typename Sample> struct PeakHold {
  std::string name = "Run-length encoding";

  Sample neutral = 0;
  Delay<Sample> delay;
  RingQueue<Sample> queue;

  PeakHold(size_t size = 65536)
  {
    resize(size);
    set(1);
  }

  void resize(size_t size)
  {
    delay.resize(size);
    queue.resize(size);
  }

  void reset()
  {
    delay.reset();
    queue.reset(neutral);
  }

  void set(size_t frames) { delay.setFrames(std::min(frames, delay.buf.size() - 1)); }

  Sample process(Sample x0)
  {
    while (!queue.empty()) {
      if (queue.back() >= x0) break;
      queue.pop_back();
    }
    queue.push_back(x0);
    if (delay.process(x0) == queue.front()) queue.pop_front();
    return queue.front();
  }

  Sample operator()(Sample v) { return process(v); }
};

} // namespace Uhhyou

namespace SignalsmithAudio {

/** Peak-hold filter.
        \diagram{peak-hold.svg}

        The size is variable, and can be changed instantly with `.set()`, or by using
   `.push()`/`.pop()` in an unbalanced way.

        This has complexity O(1) every sample when the length remains constant (balanced
   `.push()`/`.pop()`, or using `filter(v)`), and amortised O(1) complexity otherwise.  To
   avoid allocations while running, it pre-allocates a vector (not a `std::deque`) which
   determines the maximum length.
*/
template<typename Sample> class PeakHold {
  static constexpr Sample lowest = std::numeric_limits<Sample>::lowest();
  int bufferMask;
  std::vector<Sample> buffer;
  int backIndex = 0, middleStart = 0, workingIndex = 0, middleEnd = 0, frontIndex = 0;
  Sample frontMax = lowest, workingMax = lowest, middleMax = lowest;

public:
  std::string name = "Constant-time";

  PeakHold(int maxLength) { resize(maxLength); }

  int size() { return frontIndex - backIndex; }

  void resize(int maxLength)
  {
    int bufferLength = 1;
    while (bufferLength <= maxLength) bufferLength *= 2;
    buffer.resize(bufferLength);
    bufferMask = bufferLength - 1;

    frontIndex = backIndex + maxLength;
    reset();
  }

  void reset(Sample fill = lowest)
  {
    int prevSize = size();
    buffer.assign(buffer.size(), fill);
    frontMax = workingMax = middleMax = lowest;
    middleEnd = workingIndex = frontIndex = 0;
    middleStart = middleEnd - (prevSize / 2);
    backIndex = frontIndex - prevSize;
  }

  /** Sets the size immediately.
  Must be `0 <= newSize <= maxLength` (see constructor and `.resize()`).

  Shrinking doesn't destroy information, and if you expand again (with
  `preserveCurrentPeak=false`), you will get the same output as before shrinking.
  Expanding when `preserveCurrentPeak` is enabled is destructive, re-writing its history
  such that the current output value is unchanged.*/
  void set(int newSize, bool preserveCurrentPeak = false)
  {
    while (size() < newSize) {
      Sample &backPrev = buffer[backIndex & bufferMask];
      --backIndex;
      Sample &back = buffer[backIndex & bufferMask];
      back = preserveCurrentPeak ? backPrev : std::max(back, backPrev);
    }
    while (size() > newSize) {
      pop();
    }
  }

  void push(Sample v)
  {
    buffer[frontIndex & bufferMask] = v;
    ++frontIndex;
    frontMax = std::max(frontMax, v);
  }

  void pop()
  {
    if (backIndex == middleStart) {
      // Move along the maximums
      workingMax = lowest;
      middleMax = frontMax;
      frontMax = lowest;

      int prevFrontLength = frontIndex - middleEnd;
      int prevMiddleLength = middleEnd - middleStart;
      if (prevFrontLength <= prevMiddleLength + 1) {
        // Swap over simply
        middleStart = middleEnd;
        middleEnd = frontIndex;
        workingIndex = middleEnd;
      } else {
        // The front is longer than the middle - only happens if unbalanced
        // We don't move *all* of the front over, keeping half the surplus in the front
        int middleLength = (frontIndex - middleStart) / 2;
        middleStart = middleEnd;
        middleEnd += middleLength;

        // Working index is close enough that it will be finished by the time the back is
        // empty
        int backLength = middleStart - backIndex;
        int workingLength = std::min(backLength, middleEnd - middleStart);
        workingIndex = middleStart + workingLength;

        // Since the front was not completely consumed, we re-calculate the front's
        // maximum
        for (int i = middleEnd; i != frontIndex; ++i) {
          frontMax = std::max(frontMax, buffer[i & bufferMask]);
        }
        // The index might not start at the end of the working block - compute the last
        // bit immediately
        for (int i = middleEnd - 1; i != workingIndex - 1; --i) {
          buffer[i & bufferMask] = workingMax
            = std::max(workingMax, buffer[i & bufferMask]);
        }
      }

      // Is the new back (previous middle) empty? Only happens if unbalanced
      if (backIndex == middleStart) {
        // swap over again (front's empty, no change)
        workingMax = lowest;
        middleMax = frontMax;
        frontMax = lowest;
        middleStart = workingIndex = middleEnd;

        if (backIndex == middleStart) {
          --backIndex; // Only happens if you pop from an empty list - fail nicely
        }
      }

      buffer[frontIndex & bufferMask]
        = lowest; // In case of length 0, when everything points at this value
    }

    ++backIndex;
    if (workingIndex != middleStart) {
      --workingIndex;
      buffer[workingIndex & bufferMask] = workingMax
        = std::max(workingMax, buffer[workingIndex & bufferMask]);
    }
  }
  Sample read()
  {
    Sample backMax = buffer[backIndex & bufferMask];
    return std::max(backMax, std::max(middleMax, frontMax));
  }

  // For simple use as a constant-length filter
  Sample operator()(Sample v)
  {
    push(v);
    pop();
    return read();
  }
};

} // namespace SignalsmithAudio

void render()
{
  SoundFile sound("snd/input.wav");
  SignalsmithAudio::PeakHold<double> hold{31};
  std::vector<float> peak(sound.frames);
  for (size_t i = 0; i < peak.size(); ++i) peak[i] = (float)hold(sound.data[i]);
  writeWave("snd/output_peak.wav", peak, 48000);
}

#define ROW_FORMAT_STR "|{:20}|{:20}|{:20}|{:20}|\n"

template<typename PeakHold, typename Sample> void benchSerial(int seed)
{
  constexpr Sample sampleRate = Sample(48000);
  constexpr int holdSize = 11;
  constexpr size_t nSample = 200 * size_t(sampleRate);

  PeakHold hold{65536};
  hold.set(holdSize);
  auto rng = std::minstd_rand(seed);
  std::uniform_real_distribution<Sample> dist(Sample(-1), Sample(1));
  double sumElapsed = 0.0;
  for (size_t i = 0; i < nSample; ++i) {
    auto start = std::chrono::steady_clock::now();

    hold(dist(rng));

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }

  std::cout << std::format(ROW_FORMAT_STR, hold.name, sumElapsed, nSample, hold(0));
}

void benchmark()
{
  std::cout << std::format(
    ROW_FORMAT_STR, "Filter Type", "Time Elapsed [ms]", "Sample Count", "Last Output");

  constexpr int seed = 4489;
  benchSerial<Uhhyou::PeakHold<double>, double>(seed);
  benchSerial<SignalsmithAudio::PeakHold<double>, double>(seed);
}

int main()
{
  // render();
  benchmark();
  return 0;
}
