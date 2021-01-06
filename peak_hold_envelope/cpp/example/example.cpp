#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

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
Edge cases:
- All output become 0 when `setFrames(0)`.
- Bypass input when `setFrames(1)`.
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

int main()
{
  constexpr size_t length = 64;

  PeakHold<float> hold;
  hold.setFrames(4);

  std::minstd_rand rng(0);
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  for (size_t idx = 0; idx < length; ++idx) {
    auto input = dist(rng);
    std::cout << input << ", " << hold.process(input) << "\n";
  }

  return 0;
}
