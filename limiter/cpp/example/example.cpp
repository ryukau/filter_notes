#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
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

template<typename Sample> struct PeakHold {
  Sample neutral = 0;
  Delay<Sample> delay;
  RingQueue<Sample> hold;

  PeakHold(size_t size = 65536)
  {
    resize(size);
    setFrames(1);
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

template<typename Sample> struct DoubleAverageFilter {
  Sample sum1 = 0;
  Sample sum2 = 0;
  Sample buf = 0; // 出力が 1 サンプル前にずれるのを補正するディレイのバッファ。
  size_t halfDelayFrames = 0;
  Delay<Sample> delay1;
  Delay<Sample> delay2;

  void resize(size_t size)
  {
    delay1.resize(size / 2);
    delay2.resize(size / 2);
  }

  void reset()
  {
    sum1 = 0;
    sum2 = 0;
    buf = 0;
    delay1.reset();
    delay2.reset();
  }

  void setFrames(size_t frames)
  {
    halfDelayFrames = frames / 2;
    delay1.setFrames(halfDelayFrames);
    delay2.setFrames(halfDelayFrames);
  }

  Sample process(const Sample input)
  {
    sum1 += buf - delay1.process(buf);
    auto out1 = sum1 / halfDelayFrames;

    sum2 += out1 - delay2.process(out1);
    auto out2 = sum2 / halfDelayFrames;

    buf = input;
    return out2;
  }
};

template<typename Sample> struct Limiter {
  static constexpr Sample fixedGain = Sample(0.9965520801347684); // -0.03dB.
  static constexpr Sample releaseConstant = Sample(1e-5); // 適当な小さい値。

  Sample threshold = Sample(0.1);
  Sample gain = Sample(1);
  Sample release = 0; // 指数関数的増加のためのリリース係数。
  size_t attackFrames = 0;

  PeakHold<Sample> hold;
  DoubleAverageFilter<Sample> smoother;
  Delay<Sample> lookaheadDelay;

  size_t latency() { return attackFrames; }

  void resize(size_t size)
  {
    size += size % 2;
    hold.resize(size);
    smoother.resize(size);
    lookaheadDelay.resize(size);
  }

  void reset()
  {
    gain = Sample(1);
    hold.reset();
    smoother.reset();
    lookaheadDelay.reset();
  }

  void prepare(
    Sample sampleRate,
    Sample attackSeconds,
    Sample sustainSeconds,
    Sample releaseSeconds,
    Sample threshold)
  {
    auto prevAttack = attackFrames;
    attackFrames = size_t(sampleRate * attackSeconds);
    attackFrames += attackFrames % 2; // DoubleAverageFilter のために 2 の倍数にする。
    if (prevAttack != attackFrames) reset();

    release
      = std::pow(Sample(1 / releaseConstant), Sample(1 / (releaseSeconds * sampleRate)));

    this->threshold = threshold;

    hold.setFrames(attackFrames + size_t(sampleRate * sustainSeconds));
    smoother.setFrames(attackFrames);
    lookaheadDelay.setFrames(attackFrames);
  }

  inline Sample applyCharacteristicCurve(Sample x0)
  {
    return x0 > threshold ? threshold / x0 : Sample(1);
  }

  inline Sample softClip(Sample x0, Sample ratio)
  {
    const auto absed = std::fabs(x0);

    const auto a1 = threshold * ratio;
    if (absed <= a1) return x0;

    const auto a2 = 2 * threshold - a1;
    if (absed >= a2) return threshold;

    return std::copysign(
      threshold + (a2 - absed) * (a2 - absed) * Sample(0.25) / (a1 - threshold), x0);
  }

  Sample process(const Sample input)
  {
    auto holdGain = hold.process(std::fabs(input));
    auto candidate = applyCharacteristicCurve(holdGain);
    gain = std::min(gain * release, candidate);

    auto smoothed = smoother.process(gain);
    auto delayed = lookaheadDelay.process(input);
    return softClip(smoothed * delayed, fixedGain);
  }
};

int main()
{
  constexpr size_t sampleRate = 48000;
  constexpr float maxPeak = 10.0f;

  // テスト信号の生成。
  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_real_distribution<float> dist(-maxPeak, maxPeak);
  std::vector<float> input(sampleRate);
  for (size_t i = 0; i < input.size(); ++i) input[i] = dist(rng);

  // Limiter の使用例。
  Limiter<float> limiter;
  limiter.resize(65536);
  limiter.prepare(sampleRate, 0.002f, 0.002f, 0.1f, 0.5f);

  std::vector<float> output(input.size());
  for (size_t i = 0; i < input.size(); ++i) output[i] = limiter.process(input[i]);

  // 出力がしきい値以下に制限されているか確認。
  float max = 0;
  for (const auto &value : output) {
    auto absed = std::fabs(value);
    if (absed > limiter.threshold && absed > max) max = absed;
  }

  if (max == 0) {
    std::cout << "Limiting succeeded.\n";
  } else {
    std::cout << "Limiting failed.\n"
              << std::setprecision(std::numeric_limits<float>::digits10 + 1)
              << "threshold: " << limiter.threshold << "\n"
              << "max      : " << max << "\n";
  }
}
