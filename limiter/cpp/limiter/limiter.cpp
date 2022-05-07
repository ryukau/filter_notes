#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <random>
#include <sndfile.h>
#include <string>
#include <vector>

void writeWave(
  std::string &&filename, std::vector<float> &buffer, const size_t &samplerate)
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
Peak hold for limiter.
- Only outputs 0 when `setFrames(0)`.
- Bypass input when `setFrames(1)`.
*/
template<typename Sample> struct PeakHold {
  Sample neutral = 0;
  Delay<Sample> delay;
  RingQueue<Sample> queue;

  PeakHold(size_t size = 65536)
  {
    resize(size);
    setFrames(1);
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

  void setFrames(size_t frames)
  {
    delay.setFrames(std::min(frames, delay.buf.size() - 1));
  }

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
};

template<typename Sample> struct DoubleAverageFilter {
  Sample denom = Sample(1);
  Sample sum1 = 0;
  Sample sum2 = 0;
  Sample buf = 0;
  Delay<Sample> delay1;
  Delay<Sample> delay2;

  void resize(size_t size)
  {
    delay1.resize(size / 2 + 1);
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
    auto &&half = frames / 2;
    denom = 1 / Sample((half + 1) * half);
    delay1.setFrames(half + 1);
    delay2.setFrames(half);
  }

  inline Sample add(Sample lhs, Sample rhs)
  {
    bool &&swapped = lhs < rhs;
    if (swapped) std::swap(lhs, rhs);
    int expL;
    std::frexp(lhs, &expL);
    auto &&cut = std::ldexp(float(1), expL - std::numeric_limits<Sample>::digits);
    auto &&rounded = rhs - std::fmod(rhs, cut);
    return lhs + rounded;
  }

  inline Sample sub(Sample lhs, Sample rhs)
  {
    bool &&swapped = lhs < rhs;
    if (swapped) std::swap(lhs, rhs);
    int expL;
    std::frexp(lhs, &expL);
    auto &&cut = std::ldexp(float(1), expL - std::numeric_limits<Sample>::digits + 1);
    auto &&rounded = rhs - std::fmod(rhs, cut);
    return swapped ? rounded - lhs : lhs - rounded;
  }

  // Only tested where `std::numeric_limits<float>::round_style == std::round_to_nearest`.
  // add() limits output. Useful.
  // sub() breaks limiting. Useless otherwise adding some strange FX.
  Sample process(Sample input)
  {
    input *= denom;

    sum1 = add(sum1, input);
    // sum1 += input;
    Sample d1 = delay1.process(input);
    // sum1 = sub(sum1, d1);
    sum1 = std::fmax(Sample(0), sum1 - d1);

    sum2 = add(sum2, sum1);
    // sum2 += sum1;
    Sample d2 = delay2.process(sum1);
    // sum2 = sub(sum2, d2);
    sum2 = std::fmax(Sample(0), sum2 - d2);

    auto output = buf;
    buf = sum2;
    return output;
  }

  Sample processNaive(Sample input)
  {
    input *= denom;
    sum1 += input - delay1.process(input);
    sum2 += sum1 - delay2.process(sum1);
    auto output = buf;
    buf = sum2;
    return output;
  }
};

template<typename Sample, typename Int> struct DoubleAverageFilterInt {
  static constexpr Sample scale = Sample(1 << std::numeric_limits<Sample>::digits);
  // static constexpr Sample scale = Sample(1 << 8); // This somehow works.

  Sample denom = 1;
  Int sum1 = 0;
  Int sum2 = 0;
  Int buf = 0;
  Delay<Int> delay1;
  Delay<Int> delay2;

  void resize(size_t size)
  {
    delay1.resize(size / 2 + 1);
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
    auto &&half = frames / 2;
    denom = 1 / Sample((half + 1) * half * scale);
    delay1.setFrames(half + 1);
    delay2.setFrames(half);
  }

  inline void add(Int &sum, Int &x)
  {
    sum += x;
    if (sum >= x) return; // No overflow case.
    x += ~sum;
    sum = std::numeric_limits<Int>::max();
  }

  inline void sub(Int &sum, Int x) { sum = sum >= x ? sum - x : 0; }

  // Assuming input is in [0.0, 1.0].
  Sample process(Sample input)
  {
    if (input <= std::numeric_limits<Sample>::epsilon()) input = 0;

    Int sum0 = Int(scale * input);

    // Splitting subtraction to avoid unsigned negative overflow.
    add(sum1, sum0);
    auto d1 = delay1.process(sum0);
    sub(sum1, d1);

    add(sum2, sum1);
    auto d2 = delay2.process(sum1);
    sub(sum2, d2);

    auto output = buf == std::numeric_limits<Int>::max()
      ? Sample(0)
      : Sample(buf) * denom - std::numeric_limits<Sample>::epsilon();
    buf = sum2;
    return std::max(Sample(0), output);
  }

  // Assuming input is in [0.0, 1.0].
  Sample processNaive(Sample input)
  {
    if (input <= std::numeric_limits<Sample>::epsilon()) input = 0;

    Int sum0 = Int(scale * input);

    // Splitting subtraction to avoid unsigned negative overflow.
    sum1 += sum0;
    sum1 -= delay1.process(sum0);

    sum2 += sum1;
    sum2 -= delay2.process(sum1);

    auto output = Sample(buf) * denom - std::numeric_limits<Sample>::epsilon();
    buf = sum2;
    return std::max(Sample(0), output);
  }
};

template<typename Sample> struct Limiter {
  static constexpr Sample fixedGain = Sample(0.9965520801347684); // -0.03dB.
  static constexpr Sample releaseConstant = Sample(1e-5); // Small number close to 0.

  Sample threshold = Sample(1);
  Sample gain = Sample(1);
  Sample release = 0; // Release increament per frame.
  size_t attackFrames = 0;

  PeakHold<Sample> hold;
  DoubleAverageFilter<double> smoother;
  // DoubleAverageFilterInt<float, uint_fast64_t> smoother;
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
    attackFrames += attackFrames % 2; // Fix for DoubleAverageFilter.
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

  Sample peak = 0;     // debug
  Sample smoothed = 0; // debug
  Sample delayed = 0;  // debug

  Sample process(Sample input)
  {
    peak = hold.process(std::fabs(input));
    auto &&candidate = applyCharacteristicCurve(peak);
    // gain = candidate;
    gain = std::min(gain * release, candidate);

    smoothed = Sample(smoother.process(gain));
    delayed = lookaheadDelay.process(input);

    // return softClip(smoothed * delayed, fixedGain);
    return smoothed * delayed;
  }
};

constexpr float sampleRate = 48000.0f;
namespace fs = std::filesystem;

void checkThreshold(
  const float threshold, std::string stem, const std::vector<float> &wav)
{
  float max = 0;
  for (const auto &value : wav) {
    auto absed = std::fabs(value);
    if (absed > threshold && absed > max) max = absed;
  }

  if (max == 0) return;
  std::cout << "Limiting failed on: " << stem << "\n"
            << std::setprecision(std::numeric_limits<float>::digits10 + 1)
            << "threshold: " << threshold << "\n"
            << "max      : " << max << "\n";
}

void test(fs::path wavPath, fs::path &outDir)
{
  SoundFile sound(wavPath.string());

  Limiter<float> limiter;
  limiter.resize(2 * size_t(sampleRate));
  limiter.reset();
  limiter.prepare(sampleRate, 64.0f / sampleRate, 0.001f, 0.0f, 1.0f);

  std::vector<float> wav(sound.frames);
  std::vector<float> peak(sound.frames);
  std::vector<float> smoothed(sound.frames);
  std::vector<float> delayed(sound.frames);
  for (size_t i = 0; i < wav.size(); ++i) {
    wav[i] = limiter.process(sound.data[i]);
    peak[i] = limiter.peak;
    smoothed[i] = limiter.smoothed;
    delayed[i] = limiter.delayed;
  }
  checkThreshold(limiter.threshold, wavPath.stem().string(), wav);

  std::string filename = (outDir / wavPath.stem()).string();
  writeWave(filename + "_out.wav", wav, size_t(sampleRate));
  writeWave(filename + "_peak.wav", peak, size_t(sampleRate));
  writeWave(filename + "_smoothed.wav", smoothed, size_t(sampleRate));
  writeWave(filename + "_delayed.wav", delayed, size_t(sampleRate));
}

void testHoldAndSmoother(fs::path wavPath, const fs::path &outDir)
{
  SoundFile sound(wavPath.string());

  PeakHold<float> hold;
  DoubleAverageFilter<float> filter;

  hold.setFrames(32);
  filter.setFrames(32);

  std::vector<float> peak(sound.frames);
  std::vector<float> wav(sound.frames);
  for (size_t i = 0; i < wav.size(); ++i) {
    peak[i] = hold.process(sound.data[i]);
    wav[i] = filter.process(peak[i]);
  }

  writeWave((outDir / wavPath.stem()).string() + "_peak.wav", peak, size_t(sampleRate));
  writeWave((outDir / wavPath.stem()).string() + "_smooth.wav", wav, size_t(sampleRate));
}

template<typename Sample> inline Sample testSub(Sample lhs, Sample rhs)
{
  bool &&swapped = lhs < rhs;
  if (swapped) std::swap(lhs, rhs);
  int expL;
  std::frexp(lhs, &expL);
  auto &&cut = std::ldexp(float(1), expL - std::numeric_limits<Sample>::digits + 1);
  auto &&rem = rhs - std::fmod(rhs, cut);
  return swapped ? rem - lhs : lhs - rem;
}

void printSub()
{
  auto eps = std::numeric_limits<float>::epsilon();
  auto f = std::ldexp(float(0x7fffff), -27);
  auto g = std::ldexp(float(0x7fffff), -23);
  auto h = g - std::fmod(g, std::ldexp(float(1), -27));
  auto sub = testSub(f, g);
  std::cout << std::setprecision(std::numeric_limits<float>::digits10 + 1)
            << std::hexfloat //
            << "       f : " << f << "\n"
            << "   f - g : " << f - g << "\n"
            << "   f - h : " << f - h << "\n"
            << "     sub : " << sub << "\n"
            << "   f - i : " << (f - std::ldexp(float(0xf), -27)) - g << "\n";
}

template<typename Sample> inline Sample testAdd(Sample lhs, Sample rhs)
{
  bool &&swapped = lhs < rhs;
  if (swapped) std::swap(lhs, rhs);
  int expL;
  std::frexp(lhs, &expL);
  auto &&cut = std::ldexp(float(1), expL - std::numeric_limits<Sample>::digits);
  auto &&rem = rhs - std::fmod(rhs, cut);
  return lhs + rem;
}

void printAdd()
{
  auto eps = std::numeric_limits<float>::epsilon();
  auto f = std::ldexp(float(0x1), 1);
  auto g = std::ldexp(float(0x7fffff), -27);
  auto h = g - std::fmod(g, std::ldexp(float(1), -22));
  auto add = testAdd(g, f);
  std::cout << std::setprecision(std::numeric_limits<float>::digits10 + 3)
            << std::hexfloat //
            << "       f : " << f << "\n"
            << "       g : " << g << "\n"
            << "       h : " << h << "\n"
            << "     eps : " << f * eps << "\n"
            << "   f + g : " << f + g << "\n"
            << "   f + h : " << f + h << "\n"
            << "     add : " << add << "\n"
            << "   f + i : " << f + std::ldexp(float(0x7fffe0), -27) << "\n";
}

int main()
{
  // printSub();
  // printAdd();
  // std::cout << std::numeric_limits<float>::round_style << "\n"
  //           << std::numeric_limits<double>::round_style << "\n"
  //           << std::numeric_limits<long double>::round_style << "\n";
  // return 0;

  fs::path outDir("output");
  fs::create_directories(outDir);

  testHoldAndSmoother("snd/input.wav", outDir);

  // test("data/sincrack.wav", outDir);
  // test("data/moduloshaper.wav", outDir);
  // test("data/signal.wav", outDir);
  // test("data/oracleengine.wav", outDir);
  // test("data/sustest2.wav", outDir);
  // test("data/pulse.wav", outDir);
  // test("data/pulse2.wav", outDir);

  return 0;
}
