#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

constexpr double pi = 3.14159265358979323846264338;
constexpr double twopi = 2.0 * pi;

template<typename Sample, uint8_t size> struct BufferForLoop {
  std::array<Sample, size> buf{};
  std::string name;

  BufferForLoop() : name("For" + std::to_string(size)) {}

  Sample process(Sample input)
  {
    for (uint8_t i = 0; i < size - 1; ++i) buf[i] = buf[i + 1];
    buf.back() = input;
    return buf[0];
  }
};

template<typename Sample, uint8_t size> struct BufferRotate {
  std::array<Sample, size> buf{};
  std::string name;

  BufferRotate() : name("Rot" + std::to_string(size)) {}

  Sample process(Sample input)
  {
    std::rotate(buf.begin(), buf.begin() + 1, buf.end());
    buf.back() = input;
    return buf[0];
  }
};

template<typename Sample>
std::vector<Sample> generateSin(size_t samples, Sample samplerate, Sample frequency)
{
  std::vector<Sample> buf;
  buf.reserve(samples);

  Sample delta = frequency / samplerate;
  Sample phase = 0;
  for (size_t i = 0; i < samples; ++i) {
    buf.push_back(std::sin(Sample(twopi) * phase));
    phase += delta;
    phase -= std::floor(phase);
  }

  return buf;
}

template<typename Buf, typename Sample> void bench(std::vector<Sample> &signal)
{
  Buf buf;

  double sumElapsed = 0.0;
  float sumSignal = 0.0f;
  for (size_t n = 0; n < nLoop; ++n) {
    for (size_t i = 0; i < signal.size(); ++i) {
      auto start = std::chrono::steady_clock::now();
      sumSignal += buf.process(signal[i]);
      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end - start;
      sumElapsed += elapsed.count();
    }
  }

  std::cout << buf.name << ": " << nLoop << " * " << signal.size() << "[sample], "
            << sumElapsed << "[ms], " << sumSignal << "\n";
}

constexpr size_t samplerate = 48000;
constexpr size_t nLoop = 32;

int main()
{
  auto signal = generateSin<float>(samplerate, samplerate, 1000.0f);

  bench<BufferForLoop<float, 8>>(signal);
  bench<BufferForLoop<float, 16>>(signal);
  bench<BufferForLoop<float, 24>>(signal);
  bench<BufferForLoop<float, 32>>(signal);
  bench<BufferForLoop<float, 48>>(signal);
  bench<BufferForLoop<float, 64>>(signal);

  bench<BufferRotate<float, 8>>(signal);
  bench<BufferRotate<float, 16>>(signal);
  bench<BufferRotate<float, 24>>(signal);
  bench<BufferRotate<float, 32>>(signal);
  bench<BufferRotate<float, 48>>(signal);
  bench<BufferRotate<float, 64>>(signal);

  return 0;
}
