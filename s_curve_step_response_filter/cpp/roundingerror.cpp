#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

union SPFloat {
  float fp;
  int32_t bin;
};

void test1()
{
  SPFloat s = {0.1f};

  SPFloat a = s;
  for (size_t i = 0; i < 23; ++i) a.fp /= 2;

  SPFloat add = {a.fp + s.fp};
  SPFloat sub = {a.fp - s.fp};
  SPFloat s1 = {s.fp + (a.fp - s.fp)};
  SPFloat s2 = {1.0f - (s.fp + ((1.0f - a.fp) - s.fp))};

  std::cout << std::hex << std::hexfloat //
            << "a  : " << a.bin << ", " << a.fp << "\n"
            << "s  : " << s.bin << ", " << s.fp << "\n"
            << "s1 : " << s1.bin << ", " << s1.fp << "\n"
            << "s2 : " << s2.bin << ", " << s2.fp << "\n"
            << "a+s: " << add.bin << ", " << add.fp << "\n"
            << "a-s: " << sub.bin << ", " << sub.fp << "\n";
}

template<typename Sample> class Delay {
public:
  std::vector<Sample> buffer;
  int32_t wptr = 0;
  int32_t rptr = 0;

  Delay(Sample size = 65536) : buffer(size) {}

  void resize(uint32_t size)
  {
    buffer.resize(size);
    wptr = 0;
    rptr = 0;
  }

  void reset() { std::fill(buffer.begin(), buffer.end(), Sample(0)); }

  void setFrames(uint32_t delayFrames)
  {
    if (delayFrames >= buffer.size()) delayFrames = buffer.size();
    rptr = wptr - int32_t(delayFrames);
    if (rptr < 0) rptr += int32_t(buffer.size());
  }

  Sample process(Sample input)
  {
    wptr = (wptr + 1) % int32_t(buffer.size());
    buffer[wptr] = input;

    rptr = (rptr + 1) % int32_t(buffer.size());
    return buffer[rptr];
  }
};

template<typename Sample> class MovingAverage {
public:
  Sample sum = 0;
  SPFloat out{0.0f};
  uint32_t size = 0;
  Delay<Sample> delay;

  MovingAverage(uint32_t size)
  {
    reset();
    setSize(size);
  }

  void reset()
  {
    sum = 0;
    delay.reset();
  }

  void setSize(uint32_t size)
  {
    delay.setFrames(size);
    this->size = size;
  }

  void process1(Sample input)
  {
    sum += input - delay.process(input);
    out.fp = sum / size;
  }

  void process2(Sample input)
  {
    input = Sample(1) - input;
    sum += input - delay.process(input);
    out.fp = Sample(1) - sum / size;
  }
};

void test2()
{
  constexpr size_t nLoop = 400;
  constexpr uint32_t size = 100;

  MovingAverage<float> av1(size);
  MovingAverage<float> av2(size);
  // av2.sum = av2.size;

  SPFloat big{0.7f};
  SPFloat small{big};
  for (size_t i = 0; i < 23; ++i) small.fp /= 2;

  for (size_t i = 0; i < nLoop; ++i) {
    float input = (i / size) % 2 == 0 ? big.fp : small.fp;
    av1.process1(input);
    av2.process2(input);
  }

  SPFloat sum{(big.fp * (nLoop % size) + small.fp * (size - nLoop % size)) / size};

  std::cout << std::hex << std::hexfloat //
            << "av1: " << av1.out.bin << ", " << av1.out.fp << "\n"
            << "av2: " << av2.out.bin << ", " << av2.out.fp << "\n"
            << "sum: " << sum.bin << ", " << sum.fp << "\n"
            << "big: " << big.bin << ", " << big.fp << "\n"
            << "sml: " << small.bin << ", " << small.fp << "\n";
}

int main()
{
  // test1();
  test2();
  return 0;
}
