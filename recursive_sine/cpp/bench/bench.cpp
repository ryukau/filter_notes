#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <sndfile.h>
#include <string>
#include <vector>

constexpr double pi = 3.14159265358979323846264338;
constexpr double twopi = 2.0 * pi;

void writeWave(std::string filename, std::vector<float> &buffer, const size_t samplerate)
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

// Returns sin.
template<typename T> inline T poly8sincos(T theta, T *cos = nullptr)
{
  constexpr std::array<T, 9> a{
    T(0.7071067811865476),    T(-1.1107203012220528),   T(-0.872357956633633),
    T(0.4567586508845471),    T(0.1793703347937857),    T(-0.05630944123858825),
    T(-0.014746129332847833), T(0.0031986290126005882), T(0.0006324014272427058),
  };

  T normalized = theta * (T(2) / std::numbers::pi_v<T>);
  T floored = std::floor(normalized);
  T x = normalized - floored - T(0.5);

  T A = a[0] + x * x * (a[2] + x * x * (a[4] + x * x * (a[6] + x * x * a[8])));
  T B = x * (a[1] + x * x * (a[3] + x * x * (a[5] + x * x * (a[7]))));

  if (floored < T(1)) {
    if (cos != nullptr) *cos = A + B;
    return A - B;
  } else if (floored < T(2)) {
    if (cos != nullptr) *cos = -A + B;
    return A + B;
  } else if (floored < T(3)) {
    if (cos != nullptr) *cos = -A - B;
    return -A + B;
  }
  if (cos != nullptr) *cos = A - B;
  return -A - B;
}

template<typename T> inline T polysin(T theta)
{
  constexpr std::array<T, 5> a{
    T(1.5707961131825185),    T(-0.645961952025485),     T(0.0796846526139306),
    T(-0.004667901881470699), T(0.00014908811048162785),
  };

  T normalized = theta * (T(2) / std::numbers::pi_v<T>);
  T floored = std::floor(normalized);
  T x = normalized - floored;

#define POLY_SIN                                                                         \
  (x * (a[0] + x * x * (a[1] + x * x * (a[2] + x * x * (a[3] + x * x * a[4])))))

  if (floored < T(1)) {
    return POLY_SIN;
  } else if (floored < T(2)) {
    x = T(1) - x;
    return POLY_SIN;
  } else if (floored < T(3)) {
    return -POLY_SIN;
  }
  x = T(1) - x;
  return -POLY_SIN;

#undef POLY_SIN
}

template<typename T> inline T polycos(T theta)
{
  constexpr std::array<T, 6> a{
    T(1.0),
    T(-1.2337005444481282),
    T(0.25366942754733574),
    T(-0.02086304575763477),
    T(0.0009181995553362909),
    T(-2.4036897062037505e-05),
  };

  T normalized = theta * (T(2) / std::numbers::pi_v<T>);
  T floored = std::floor(normalized);
  T x = normalized - floored;

#define POLY_COS                                                                         \
  (a[0]                                                                                  \
   + x * x * (a[1] + x * x * (a[2] + x * x * (a[3] + x * x * (a[4] + x * x * a[5])))))

  if (floored < T(1)) {
    return POLY_COS;
  } else if (floored < T(2)) {
    x = T(1) - x;
    return -POLY_COS;
  } else if (floored < T(3)) {
    return -POLY_COS;
  }
  x = T(1) - x;
  return POLY_COS;

#undef POLY_COS
}

// `x` in [0, 1]. Theta = 2*pi*x.
template<typename T> inline T normalizedSin(T x)
{
  constexpr std::array<T, 13> a{
    T(6.283185288183828),    T(7.354613470558787e-05), T(-41.344003959896874),
    T(0.031852327445557486), T(81.34625938441181),     T(1.383126951297819),
    T(-81.83560590105544),   T(13.594671970111168),    T(16.042944105268674),
    T(35.73859930409584),    T(-49.38609862670521),    T(21.444085727807128),
    T(-3.299090108871889),
  };
  return
    x * (a[0] + //
    x * (a[1] + //
    x * (a[2] + //
    x * (a[3] + //
    x * (a[4] + //
    x * (a[5] + //
    x * (a[6] + //
    x * (a[7] + //
    x * (a[8] + //
    x * (a[9] + //
    x * (a[10] + //
    x * (a[11] + //
    x * (a[12])))))))))))));

  // constexpr std::array<T, 9> a{
  //   T(6.2831853071147385),  T(0.053068356589663544), T(-42.31874649410266),
  //   T(7.6665725618280485),  T(48.04106433534738),    T(89.84725663966985),
  //   T(-227.03485903958665), T(151.02316071405377),   T(-33.560702380892195),
  // };
  // return
  //   x * (a[0] + //
  //   x * (a[1] + //
  //   x * (a[2] + //
  //   x * (a[3] + //
  //   x * (a[4] + //
  //   x * (a[5] + //
  //   x * (a[6] + //
  //   x * (a[7] + //
  //   x * (a[8])))))))));
}

// `x` in [0, 1]. Theta = 2*pi*x.
template<typename T> inline T normalizedCos(T x)
{
  constexpr std::array<T, 14> a{
    T(0.999999983432531),    T(0.0000000000000),        T(-19.73955306568543),
    T(0.010070968739952127), T(64.81072963064766),      T(0.9513378079994359),
    T(-89.99919246648228),   T(14.76798999895295),      T(26.758717564851263),
    T(53.086341276427056),   T(-84.00408633950707),     T(40.029173640453855),
    T(-6.671529040037597),   T(2.8288922735269303e-08),
  };
  return a[0] +
    x * (a[1] + //
    x * (a[2] + //
    x * (a[3] + //
    x * (a[4] + //
    x * (a[5] + //
    x * (a[6] + //
    x * (a[7] + //
    x * (a[8] + //
    x * (a[9] + //
    x * (a[10] + //
    x * (a[11] + //
    x * (a[12] + //
    x * (a[13])))))))))))));
}

template<typename Sample> struct StdSin {
  Sample invFs = Sample(1) / Sample(44100);
  Sample phase = 0;

  void set(Sample sampleRate) { invFs = Sample(1) / sampleRate; }

  Sample process(Sample freq)
  {
    phase += freq * invFs;
    phase -= std::floor(phase);
    return std::sin(Sample(twopi) * phase);
  }
};

template<typename Sample> struct PolySin {
  Sample invFs = Sample(1) / Sample(44100);
  Sample phase = 0;

  void set(Sample sampleRate) { invFs = Sample(1) / sampleRate; }

  Sample process(Sample freq)
  {
    phase += freq * invFs;
    phase -= std::floor(phase);
    return polysin(Sample(twopi) * phase);
  }
};

template<typename Sample> struct PolyCos {
  Sample invFs = Sample(1) / Sample(44100);
  Sample phase = 0;

  void set(Sample sampleRate) { invFs = Sample(1) / sampleRate; }

  Sample process(Sample freq)
  {
    phase += freq * invFs;
    phase -= std::floor(phase);
    return polycos(Sample(twopi) * phase);
  }
};

template<typename Sample> struct Poly8Sin {
  Sample invFs = Sample(1) / Sample(44100);
  Sample phase = 0;

  void set(Sample sampleRate) { invFs = Sample(1) / sampleRate; }

  Sample process(Sample freq)
  {
    phase += freq * invFs;
    phase -= std::floor(phase);
    return poly8sincos(Sample(twopi) * phase);
  }
};

template<typename Sample> struct NormalizedSin {
  Sample invFs = Sample(1) / Sample(44100);
  Sample phase = 0;

  void set(Sample sampleRate) { invFs = Sample(1) / sampleRate; }

  Sample process(Sample freq)
  {
    phase += freq * invFs;
    phase -= std::floor(phase);
    return normalizedSin(phase);
  }
};

template<typename Sample> struct NormalizedCos {
  Sample invFs = Sample(1) / Sample(44100);
  Sample phase = 0;

  void set(Sample sampleRate) { invFs = Sample(1) / sampleRate; }

  Sample process(Sample freq)
  {
    phase += freq * invFs;
    phase -= std::floor(phase);
    return normalizedCos(phase);
  }
};

template<typename Sample> struct StdSinCos {
  Sample invFs = Sample(1) / Sample(44100);
  Sample phase = 0;

  void set(Sample sampleRate) { invFs = Sample(1) / sampleRate; }

  Sample process(Sample freq)
  {
    phase += freq * invFs;
    phase -= std::floor(phase);
    return std::sin(Sample(twopi) * phase) * std::cos(Sample(twopi) * phase);
  }
};

template<typename Sample> struct Poly8SinCos {
  Sample invFs = Sample(1) / Sample(44100);
  Sample phase = 0;

  void set(Sample sampleRate) { invFs = Sample(1) / sampleRate; }

  Sample process(Sample freq)
  {
    phase += freq * invFs;
    phase -= std::floor(phase);
    Sample cs;
    Sample sn = poly8sincos(Sample(twopi) * phase, &cs);
    return sn * cs;
  }
};

template<typename Sample> struct TableSin {
  static constexpr int32_t tableSize = 0x10'0000;
  std::vector<Sample> table;
  int32_t phase = 0;
  int32_t scale = 1;

  TableSin()
  {
    table.resize(tableSize);
    for (int32_t i = 0; i < table.size(); ++i)
      table[i] = std::sin(Sample(twopi) * Sample(i) / Sample(table.size()));
  }

  void set(int32_t sampleRate) { scale = tableSize / sampleRate; }

  // freq is normalized in [0, 1).
  Sample process(Sample freq)
  {
    int32_t delta = int32_t(freq * scale);
    phase += delta;
    phase &= (tableSize - 1);
    return table[phase];
  }
};

template<typename Sample, size_t nOsc> struct StdSinNote {
  Sample a4freq = Sample(440) / Sample(44100);
  std::array<Sample, nOsc> phase{};

  void set(Sample sampleRate) { a4freq = Sample(440) / sampleRate; }

  Sample process(const std::array<Sample, nOsc> &midinote)
  {
    Sample out = 0;
    for (size_t i = 0; i < nOsc; ++i) {
      auto freq = a4freq * std::exp2((midinote[i] - Sample(69)) / Sample(12));
      phase[i] += freq;
      phase[i] -= std::floor(phase[i]);
      out += std::sin(Sample(twopi) * phase[i]);
    }
    return out;
  }
};

template<typename Sample, size_t nOsc> struct Poly8SinNote {
  Sample sampleRate = 44100;
  std::array<Sample, nOsc> phase{};

  void set(Sample sampleRate) { this->sampleRate = sampleRate; }

  Sample process(const std::array<Sample, nOsc> &midinote)
  {
    Sample out = 0;
    for (size_t i = 0; i < nOsc; ++i) {
      auto freq = Sample(440) * std::exp2((midinote[i] - Sample(69)) / Sample(12));
      phase[i] += freq / sampleRate;
      phase[i] -= std::floor(phase[i]);
      out += poly8sincos(Sample(twopi) * phase[i]);
    }
    return out;
  }
};

template<typename Sample, size_t nOsc> struct NormalizedSinNote {
  Sample a4freq = Sample(440) / Sample(44100);
  std::array<Sample, nOsc> phase{};

  void set(Sample sampleRate) { a4freq = Sample(440) / sampleRate; }

  Sample process(const std::array<Sample, nOsc> &midinote)
  {
    Sample out = 0;
    for (size_t i = 0; i < nOsc; ++i) {
      auto freq = a4freq * std::exp2((midinote[i] - Sample(69)) / Sample(12));
      phase[i] += freq;
      phase[i] -= std::floor(phase[i]);
      out += normalizedSin(phase[i]);
    }
    return out;
  }
};

template<typename Sample, size_t nOsc> struct TableSinNote {
  static constexpr int32_t tableSize = 0x10'0000;
  std::vector<Sample> table;
  std::array<int32_t, nOsc> phase{};
  int32_t scale = 1;

  TableSinNote()
  {
    table.resize(tableSize);
    for (int32_t i = 0; i < table.size(); ++i)
      table[i] = std::sin(Sample(twopi) * Sample(i) / Sample(table.size()));
  }

  void set(int32_t sampleRate) { scale = tableSize / sampleRate; }

  Sample process(const std::array<Sample, nOsc> &midinote)
  {
    Sample out = 0;
    for (size_t i = 0; i < nOsc; ++i) {
      auto freq = Sample(440) * std::exp2((midinote[i] - Sample(69)) / Sample(12));
      int32_t delta = int32_t(freq * scale);
      phase[i] += delta;
      phase[i] &= (tableSize - 1);
      out += table[phase[i]];
    }
    return out;
  }
};

template<typename Sample, size_t nOsc> class MagicCircleNote {
public:
  static constexpr size_t centPerOctave = 1200;
  static constexpr size_t centPerSemitone = 100;
  static constexpr size_t nSemitone = 156;
  static constexpr auto lowestFreq = Sample(2.575215288381796); // 440 * 2^(-89 / 12).

  std::array<Sample, nSemitone * centPerSemitone> table{};
  std::array<Sample, nOsc> u{};
  std::array<Sample, nOsc> v{};

  void set(Sample sampleRate)
  {
    auto normalizedF0 = lowestFreq / sampleRate;
    for (size_t i = 0; i < table.size(); ++i) {
      auto freq = normalizedF0 * std::exp2(Sample(i) / Sample(centPerOctave));
      table[i] = 2 * std::sin(Sample(pi) * freq);
    }

    u.fill(Sample(1));
    v.fill(Sample(0));
  }

  // Range of `midinote` is [-20, 136].
  Sample process(const std::array<Sample, nOsc> &midinote)
  {
    for (size_t i = 0; i < nOsc; ++i) {
      auto index
        = std::min(size_t(centPerSemitone * (midinote[i] + 20)), table.size() - 1);
      u[i] -= table[index] * v[i];
      v[i] += table[index] * u[i];
    }
    return std::accumulate(v.begin(), v.end(), Sample(0));
  }
};

template<typename Osc> void benchFreq(std::string name)
{
  constexpr size_t nLoop = 1;
  constexpr size_t sampleRate = 48000;

  Osc osc;
  osc.set(sampleRate);

  double sumElapsed = 0.0;
  std::vector<float> wav(sampleRate);
  for (size_t n = 0; n < nLoop; ++n) {
    for (size_t i = 0; i < wav.size(); ++i) {
      auto freq = 100.0f + 50.0f * std::sin(8 * i / float(sampleRate));

      auto start = std::chrono::steady_clock::now();

      wav[i] = osc.process(freq);

      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end - start;
      sumElapsed += elapsed.count();
    }
  }
  std::cout << name << ": " << wav.size() << "[sample], " << sumElapsed << "[ms], "
            << "\n";

  auto path = "snd/" + name + ".wav";
  writeWave(path, wav, sampleRate);
}

template<typename Osc, size_t nOsc> void benchMidinote(std::string name)
{
  constexpr size_t nLoop = 1;
  constexpr size_t sampleRate = 48000;

  Osc osc;
  osc.set(sampleRate);

  double sumElapsed = 0.0;
  std::vector<float> wav(sampleRate);
  std::array<float, nOsc> midinote;
  for (size_t n = 0; n < nLoop; ++n) {
    for (size_t i = 0; i < wav.size(); ++i) {
      midinote[0] = 50.0f + 50.0f * std::sin(8 * i / float(sampleRate));
      for (size_t i = 1; i < midinote.size(); ++i) midinote[i] = midinote[i - 1] - 1.0f;

      auto start = std::chrono::steady_clock::now();

      wav[i] = osc.process(midinote) / nOsc;

      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end - start;
      sumElapsed += elapsed.count();
    }
  }
  std::cout << name << ": " << wav.size() << "[sample], " << sumElapsed << "[ms], "
            << "\n";

  auto path = "snd/" + name + ".wav";
  writeWave(path, wav, sampleRate);
}

int main()
{
  std::cout << "---sin or cos\n";
  benchFreq<TableSin<float>>("TableSin");
  benchFreq<StdSin<float>>("StdSin");
  benchFreq<PolyCos<float>>("PolyCos");
  benchFreq<PolySin<float>>("PolySin");
  benchFreq<Poly8Sin<float>>("Poly8Sin");
  benchFreq<NormalizedSin<float>>("NormalizedSin");
  benchFreq<NormalizedCos<float>>("NormalizedCos");

  std::cout << "\n---sincos\n";
  benchFreq<StdSinCos<float>>("StdSinCos");
  benchFreq<Poly8SinCos<float>>("Poly8SinCos");

  std::cout << "\n--- Multi Output\n";
  constexpr size_t nOsc = 8;
  benchMidinote<TableSinNote<float, nOsc>, nOsc>("TableSinNote");
  benchMidinote<StdSinNote<float, nOsc>, nOsc>("StdSinNote");
  benchMidinote<Poly8SinNote<float, nOsc>, nOsc>("Poly8SinNote");
  benchMidinote<NormalizedSinNote<float, nOsc>, nOsc>("NormalizedSinNote");
  benchMidinote<MagicCircleNote<float, nOsc>, nOsc>("MagicCircleNote");

  return 0;
}
