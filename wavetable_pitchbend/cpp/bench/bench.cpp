#define POCKETFFT_NO_MULTITHREADING
#include "pocketfft_hdronly.h"

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

template<typename T> class PocketFFT {
public:
  pocketfft::shape_t shape;
  pocketfft::stride_t strideR;
  pocketfft::stride_t strideC;
  pocketfft::shape_t axes;
  size_t ndata = 1;

  void setShape(pocketfft::shape_t shape)
  {
    this->shape = shape;

    strideR.resize(shape.size());
    strideC.resize(shape.size());

    size_t tmpR = sizeof(T);
    size_t tmpC = sizeof(std::complex<T>);
    for (int i = int(shape.size()) - 1; i >= 0; --i) {
      strideR[i] = tmpR;
      tmpR *= shape[i];
      strideC[i] = tmpC;
      tmpC *= shape[i];
    }

    ndata = 1;
    for (const auto &shp : shape) ndata *= shp;

    axes.resize(shape.size());
    for (size_t i = 0; i < axes.size(); ++i) axes[i] = i;
  }

  void r2c(const T *data_in, std::complex<T> *data_out, bool forward = true, T scale = 1)
  {
    pocketfft::r2c(shape, strideR, strideC, axes, forward, data_in, data_out, scale);
  }

  void c2r(const std::complex<T> *data_in, T *data_out, bool forward = false, T scale = 1)
  {
    pocketfft::c2r(
      shape, strideC, strideR, axes, forward, data_in, data_out, scale / ndata);
  }
};

/*
---python
import numpy
from scipy import signal
sos = signal.ellip(12, 0.01, 100, 0.4, "low", output="sos", fs=2)
---
*/
template<typename Sample> struct SosEllipticDecimation2 {
  constexpr static size_t nSection = 6;

  constexpr static std::array<std::array<Sample, 5>, nSection> co{{
    {Sample(0.0010527785162220168), Sample(0.0019568431992735415),
     Sample(0.0010527785162220173), Sample(-1.0423808853583867),
     Sample(0.2987544208875724)},
    {Sample(1.0), Sample(1.0976312338288472), Sample(0.9999999999999999),
     Sample(-0.9342188057739501), Sample(0.4275886413822211)},
    {Sample(1.0), Sample(0.4245180333045316), Sample(0.9999999999999998),
     Sample(-0.7883710918539022), Sample(0.6047059346890982)},
    {Sample(1.0), Sample(0.036612291733545285), Sample(1.0), Sample(-0.6673897039419319),
     Sample(0.7605822883546196)},
    {Sample(1.0), Sample(-0.15846910443551615), Sample(0.9999999999999999),
     Sample(-0.5934855881848028), Sample(0.8756817382088157)},
    {Sample(1.0), Sample(-0.23876536287134129), Sample(1.0000000000000002),
     Sample(-0.5673350256436269), Sample(0.9613812181636632)},
  }};
};

// SOS: Second order sections.
template<typename Sample, typename IIR, size_t oversample = 2> class SosFilter {
public:
  void reset()
  {
    x0.fill(0);
    x1.fill(0);
    x2.fill(0);
    y0.fill(0);
    y1.fill(0);
    y2.fill(0);
  }

  void push(Sample input)
  {
    x0[0] = input;
    for (size_t i = 1; i < IIR::nSection; ++i) x0[i] = y0[i - 1];

    for (size_t i = 0; i < IIR::nSection; ++i) {
      y0[i] = +IIR::co[i][0] * x0[i] + IIR::co[i][1] * x1[i] + IIR::co[i][2] * x2[i]
        - IIR::co[i][3] * y1[i] - IIR::co[i][4] * y2[i];
    }

    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }

  inline Sample output() { return y0[IIR::nSection - 1]; }

  Sample process(const std::array<Sample, oversample> &input)
  {
    for (const auto &value : input) push(value);
    return output();
  }

  std::array<Sample, IIR::nSection> x0{};
  std::array<Sample, IIR::nSection> x1{};
  std::array<Sample, IIR::nSection> x2{};
  std::array<Sample, IIR::nSection> y0{};
  std::array<Sample, IIR::nSection> y1{};
  std::array<Sample, IIR::nSection> y2{};
};

template<typename Sample> inline Sample frequencyToMidinote(Sample freq)
{
  return Sample(12) * std::log2(freq / Sample(440)) + Sample(69);
}

template<typename Sample> inline Sample midinoteToFrequency(Sample note)
{
  return Sample(440) * std::exp2((note - Sample(69)) / Sample(12));
}

// Range of t is in [0, 1]. Interpoltes between y1 and y2.
template<typename Sample>
inline Sample cubicInterp(Sample y0, Sample y1, Sample y2, Sample y3, Sample t)
{
  auto t2 = t * t;
  auto c0 = y1 - y2;
  auto c1 = (y2 - y0) * Sample(0.5);
  auto c2 = c0 + c1;
  auto c3 = c0 + c2 + (y3 - y1) * Sample(0.5);
  return c3 * t * t2 - (c2 + c3) * t2 + c1 * t + y1;
}

template<typename Sample>
std::vector<std::complex<Sample>> generateSawSpectrum(size_t size)
{
  std::vector<std::complex<Sample>> spectrum(size / 2 + 1);

  spectrum[0] = std::complex<Sample>(0, 0);

  Sample sign = Sample(-1);
  Sample scale = -Sample(size) / Sample(pi);
  for (size_t k = 1; k < spectrum.size(); ++k) {
    spectrum[k] = std::complex<Sample>(0, scale * sign / k);
    sign *= Sample(-1);
  }

  return spectrum;
}

template<typename Sample> class TableOsc {
public:
  Sample sampleRate = Sample(44100);
  Sample phase = Sample(0);
  Sample basenote = Sample(0);
  size_t size = 0;
  SosFilter<Sample, SosEllipticDecimation2<float>, 2> lowpass;
  std::vector<std::vector<Sample>> table;
  PocketFFT<Sample> fft;

  TableOsc(Sample sampleRate)
  {
    this->sampleRate = Sample(2) * sampleRate;

    // Numerator in log2 is added for experiment.
    // It is the lowest frequency that the wavetable can play without loss of harmonics.
    // 1 Hz in this case.
    auto exponent = size_t(std::floor(std::log2(this->sampleRate / Sample(1))));
    exponent = std::clamp<size_t>(exponent, 1, std::numeric_limits<size_t>::digits - 1);
    size = size_t(1) << exponent;
    auto spectrum = generateSawSpectrum<Sample>(size);

    pocketfft::shape_t shape{size};
    fft.setShape(shape);

    basenote = frequencyToMidinote(this->sampleRate / Sample(2 * size));

    // On cutoff calculation:
    // - `idx + 1` is the same as using `tableSize / 2` as numerator.
    // - Last +1 for DC component.
    table.resize(exponent);
    std::vector<std::complex<Sample>> spec(spectrum.size());
    for (size_t idx = 0; idx < table.size(); ++idx) {
      std::fill(spec.begin(), spec.end(), std::complex<Sample>(0, 0));
      auto cutoff = size / (size_t(1) << (idx + 1)) + 1; // +1 for DC component.
      std::copy(spectrum.begin(), spectrum.begin() + cutoff, spec.begin());

      table[idx].resize(size + 1);
      fft.c2r(spec.data(), table[idx].data());
      table[idx].back() = table[idx][0];
    }
  }

  void reset()
  {
    phase = Sample(0);
    lowpass.reset();
  }

  void debugRenderTable()
  {
    for (size_t idx = 0; idx < table.size(); ++idx) {
      std::string name = "snd/simple" + std::to_string(idx) + ".wav";
      writeWave(name, table[idx], size_t(sampleRate));
    }
  }

  Sample processSample(Sample note)
  {
    phase += std::clamp(midinoteToFrequency(note) / sampleRate, Sample(0), Sample(0.5));
    phase -= std::floor(phase);

    auto pos = Sample(size) * phase;
    auto idx = size_t(pos);
    auto frac = pos - Sample(idx);

    auto octave = std::clamp(size_t(note - basenote) / 12, size_t(0), table.size() - 1);
    auto x0 = table[octave][idx];
    auto x1 = table[octave][idx + 1];
    return x0 + frac * (x1 - x0);
  }

  Sample process(Sample note)
  {
    std::array<Sample, 2> sample;
    for (auto &value : sample) value = processSample(note);
    return lowpass.process(sample);
  }
};

template<typename Sample> class ArrayOsc {
public:
  static constexpr size_t tableSize = 8192;
  static constexpr size_t nOctave = 16;

  Sample sampleRate = Sample(44100);
  Sample phase = Sample(0);
  Sample basenote = Sample(0);
  SosFilter<Sample, SosEllipticDecimation2<float>, 2> lowpass;
  std::array<std::array<Sample, tableSize + 1>, nOctave> table;
  PocketFFT<Sample> fft;

  ArrayOsc(Sample sampleRate)
  {
    this->sampleRate = Sample(2) * sampleRate;

    auto spectrum = generateSawSpectrum<Sample>(tableSize);

    pocketfft::shape_t shape{tableSize};
    fft.setShape(shape);

    basenote = frequencyToMidinote(this->sampleRate / Sample(2 * tableSize));

    std::vector<std::complex<Sample>> spec(spectrum.size());
    for (size_t idx = 0; idx < nOctave; ++idx) {
      std::fill(spec.begin(), spec.end(), std::complex<Sample>(0, 0));
      auto cutoff = tableSize / (size_t(1) << (idx + 1)) + 1;
      std::copy(spectrum.begin(), spectrum.begin() + cutoff, spec.begin());

      fft.c2r(spec.data(), table[idx].data());
      table[idx].back() = table[idx][0];
    }
  }

  void reset()
  {
    phase = Sample(0);
    lowpass.reset();
  }

  void debugRenderTable()
  {
    std::vector<Sample> wav(table[0].size());
    for (size_t idx = 0; idx < nOctave; ++idx) {
      std::string name = "snd/array" + std::to_string(idx) + ".wav";
      std::copy(table[idx].begin(), table[idx].end(), wav.begin());
      writeWave(name, wav, size_t(sampleRate));
    }
  }

  Sample processSample(Sample note)
  {
    phase += std::clamp(midinoteToFrequency(note) / sampleRate, Sample(0), Sample(0.5));
    phase -= std::floor(phase);

    auto pos = Sample(tableSize) * phase;
    auto idx = size_t(pos);
    auto frac = pos - Sample(idx);

    auto octave = std::clamp(size_t(note - basenote) / 12, size_t(0), nOctave - 1);
    auto x0 = table[octave][idx];
    auto x1 = table[octave][idx + 1];
    return x0 + frac * (x1 - x0);
  }

  Sample process(Sample note)
  {
    std::array<Sample, 2> sample;
    for (auto &value : sample) value = processSample(note);
    return lowpass.process(sample);
  }
};

template<typename Sample> class TableOscAltInterval {
public:
  Sample sampleRate = Sample(44100);
  Sample phase = Sample(0);
  Sample basenote = Sample(0);
  Sample interval = Sample(12);
  Sample maxIdx = Sample(0);
  size_t size = 0;
  SosFilter<Sample, SosEllipticDecimation2<float>, 2> lowpass;
  std::vector<std::vector<Sample>> table;
  PocketFFT<Sample> fft;

  TableOscAltInterval(Sample sampleRate)
  {
    this->sampleRate = Sample(2) * sampleRate;

    auto exponent = size_t(std::log2(this->sampleRate / Sample(10)));
    exponent = std::clamp<size_t>(exponent, 1, std::numeric_limits<size_t>::digits - 1);
    size = size_t(1) << exponent;
    auto specSize = size / 2;
    auto spectrum = generateSawSpectrum<Sample>(size);

    pocketfft::shape_t shape{size};
    fft.setShape(shape);

    const Sample bendRange = Sample(1.5);
    const size_t nTable
      = size_t(-std::log(Sample(1) / Sample(specSize)) / std::log(bendRange));
    maxIdx = Sample(nTable - 1);

    basenote = frequencyToMidinote(this->sampleRate / Sample(size));
    interval = Sample(12) * std::log2(bendRange);

    table.resize(nTable + 1); // Last table is filled by 0.
    std::vector<std::complex<Sample>> spec(spectrum.size());
    for (size_t idx = 0; idx < table.size() - 1; ++idx) {
      std::fill(spec.begin(), spec.end(), std::complex<Sample>(0, 0));
      auto cutoff = size_t(specSize * std::pow(bendRange, -Sample(idx))) + 1;
      std::copy(spectrum.begin(), spectrum.begin() + cutoff, spec.begin());

      table[idx].resize(size + 1);
      fft.c2r(spec.data(), table[idx].data());
      table[idx].back() = table[idx][0];
    }
    table.back().resize(size + 1);
    std::fill(table.back().begin(), table.back().end(), Sample(0));
  }

  void reset()
  {
    phase = Sample(0);
    lowpass.reset();
  }

  void debugRenderTable()
  {
    for (size_t idx = 0; idx < table.size(); ++idx) {
      std::string name = "snd/altinterval" + std::to_string(idx) + ".wav";
      writeWave(name, table[idx], size_t(sampleRate));
    }
  }

  Sample processSample(Sample note)
  {
    phase += std::clamp(midinoteToFrequency(note) / sampleRate, Sample(0), Sample(0.5));
    phase -= std::floor(phase);

    auto octFloat = std::clamp((note - basenote) / interval, Sample(0), maxIdx);
    auto iTbl = size_t(octFloat);
    auto yFrac = octFloat - Sample(iTbl);

    auto pos = Sample(size) * phase;
    auto idx = size_t(pos);
    auto xFrac = pos - Sample(idx);

    auto a0 = table[iTbl][idx];
    auto a1 = table[iTbl][idx + 1];
    auto s0 = a0 + xFrac * (a1 - a0);

    auto b0 = table[iTbl + 1][idx];
    auto b1 = table[iTbl + 1][idx + 1];
    auto s1 = b0 + xFrac * (b1 - b0);

    return s0 + yFrac * (s1 - s0);
  }

  Sample process(Sample note)
  {
    std::array<Sample, 2> sample;
    for (auto &value : sample) value = processSample(note);
    return lowpass.process(sample);
  }
};

template<typename Sample> class MipmapOsc {
public:
  Sample sampleRate = Sample(44100);
  Sample phase = Sample(0);
  Sample basenote = Sample(0);
  size_t size = 0;
  SosFilter<Sample, SosEllipticDecimation2<float>, 2> lowpass;
  std::vector<std::vector<Sample>> table;
  PocketFFT<Sample> fft;

  MipmapOsc(Sample sampleRate)
  {
    this->sampleRate = Sample(2) * sampleRate;

    auto exponent = size_t(std::floor(std::log2(this->sampleRate / Sample(1))));
    if (exponent >= std::numeric_limits<size_t>::digits)
      exponent = std::numeric_limits<size_t>::digits - 1;
    size = size_t(1) << exponent;
    auto spectrum = generateSawSpectrum<Sample>(size);

    basenote = frequencyToMidinote(this->sampleRate / Sample(size));

    if (exponent <= 2) exponent = 3;
    table.resize(size_t(exponent) - 2);
    std::vector<std::complex<Sample>> spec;
    for (size_t idx = 0; idx < table.size(); ++idx) {
      auto cutoff = size / (size_t(1) << (idx + 2)) + 1;

      spec.resize(cutoff, std::complex<Sample>(0, 0));
      for (size_t j = 0; j < spec.size(); ++j)
        spec[j] = spectrum[j] / Sample(size_t(1) << (idx + 1));

      table[idx].resize(2 * spec.size() - 1);

      pocketfft::shape_t shape{table[idx].size() - 1};
      fft.setShape(shape);
      fft.c2r(spec.data(), table[idx].data());

      table[idx].back() = table[idx][0];
    }
  }

  void reset()
  {
    phase = Sample(0);
    lowpass.reset();
  }

  void debugRenderTable()
  {
    for (size_t idx = 0; idx < table.size(); ++idx) {
      std::string name = "snd/mipmap" + std::to_string(idx) + ".wav";
      writeWave(name, table[idx], size_t(sampleRate));
    }
  }

  Sample processSample(Sample note)
  {
    phase += std::clamp(midinoteToFrequency(note) / sampleRate, Sample(0), Sample(0.5));
    phase -= std::floor(phase);

    auto octave = std::clamp(size_t(note - basenote) / 12, size_t(0), table.size() - 1);
    auto pos = Sample(table[octave].size() - 1) * phase;
    auto idx = size_t(pos);
    auto frac = pos - Sample(idx);

    auto x0 = table[octave][idx];
    auto x1 = table[octave][idx + 1];
    return x0 + frac * (x1 - x0);
  }

  Sample process(Sample note)
  {
    std::array<Sample, 2> sample;
    for (auto &value : sample) value = processSample(note);
    return lowpass.process(sample);
  }
};

template<typename Sample> class LpsOsc {
public:
  Sample sampleRate = Sample(44100);
  Sample phase = Sample(0);
  Sample basenote = Sample(0);
  Sample tableSize = Sample(0);
  size_t maxNote = 0;
  std::vector<std::vector<Sample>> table;
  PocketFFT<Sample> fft;

  LpsOsc(Sample sampleRate)
  {
    this->sampleRate = sampleRate;

    auto minFreq = midinoteToFrequency<Sample>(0);
    auto expFloat = std::log2(this->sampleRate / std::floor(minFreq));
    auto exponent = size_t(std::ceil(expFloat));
    if (exponent >= std::numeric_limits<size_t>::digits)
      exponent = std::numeric_limits<size_t>::digits - 1;
    auto sizeInt = size_t(1) << exponent;
    tableSize = Sample(sizeInt);
    auto spectrum = generateSawSpectrum<Sample>(sizeInt);

    pocketfft::shape_t shape{sizeInt};
    fft.setShape(shape);

    maxNote = size_t(std::ceil(frequencyToMidinote<Sample>(20000)));
    table.resize(maxNote + 1);
    std::vector<std::complex<Sample>> spec(spectrum.size());
    for (size_t idx = 0; idx < table.size(); ++idx) {
      auto freq = midinoteToFrequency(Sample(idx));
      auto cutoff = size_t(sizeInt / 2 * minFreq / freq);
      std::copy(spectrum.begin(), spectrum.begin() + cutoff, spec.begin());
      std::fill(spec.begin() + cutoff, spec.end(), std::complex<Sample>(0, 0));

      table[idx].resize(sizeInt + 1);
      fft.c2r(spec.data(), table[idx].data());
      table[idx].back() = table[idx][0];
    }
    std::fill(table.back().begin(), table.back().end(), Sample(0));
  }

  void reset() { phase = Sample(0); }

  void debugRenderTable()
  {
    for (size_t idx = 0; idx < table.size(); ++idx) {
      std::string name = "snd/lpsosc" + std::to_string(idx) + ".wav";
      writeWave(name, table[idx], size_t(sampleRate));
    }
  }

  inline void processPhase(Sample note)
  {
    auto tick = midinoteToFrequency(note) * tableSize / sampleRate;
    if (tick >= tableSize || tick < 0) tick = 0;

    phase += tick;
    if (phase >= tableSize) phase -= tableSize;
  }

  Sample processLow(Sample note)
  {
    auto idx = size_t(phase);
    const auto &a0 = table[0][idx];
    const auto &a1 = table[0][idx + 1];
    auto fracX = phase - Sample(size_t(phase));
    return a0 + fracX * (a1 - a0);
  }

  Sample processFull(Sample note)
  {
    processPhase(note);

    if (note < 0) return processLow(note);

    auto nn = size_t(note);
    if (nn > maxNote) nn = maxNote - 1;

    auto idx = size_t(phase);
    const auto &a0 = table[nn][idx];
    const auto &a1 = table[nn][idx + 1];
    const auto &b0 = table[nn + 1][idx];
    const auto &b1 = table[nn + 1][idx + 1];

    auto fracX = phase - Sample(size_t(phase));
    auto x0 = a0 + fracX * (a1 - a0);
    auto x1 = b0 + fracX * (b1 - b0);

    return x0 + (note - Sample(nn)) * (x1 - x0);
  }

  Sample process(Sample note)
  {
    processPhase(note);

    auto nn = std::clamp<size_t>(size_t(note), 0, maxNote);

    auto idx = size_t(phase);
    const auto &a0 = table[nn][idx];
    const auto &a1 = table[nn][idx + 1];
    auto fracX = phase - Sample(size_t(phase));
    return a0 + fracX * (a1 - a0);
  }
};

template<typename Sample> class CpsOsc {
public:
  Sample sampleRate = Sample(44100);
  Sample phase = Sample(0);
  Sample basenote = Sample(0);
  Sample tableSize = Sample(0);
  size_t maxNote = 0;
  std::vector<std::vector<Sample>> table;
  PocketFFT<Sample> fft;

  CpsOsc(Sample sampleRate)
  {
    this->sampleRate = sampleRate;

    auto minFreq = midinoteToFrequency<Sample>(0);
    auto expFloat = std::log2(this->sampleRate / std::floor(minFreq));
    auto exponent = size_t(std::ceil(expFloat));
    if (exponent >= std::numeric_limits<size_t>::digits)
      exponent = std::numeric_limits<size_t>::digits - 1;
    auto sizeInt = size_t(1) << exponent;
    tableSize = Sample(sizeInt);
    auto spectrum = generateSawSpectrum<Sample>(sizeInt);

    pocketfft::shape_t shape{sizeInt};
    fft.setShape(shape);

    maxNote = size_t(std::ceil(frequencyToMidinote<Sample>(20000)));
    table.resize(maxNote + 1);
    std::vector<std::complex<Sample>> spec(spectrum.size());
    size_t back1 = sizeInt + 2;
    size_t back2 = sizeInt + 1;
    size_t back3 = sizeInt;
    for (size_t idx = 0; idx < table.size(); ++idx) {
      auto freq = midinoteToFrequency(Sample(idx));
      auto cutoff = size_t(sizeInt / 2 * minFreq / freq);
      std::copy(spectrum.begin(), spectrum.begin() + cutoff, spec.begin());
      std::fill(spec.begin() + cutoff, spec.end(), std::complex<Sample>(0, 0));

      table[idx].resize(sizeInt + 3);
      fft.c2r(spec.data(), table[idx].data());

      std::rotate(table[idx].rbegin(), table[idx].rbegin() + 1, table[idx].rend());

      table[idx][0] = table[idx][back3];
      table[idx][back2] = table[idx][1];
      table[idx][back1] = table[idx][2];
    }
    std::fill(table.back().begin(), table.back().end(), Sample(0));
  }

  void reset() { phase = Sample(0); }

  void debugRenderTable()
  {
    for (size_t idx = 0; idx < table.size(); ++idx) {
      std::string name = "snd/cpsosc" + std::to_string(idx) + ".wav";
      writeWave(name, table[idx], size_t(sampleRate));
    }
  }

  inline void processPhase(Sample note)
  {
    auto tick = midinoteToFrequency(note) * tableSize / sampleRate;
    if (tick >= tableSize || tick < 0) tick = 0;

    phase += tick;
    if (phase >= tableSize) phase -= tableSize;
  }

  Sample process(Sample note)
  {
    processPhase(note);
    auto nn = note < 0 ? 0 : size_t(note);
    if (nn >= maxNote) nn = maxNote - 1;
    auto idx = size_t(phase);
    return cubicInterp(
      table[nn][idx], table[nn][idx + 1], table[nn][idx + 2], table[nn][idx + 3],
      phase - Sample(size_t(phase)));
  }
};

constexpr size_t nLoop = 1;

template<typename Osc> void bench(std::string name)
{
  constexpr float sampleRate = 48000.0f;

  Osc osc(sampleRate);
  // osc.debugRenderTable();

  double sumElapsed = 0.0;
  std::vector<float> wav;
  wav.resize(8 * size_t(sampleRate));
  for (size_t loop = 0; loop < nLoop; ++loop) {
    osc.reset();
    for (size_t i = 0; i < wav.size(); ++i) {
      auto note = 128.0f * i / float(wav.size());

      auto start = std::chrono::steady_clock::now();

      wav[i] = 0.25f * osc.process(note);

      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end - start;
      sumElapsed += elapsed.count();
    }
  }
  std::cout << name << ": " << wav.size() << "[sample], " << sumElapsed << "[ms], "
            << "\n";

  std::string filename = "snd/chirp_" + name + ".wav";
  writeWave(filename, wav, size_t(sampleRate));
}

template<typename Osc> void testOutOfRangePitch(std::string name)
{
  constexpr float sampleRate = 48000.0f;

  Osc osc(sampleRate);
  osc.reset();

  std::vector<float> wav;
  wav.resize(8 * size_t(sampleRate));
  for (size_t i = 0; i < wav.size(); ++i) {
    auto note = 300.0f * i / float(wav.size()) - 100;
    wav[i] = 0.25f * osc.process(note);
  }
  std::cout << name << ": "
            << "\n";

  std::string filename = "snd/testPitch_" + name + ".wav";
  writeWave(filename, wav, size_t(sampleRate));
}

int main()
{
  std::cout << "--- Warm up\n";
  bench<TableOsc<float>>("Simple");

  std::cout << "\n--- Benchmark\n";
  bench<TableOsc<float>>("Simple");
  bench<ArrayOsc<float>>("Array");
  bench<TableOscAltInterval<float>>("AltInterval");
  bench<MipmapOsc<float>>("Mipmap");
  bench<LpsOsc<float>>("Lpsosc");
  bench<CpsOsc<float>>("Cpsosc");

  std::cout << "\n--- Test out of range pitch\n";
  testOutOfRangePitch<TableOsc<float>>("Simple");
  testOutOfRangePitch<ArrayOsc<float>>("Array");
  testOutOfRangePitch<TableOscAltInterval<float>>("AltInterval");
  testOutOfRangePitch<MipmapOsc<float>>("Mipmap");
  testOutOfRangePitch<LpsOsc<float>>("Lpsosc");
  testOutOfRangePitch<CpsOsc<float>>("Cpsosc");

  return 0;
}
