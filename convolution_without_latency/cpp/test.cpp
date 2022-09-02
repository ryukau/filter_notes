/**
FFTConvolutionFilter: 192000[sample], Average11.7556[ms], Spike 1.2079[ms]
ImmediateConvolver: 192000[sample], Average50.0951[ms], Spike 0.9564[ms]
NaiveSplitConvolver: 192000[sample], Average55.6657[ms], Spike 0.415[ms]
NaiveSplitConvolver: 192000[sample], Average67.6343[ms], Spike 0.1786[ms]
NaiveSplitConvolver: 192000[sample], Average98.417[ms], Spike 0.0882[ms]
NaiveSplitConvolver: 192000[sample], Average148.647[ms], Spike 0.0455[ms]
NaiveSplitConvolver: 192000[sample], Average259.182[ms], Spike 0.0413[ms]
NaiveSplitConvolver: 192000[sample], Average541.798[ms], Spike 0.0366[ms]
*/

#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <numeric>
#include <sndfile.h>
#include <vector>

constexpr double halfpi = 1.57079632679489661923;
constexpr double pi = 3.14159265358979323846;
constexpr double twopi = 6.28318530717958647692;

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
  if (sf_write_float(file, &buffer[0], length) != sf_count_t(length))
    std::cout << sf_strerror(file) << std::endl;

  if (sf_close(file) != 0) {
    std::cout << "Error: sf_close failed." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void writeWave(
  std::string filename, float *buffer, size_t bufferSize, const size_t &samplerate)
{
  SF_INFO sfinfo;
  memset(&sfinfo, 0, sizeof(sfinfo));
  sfinfo.samplerate = int(samplerate);
  sfinfo.frames = bufferSize;
  sfinfo.channels = 1;
  sfinfo.format = (SF_FORMAT_WAV | SF_FORMAT_FLOAT);

  SNDFILE *file = sf_open(filename.c_str(), SFM_WRITE, &sfinfo);
  if (!file) {
    std::cout << "Error: sf_open failed." << std::endl;
    exit(EXIT_FAILURE);
  }

  size_t length = sfinfo.channels * bufferSize;
  if (sf_write_float(file, buffer, length) != sf_count_t(length))
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

template<size_t nTap> class FFTConvolutionFilter {
private:
  static constexpr size_t half = nTap;
  static constexpr size_t bufSize = 2 * half;
  static constexpr size_t spcSize = nTap + 1; // spc = spectrum.

  std::array<float *, 2> buf;
  std::complex<float> *spc;

  float *coefficient;
  std::complex<float> *fir;

  float *flt; // filtered.

  fftwf_plan firPlan;
  std::array<fftwf_plan, 2> forwardPlan;
  fftwf_plan inversePlan;

  size_t front = 0;
  std::array<size_t, 2> wptr{};
  size_t rptr = 0;

  void debugWriteComplex(std::string name, std::complex<float> *arr, size_t size)
  {
    std::vector<float> real(size);
    std::vector<float> imag(size);

    for (size_t i = 0; i < size; ++i) {
      real[i] = arr[i].real();
      imag[i] = arr[i].imag();
    }

    std::string realname = "snd/" + name + "_real.wav";
    std::string imagname = "snd/" + name + "_imag.wav";
    writeWave(realname, real, 48000);
    writeWave(imagname, imag, 48000);
  }

  void debugDumpAll()
  {
    // writeWave("snd/buf.wav", buf, bufSize, 48000);
    // debugWriteComplex("spc", spc, spcSize);
    // writeWave("snd/flt0.wav", flt[0], bufSize, 48000);
    // writeWave("snd/flt1.wav", flt[1], bufSize, 48000);
    writeWave("snd/coefficient.wav", coefficient, bufSize, 48000);
    // debugWriteComplex("fir", fir, spcSize);
  }

public:
  FFTConvolutionFilter()
  {
    static_assert(
      nTap && ((nTap & (nTap - 1)) == 0),
      "FFTConvolutionFilter: nTap must be power of 2.");

    for (auto &bf : buf) bf = (float *)fftwf_malloc(sizeof(float) * bufSize);
    spc = (std::complex<float> *)fftwf_malloc(sizeof(std::complex<float>) * spcSize);
    flt = (float *)fftwf_malloc(sizeof(float) * bufSize);

    coefficient = (float *)fftwf_malloc(sizeof(float) * bufSize);
    std::fill(coefficient, coefficient + bufSize, float(0));

    fir = (std::complex<float> *)fftwf_malloc(sizeof(std::complex<float>) * spcSize);
    firPlan = fftwf_plan_dft_r2c_1d(
      bufSize, coefficient, reinterpret_cast<fftwf_complex *>(fir), FFTW_ESTIMATE);

    for (size_t idx = 0; idx < forwardPlan.size(); ++idx) {
      forwardPlan[idx] = fftwf_plan_dft_r2c_1d(
        bufSize, buf[idx], reinterpret_cast<fftwf_complex *>(spc), FFTW_ESTIMATE);
    }
    inversePlan = fftwf_plan_dft_c2r_1d(
      bufSize, reinterpret_cast<fftwf_complex *>(spc), flt, FFTW_ESTIMATE);

    reset();
  }

  ~FFTConvolutionFilter()
  {
    fftwf_destroy_plan(firPlan);
    for (auto &fp : forwardPlan) fftwf_destroy_plan(fp);
    fftwf_destroy_plan(inversePlan);

    for (auto &bf : buf) fftwf_free(bf);
    fftwf_free(spc);

    fftwf_free(coefficient);
    fftwf_free(fir);

    fftwf_free(flt);
  }

  inline size_t latency()
  {
    // `(half - 1)` is the latency of FFT buffering.
    // `(half / 2 - 1)` is the latency of FIR filter specific to `refreshFir()`.
    return (half - 1) + (half / 2 - 1);
  }

  void refreshFir(float sampleRate, float cutoffHz, bool isHighpass)
  {
    const auto nyquist = sampleRate / float(2);
    if (cutoffHz > nyquist) cutoffHz = nyquist;

    bool isEven = (half / 2 & 1) == 0;
    size_t end = half;
    if (isEven) --end; // Always use odd length FIR.

    auto mid = float(end - 1) / float(2);
    auto cutoff = float(twopi) * cutoffHz / sampleRate;
    for (size_t idx = 0; idx < end; ++idx) {
      float m = float(idx) - mid;
      float x = cutoff * m;
      coefficient[idx] = x == 0 ? float(1) : std::sin(x) / (x);
    }

    // Apply Nuttall window.
    float tpN = float(twopi) / float(end - 1);
    for (size_t n = 0; n < end; ++n) {
      auto c0 = float(0.3635819);
      auto c1 = float(0.4891775) * std::cos(tpN * n);
      auto c2 = float(0.1365995) * std::cos(tpN * n * float(2));
      auto c3 = float(0.0106411) * std::cos(tpN * n * float(3));
      coefficient[n] *= c0 - c1 + c2 - c3;
    }

    // Normalize to fix FIR scaling.
    float sum = std::accumulate(coefficient, coefficient + half, float(0));
    for (size_t idx = 0; idx < end; ++idx) coefficient[idx] /= sum;

    if (isHighpass) {
      for (size_t idx = 0; idx < end; ++idx) coefficient[idx] = -coefficient[idx];
      coefficient[size_t(mid)] += float(1);
    }

    // Normalize for FFT scaling.
    for (size_t idx = 0; idx < end; ++idx) coefficient[idx] /= float(bufSize);

    fftwf_execute(firPlan);
    reset();

    debugDumpAll();
  }

  void reset()
  {
    front = 0;
    wptr[0] = half;
    wptr[1] = 0;
    rptr = half;

    for (auto &bf : buf) std::fill(bf, bf + bufSize, float(0));
    std::fill(spc, spc + spcSize, std::complex<float>(0, 0));
    std::fill(flt, flt + bufSize, float(0));
  }

  float process(float input)
  {
    buf[0][wptr[0]] = input;
    buf[1][wptr[1]] = input;

    for (auto &w : wptr) {
      if (++w >= bufSize) w = 0;
    }

    if (wptr[front] == 0) {
      fftwf_execute(forwardPlan[front]);
      for (size_t i = 0; i < spcSize; ++i) spc[i] *= fir[i];
      fftwf_execute(inversePlan);

      front ^= 1;
    }

    if (++rptr >= bufSize) rptr = half;
    return flt[rptr];
  }
};

template<typename Sample, size_t nTap> class DirectConvolver {
private:
  std::array<Sample, nTap> co{};
  std::array<Sample, nTap> buf{};

public:
  void setFir(std::vector<float> &source)
  {
    if (source.size() < nTap) {
      std::cerr << "Error in DirectConvolver: source.size() is less than nTap.\n";
      return;
    }

    std::copy(source.begin(), source.begin() + nTap, co.begin());
  }

  void reset() { buf.fill({}); }

  Sample process(Sample input)
  {
    std::rotate(buf.rbegin(), buf.rbegin() + 1, buf.rend());
    buf[0] = input;

    Sample output = 0;
    for (size_t n = 0; n < nTap; ++n) output += buf[n] * co[n];
    return output;
  }
};

class OverlapSaveConvolver {
private:
  static constexpr size_t nBuffer = 2;

  size_t half = 1;
  size_t bufSize = 2;
  size_t spcSize = 1; // spc = spectrum.

  std::array<float *, nBuffer> buf;
  std::complex<float> *spc;
  std::complex<float> *fir;
  float *flt; // filtered.

  std::array<fftwf_plan, nBuffer> forwardPlan;
  fftwf_plan inversePlan;

  size_t front = 0;
  std::array<size_t, nBuffer> wptr{};
  size_t rptr = 0;
  size_t offset = 0;

public:
  void init(size_t nTap, size_t delay = 0)
  {
    offset = delay;

    half = nTap;
    bufSize = 2 * half;
    spcSize = nTap + 1;

    for (size_t idx = 0; idx < nBuffer; ++idx) {
      buf[idx] = (float *)fftwf_malloc(sizeof(float) * bufSize);
    }
    spc = (std::complex<float> *)fftwf_malloc(sizeof(std::complex<float>) * spcSize);
    flt = (float *)fftwf_malloc(sizeof(float) * bufSize);

    fir = (std::complex<float> *)fftwf_malloc(sizeof(std::complex<float>) * spcSize);
    std::fill(fir, fir + spcSize, std::complex<float>(0, 0));

    for (size_t idx = 0; idx < nBuffer; ++idx) {
      forwardPlan[idx] = fftwf_plan_dft_r2c_1d(
        int(bufSize), buf[idx], reinterpret_cast<fftwf_complex *>(spc), FFTW_ESTIMATE);
    }
    inversePlan = fftwf_plan_dft_c2r_1d(
      int(bufSize), reinterpret_cast<fftwf_complex *>(spc), flt, FFTW_ESTIMATE);
  }

  ~OverlapSaveConvolver()
  {
    for (auto &fp : forwardPlan) fftwf_destroy_plan(fp);
    fftwf_destroy_plan(inversePlan);

    for (auto &bf : buf) fftwf_free(bf);
    fftwf_free(spc);
    fftwf_free(fir);
    fftwf_free(flt);
  }

  void setFir(std::vector<float> &source, size_t start, size_t end)
  {
    float *coefficient = (float *)fftwf_malloc(sizeof(float) * bufSize);
    std::copy(source.begin() + start, source.begin() + end, coefficient);
    std::fill(coefficient + half, coefficient + bufSize, float(0));

    // FFT scaling.
    for (size_t idx = 0; idx < half; ++idx) coefficient[idx] /= float(bufSize);

    auto firPlan = fftwf_plan_dft_r2c_1d(
      int(bufSize), coefficient, reinterpret_cast<fftwf_complex *>(fir), FFTW_ESTIMATE);
    fftwf_execute(firPlan);

    fftwf_destroy_plan(firPlan);
    fftwf_free(coefficient);
  }

  void reset()
  {
    wptr[0] = half + offset;
    wptr[1] = offset;
    for (auto &w : wptr) w %= bufSize;
    front = wptr[1] < wptr[0] ? 0 : 1;
    rptr = half + offset % half;

    for (size_t idx = 0; idx < nBuffer; ++idx) {
      std::fill(buf[idx], buf[idx] + bufSize, float(0));
    }
    std::fill(spc, spc + spcSize, std::complex<float>(0, 0));
    std::fill(flt, flt + bufSize, float(0));
  }

  float process(float input)
  {
    buf[0][wptr[0]] = input;
    buf[1][wptr[1]] = input;

    for (auto &w : wptr) {
      if (++w >= bufSize) w = 0;
    }

    if (wptr[front] == 0) {
      fftwf_execute(forwardPlan[front]);
      for (size_t i = 0; i < spcSize; ++i) spc[i] *= fir[i];
      fftwf_execute(inversePlan);

      front ^= 1;
    }

    if (++rptr >= bufSize) rptr = half;
    return flt[rptr];
  }
};

class OverlapAddConvolver {
private:
  static constexpr size_t nBuffer = 2;

  size_t half = 1;
  size_t bufSize = 2;
  size_t spcSize = 1; // spc = spectrum.

  float *buf;
  std::complex<float> *spc;
  std::complex<float> *fir;
  std::array<float *, nBuffer> flt; // filtered.

  fftwf_plan forwardPlan;
  std::array<fftwf_plan, nBuffer> inversePlan;

  size_t front = 0;
  size_t wptr = 0;
  std::array<size_t, nBuffer> rptr{};

public:
  void init(size_t nTap)
  {
    half = nTap;
    bufSize = 2 * half;
    spcSize = nTap + 1;

    buf = (float *)fftwf_malloc(sizeof(float) * bufSize);
    spc = (std::complex<float> *)fftwf_malloc(sizeof(std::complex<float>) * spcSize);
    for (size_t idx = 0; idx < nBuffer; ++idx) {
      flt[idx] = (float *)fftwf_malloc(sizeof(float) * bufSize);
    }

    fir = (std::complex<float> *)fftwf_malloc(sizeof(std::complex<float>) * spcSize);
    std::fill(fir, fir + spcSize, std::complex<float>(0, 0));

    forwardPlan = fftwf_plan_dft_r2c_1d(
      int(bufSize), buf, reinterpret_cast<fftwf_complex *>(spc), FFTW_ESTIMATE);
    for (size_t idx = 0; idx < nBuffer; ++idx) {
      inversePlan[idx] = fftwf_plan_dft_c2r_1d(
        int(bufSize), reinterpret_cast<fftwf_complex *>(spc), flt[idx], FFTW_ESTIMATE);
    }
  }

  ~OverlapAddConvolver()
  {
    fftwf_destroy_plan(forwardPlan);
    for (auto &ip : inversePlan) fftwf_destroy_plan(ip);

    fftwf_free(buf);
    fftwf_free(spc);
    fftwf_free(fir);
    for (auto &fl : flt) fftwf_free(fl);
  }

  void setFir(std::vector<float> &source, size_t start, size_t end)
  {
    float *coefficient = (float *)fftwf_malloc(sizeof(float) * bufSize);
    std::copy(source.begin() + start, source.begin() + end, coefficient);
    std::fill(coefficient + half, coefficient + bufSize, float(0));

    // FFT scaling.
    for (size_t idx = 0; idx < half; ++idx) coefficient[idx] /= float(bufSize);

    auto firPlan = fftwf_plan_dft_r2c_1d(
      int(bufSize), coefficient, reinterpret_cast<fftwf_complex *>(fir), FFTW_ESTIMATE);
    fftwf_execute(firPlan);

    fftwf_destroy_plan(firPlan);
    fftwf_free(coefficient);
  }

  void reset()
  {
    front = 0;
    wptr = 0;
    rptr[0] = 1;
    rptr[1] = 1 + half;

    std::fill(buf, buf + bufSize, float(0));
    std::fill(spc, spc + spcSize, std::complex<float>(0, 0));
    for (size_t idx = 0; idx < nBuffer; ++idx) {
      std::fill(flt[idx], flt[idx] + bufSize, float(0));
    }
  }

  float process(float input)
  {
    buf[wptr] = input;

    if (++wptr >= half) {
      wptr = 0;
      front ^= 1;

      fftwf_execute(forwardPlan);
      for (size_t i = 0; i < spcSize; ++i) spc[i] *= fir[i];
      fftwf_execute(inversePlan[front]);
    }

    float output = flt[0][rptr[0]] + flt[1][rptr[1]];
    for (auto &r : rptr) {
      if (++r >= bufSize) r = 0;
    }
    return output;
  }
};

inline std::vector<float>
getNuttallFir(size_t nTap, float sampleRate, float cutoffHz, bool isHighpass)
{
  const auto nyquist = sampleRate / float(2);
  if (cutoffHz > nyquist) cutoffHz = nyquist;

  bool isEven = (nTap / 2 & 1) == 0;
  size_t end = nTap;
  if (isEven) --end; // Always use odd length FIR.

  std::vector<float> coefficient(nTap);

  auto mid = float(end - 1) / float(2);
  auto cutoff = float(twopi) * cutoffHz / sampleRate;
  for (size_t idx = 0; idx < end; ++idx) {
    float m = float(idx) - mid;
    float x = cutoff * m;
    coefficient[idx] = x == 0 ? float(1) : std::sin(x) / (x);
  }

  // Apply Nuttall window.
  float tpN = float(twopi) / float(end - 1);
  for (size_t n = 0; n < end; ++n) {
    auto c0 = float(0.3635819);
    auto c1 = float(0.4891775) * std::cos(tpN * n);
    auto c2 = float(0.1365995) * std::cos(tpN * n * float(2));
    auto c3 = float(0.0106411) * std::cos(tpN * n * float(3));
    coefficient[n] *= c0 - c1 + c2 - c3;
  }

  // Normalize to fix FIR scaling.
  float sum = std::accumulate(coefficient.begin(), coefficient.end(), float(0));
  for (size_t idx = 0; idx < end; ++idx) coefficient[idx] /= sum;

  if (isHighpass) {
    for (size_t idx = 0; idx < end; ++idx) coefficient[idx] = -coefficient[idx];
    coefficient[size_t(mid)] += float(1);
  }

  return coefficient;
}

/**
FFT convolver without latency.

It must be `lengthInPow2 > minBlockSizeInPowerOfTwo`.

Signal to noise ratio is around -120 dB (1 : 1e-6). Slightly worse than direct
overlap-add.

This is an impelmentation of minimum computation cost solution in following paper:
- William G. Gardner, 1993-11-11, "Efficient Convolution Without Latency"
*/
template<size_t lengthInPow2, size_t minBlockSizeInPow2 = 4> class ImmediateConvolver {
private:
  static constexpr size_t nTap = size_t(1) << lengthInPow2;
  static constexpr size_t nFftConvolver = lengthInPow2 - minBlockSizeInPow2;

  DirectConvolver<float, size_t(1) << minBlockSizeInPow2> firstConvolver;
  // std::array<OverlapAddConvolver, nFftConvolver> fftConvolver;
  std::array<OverlapSaveConvolver, nFftConvolver> fftConvolver;
  std::array<float, nFftConvolver> sumBuffer{};

public:
  ImmediateConvolver()
  {
    static_assert(
      lengthInPow2 > minBlockSizeInPow2,
      "ImmediateConvolver: lengthInPow2 must be greater than minBlockSizeInPow2.");

    for (size_t idx = 0; idx < nFftConvolver; ++idx) {
      fftConvolver[idx].init(size_t(1) << (minBlockSizeInPow2 + idx));
    }

    reset();
  }

  void refreshFir(float sampleRate, float cutoffHz, bool isHighpass)
  {
    auto coefficient = getNuttallFir(nTap, sampleRate, cutoffHz, isHighpass);

    // Set FIR coefficients.
    firstConvolver.setFir(coefficient);
    for (size_t idx = 0; idx < nFftConvolver; ++idx) {
      fftConvolver[idx].setFir(
        coefficient, size_t(1) << (minBlockSizeInPow2 + idx),
        size_t(1) << (minBlockSizeInPow2 + idx + 1));
    }
  }

  void setFir(std::vector<float> &source)
  {
    if (source.size() < nTap) source.resize(nTap);

    firstConvolver.setFir(source);
    for (size_t idx = 0; idx < nFftConvolver; ++idx) {
      size_t start = size_t(1) << (minBlockSizeInPow2 + idx);
      size_t end = size_t(1) << (minBlockSizeInPow2 + idx + 1);
      fftConvolver[idx].setFir(source, start, end);
    }
  }

  void reset()
  {
    firstConvolver.reset();
    for (auto &conv : fftConvolver) conv.reset();
    sumBuffer.fill({});
  }

  float process(float input)
  {
    float output = std::accumulate(sumBuffer.begin(), sumBuffer.end(), float(0));

    for (size_t idx = 0; idx < nFftConvolver; ++idx)
      sumBuffer[idx] = fftConvolver[idx].process(input);

    return output + firstConvolver.process(input);
  }
};

class FixedIntDelayVector {
public:
  std::vector<float> buf{};
  size_t ptr = 0;

  void resize(size_t size) { buf.resize(size); }
  void reset(float value = 0) { std::fill(buf.begin(), buf.end(), value); }

  float process(float input)
  {
    if (++ptr >= buf.size()) ptr -= buf.size();
    auto output = buf[ptr];
    buf[ptr] = input;
    return output;
  }
};

/**
A variation of convolver without latency.

SplitConvolver splits filter kernel into several blocks, then compute the blocks in
different timings. This is a mitigation of CPU load spikes that's caused by FFT.

Memory usage can be reduced by sharing input buffer.
*/
template<size_t nBlock, size_t blockSizeInPow2, size_t minBlockSizeInPow2 = 4>
class SplitConvolver {
public:
  static constexpr size_t nTap = nBlock * (size_t(1) << blockSizeInPow2);
  static constexpr size_t nFftConvolver = nBlock - 2;
  static constexpr size_t blockSize = size_t(1) << blockSizeInPow2;

  ImmediateConvolver<blockSizeInPow2, minBlockSizeInPow2> immediateConvolver;
  OverlapSaveConvolver firstFftConvolver;
  std::array<OverlapSaveConvolver, nFftConvolver> fftConvolver;
  std::array<FixedIntDelayVector, nFftConvolver> outputDelay;
  float buf = 0;

  SplitConvolver()
  {
    firstFftConvolver.init(blockSize, blockSize / nBlock);
    for (size_t idx = 0; idx < nFftConvolver; ++idx) {
      size_t offset = (idx + 2) * blockSize / nBlock;
      fftConvolver[idx].init(blockSize, offset);
      outputDelay[idx].resize((idx + 1) * blockSize + 1);
    }
  }

  void reset()
  {
    immediateConvolver.reset();

    firstFftConvolver.reset();
    for (size_t idx = 0; idx < nFftConvolver; ++idx) {
      fftConvolver[idx].reset();
      outputDelay[idx].reset();
    }
    buf = 0;
  }

  void refreshFir(float sampleRate, float cutoffHz, bool isHighpass)
  {
    auto coefficient = getNuttallFir(nTap, sampleRate, cutoffHz, isHighpass);

    immediateConvolver.setFir(coefficient);
    firstFftConvolver.setFir(coefficient, blockSize, 2 * blockSize);
    for (size_t idx = 0; idx < nFftConvolver; ++idx) {
      fftConvolver[idx].setFir(coefficient, (idx + 2) * blockSize, (idx + 3) * blockSize);
    }
  }

  void setFir(std::vector<float> &coefficient)
  {
    immediateConvolver.setFir(coefficient);
    firstFftConvolver.setFir(coefficient, blockSize, 2 * blockSize);
    for (size_t idx = 0; idx < nFftConvolver; ++idx) {
      fftConvolver[idx].setFir(coefficient, (idx + 2) * blockSize, (idx + 3) * blockSize);
    }
  }

  float process(float input)
  {
    auto output = immediateConvolver.process(input);
    output += buf;
    buf = firstFftConvolver.process(input);
    for (size_t idx = 0; idx < nFftConvolver; ++idx) {
      auto value = fftConvolver[idx].process(input);
      output += outputDelay[idx].process(value);
    }
    return output;
  }
};

template<typename Convolver> void test(std::string name, size_t nLoop = 1)
{
  const float sampleRate = 48000.0f;
  const size_t length = size_t(1.0f * sampleRate);

  std::vector<float> input(length);
  std::vector<float> output(length);

  std::fill(input.begin(), input.end(), float(1));

  Convolver convolver;

  double sumElapsed = 0.0;
  for (size_t n = 0; n < nLoop; ++n) {
    convolver.reset();

    for (size_t i = 0; i < output.size(); ++i) {
      auto start = std::chrono::steady_clock::now();

      output[i] = convolver.process(input[i]);

      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end - start;
      sumElapsed += elapsed.count();
    }
  }
  std::cout << name << ": " << output.size() << "[sample], " << sumElapsed / double(nLoop)
            << "[ms]\n";

  std::string filename = "snd/" + name + ".wav";
  writeWave(filename, output, size_t(sampleRate));
}

template<typename Convolver>
void testFFTConvolutionFilter(std::string name, size_t nLoop = 1)
{
  const float sampleRate = 48000.0f;
  const size_t length = size_t(4.0f * sampleRate);

  std::vector<float> input(length);
  std::vector<float> output(length);
  std::vector<float> time(length);

  std::fill(input.begin(), input.begin() + length / 2, float(1));
  std::fill(input.begin() + length / 2, input.end(), float(0));

  Convolver convolver;
  convolver.refreshFir(sampleRate, 1000.0f, true);

  double sumElapsed = 0.0;
  double timeSpike = 0.0;

  // double phase = 0.0;
  // double delta = 100.0 / sampleRate;
  for (size_t n = 0; n < nLoop; ++n) {
    convolver.reset();

    for (size_t i = 0; i < output.size(); ++i) {
      // float sn = std::sin(float(twopi) * float(phase));
      // phase += delta;
      // phase -= std::floor(phase);

      auto start = std::chrono::steady_clock::now();

      output[i] = convolver.process(input[i]);
      // output[i] = convolver.process(sn);

      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end - start;
      sumElapsed += elapsed.count();
      if (timeSpike < elapsed.count()) timeSpike = elapsed.count();
      time[i] = float(elapsed.count());
    }
  }
  std::cout << name << ": " << output.size() << "[sample], Average"
            << sumElapsed / double(nLoop) << "[ms], Spike " << timeSpike << "[ms]\n";

  std::string filename = "snd/" + name + ".wav";
  writeWave(filename, output, size_t(sampleRate));

  filename = "snd/time_" + name + ".wav";
  writeWave(filename, time, size_t(sampleRate));
}

void testAsymmectricFIR()
{
  // using Convolver = DirectConvolver<float, 32768>;
  // using Convolver = OverlapSaveConvolver;
  using Convolver = SplitConvolver<16, 11>;
  const float sampleRate = 48000.0f;

  Convolver convolver;
  SoundFile snd("data/asym_fir.wav");

  // // OverlapSaveConvolver specific.
  // convolver.init(snd.data.size());
  // convolver.setFir(snd.data, 0, snd.data.size());

  convolver.setFir(snd.data);
  convolver.reset();

  std::vector<float> data((size_t)(2 * sampleRate));
  data[0] = 1;

  for (size_t i = 0; i < data.size(); ++i) data[i] = convolver.process(data[i]);

  writeWave("snd/asym_fir_cpp.wav", data, size_t(sampleRate));
}

int main()
{
  constexpr size_t lengthInPow2 = 15;
  constexpr size_t nLoop = 1;

  std::cout << "dry run\n";
  testFFTConvolutionFilter<ImmediateConvolver<lengthInPow2>>("ImmediateConvolver");

  std::cout << "benchmark\n";
  testFFTConvolutionFilter<FFTConvolutionFilter<size_t(1) << lengthInPow2>>(
    "FFTConvolutionFilter", nLoop);
  testFFTConvolutionFilter<ImmediateConvolver<lengthInPow2>>("ImmediateConvolver", nLoop);

  testFFTConvolutionFilter<SplitConvolver<4, 13>>("SplitConvolver4_13", nLoop);
  testFFTConvolutionFilter<SplitConvolver<8, 12>>("SplitConvolver8_12", nLoop);
  testFFTConvolutionFilter<SplitConvolver<16, 11>>("SplitConvolver16_11", nLoop);
  testFFTConvolutionFilter<SplitConvolver<32, 10>>("SplitConvolver32_10", nLoop);
  testFFTConvolutionFilter<SplitConvolver<64, 9>>("SplitConvolver64_9", nLoop);
  testFFTConvolutionFilter<SplitConvolver<128, 8>>("SplitConvolver128_8", nLoop);

  testAsymmectricFIR();

  return 0;
}
