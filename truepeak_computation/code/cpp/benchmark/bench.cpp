#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include <random>
#include <sndfile.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

constexpr double pi = 3.14159265358979323846264338;
constexpr double twopi = 2.0 * pi;

template<typename T> T toDecibel(T x) { return T(20) * std::log10(x); }

void writeWave(const char *filename, std::vector<float> &buffer, const size_t &samplerate)
{
  SF_INFO sfinfo;
  memset(&sfinfo, 0, sizeof(sfinfo));
  sfinfo.samplerate = int(samplerate);
  sfinfo.frames = buffer.size();
  sfinfo.channels = 1;
  sfinfo.format = (SF_FORMAT_WAV | SF_FORMAT_FLOAT);

  SNDFILE *file = sf_open(filename, SFM_WRITE, &sfinfo);
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

  SoundFile(const char *path)
  {
    SF_INFO sfinfo;
    memset(&sfinfo, 0, sizeof(sfinfo));

    SNDFILE *file = sf_open(path, SFM_READ, &sfinfo);
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

// 12.2 sample delay.
template<typename Sample> struct ThiranAP0 {
  constexpr static size_t nSection = 6;

  constexpr static std::array<std::array<Sample, 5>, nSection> co{{
    {Sample(-4.727394451998352e-07), Sample(4.858748704357417e-07),
     Sample(-9.312803872380755e-06), Sample(-0.0521727803026907),
     Sample(0.050762310865564084)},
    {Sample(1.0), Sample(-4.442270342973251), Sample(18.48262423641914),
     Sample(-0.24034846384097586), Sample(0.05410487099713843)},
    {Sample(1.0), Sample(-7.436178165899102), Sample(18.277812772383758),
     Sample(-0.2338673463646018), Sample(0.0)},
    {Sample(1.0), Sample(1.7134397024528651), Sample(12.808242419750542),
     Sample(-0.4068417954873905), Sample(0.05471114145073639)},
    {Sample(1.0), Sample(3.4981346117372594), Sample(5.803790586563618),
     Sample(0.13377633295031668), Sample(0.07807472463653542)},
    {Sample(1.0), Sample(-4.275928279619698), Sample(0.0), Sample(0.6027327415699312),
     Sample(0.17230118576557518)},
  }};
};

// 12.4 sample delay.
template<typename Sample> struct ThiranAP1 {
  constexpr static size_t nSection = 6;

  constexpr static std::array<std::array<Sample, 5>, nSection> co{{
    {Sample(-1.3805232842433398e-06), Sample(1.5219474676742606e-06),
     Sample(-2.3938570964448888e-05), Sample(-0.06357720642282934),
     Sample(0.05766941085554396)},
    {Sample(1.0), Sample(-4.02886282984051), Sample(14.72542178017142),
     Sample(-0.2735991464275419), Sample(0.06790976957593897)},
    {Sample(1.0), Sample(-6.654709260511208), Sample(14.527156124516953),
     Sample(-0.45808754332035756), Sample(0.06883659757138261)},
    {Sample(1.0), Sample(1.002985687612675), Sample(10.780414439513091),
     Sample(-0.2624929943256504), Sample(0.0)},
    {Sample(1.0), Sample(2.7465004855109503), Sample(4.754837636026572),
     Sample(0.09303776707660172), Sample(0.09276081226847141)},
    {Sample(1.0), Sample(-3.809625481887683), Sample(0.0), Sample(0.5776223492262271),
     Sample(0.21031212347256095)},
  }};
};

// 12.6 sample delay.
template<typename Sample> struct ThiranAP2 {
  constexpr static size_t nSection = 6;

  constexpr static std::array<std::array<Sample, 5>, nSection> co{{
    {Sample(-2.9167632147054213e-06), Sample(3.237806542173795e-06),
     Sample(-4.6867759764576194e-05), Sample(-0.06908387681505378),
     Sample(0.0622338944587217)},
    {Sample(1.0), Sample(-3.7602717797285585), Sample(12.388577634342536),
     Sample(-0.3035273209496412), Sample(0.08071951676098124)},
    {Sample(1.0), Sample(-6.128153542196332), Sample(12.232549726380213),
     Sample(-0.5009710713851159), Sample(0.08174910565402835)},
    {Sample(1.0), Sample(0.5026375541049284), Sample(9.534195796338215),
     Sample(-0.28618703542249796), Sample(0.0)},
    {Sample(1.0), Sample(2.2636594104437835), Sample(4.226231256487883),
     Sample(0.05271944953112453), Sample(0.1048856160877333)},
    {Sample(1.0), Sample(-3.4942183824773094), Sample(0.0), Sample(0.535621283612614),
     Sample(0.23661743508826213)},
  }};
};

// 12.8 sample delay.
template<typename Sample> struct ThiranAP3 {
  constexpr static size_t nSection = 6;

  constexpr static std::array<std::array<Sample, 5>, nSection> co{{
    {Sample(-5.319560569852848e-06), Sample(5.79452129833411e-06),
     Sample(-8.140073150945665e-05), Sample(-0.07118512562337313),
     Sample(0.06535027967451192)},
    {Sample(1.0), Sample(-3.5606780149659247), Sample(10.711064576149075),
     Sample(-0.33242988963901143), Sample(0.09336140146393108)},
    {Sample(1.0), Sample(-5.720747652112522), Sample(10.589722958963094),
     Sample(-0.5402169324241914), Sample(0.09443117670550202)},
    {Sample(1.0), Sample(0.11299402378611216), Sample(8.582256726889149),
     Sample(-0.30770693924232256), Sample(0.0)},
    {Sample(1.0), Sample(1.8964517011722102), Sample(3.8832043429973853),
     Sample(0.013166003695985196), Sample(0.11651946939164404)},
    {Sample(1.0), Sample(-3.249845461601337), Sample(0.0), Sample(0.488372883232914),
     Sample(0.2575192834761079)},
  }};
};

template<typename Sample> struct BS1770FIR {
  constexpr static size_t bufferSize = 12;
  constexpr static size_t intDelay = 5;

  constexpr static std::array<std::array<Sample, 12>, 4> coefficient{{
    {Sample(0.001708984375), Sample(0.010986328125), Sample(-0.0196533203125),
     Sample(0.033203125), Sample(-0.0594482421875), Sample(0.1373291015625),
     Sample(0.97216796875), Sample(-0.102294921875), Sample(0.047607421875),
     Sample(-0.026611328125), Sample(0.014892578125), Sample(-0.00830078125)},
    {Sample(-0.0291748046875), Sample(0.029296875), Sample(-0.0517578125),
     Sample(0.089111328125), Sample(-0.16650390625), Sample(0.465087890625),
     Sample(0.77978515625), Sample(-0.2003173828125), Sample(0.1015625),
     Sample(-0.0582275390625), Sample(0.0330810546875), Sample(-0.0189208984375)},
    {Sample(-0.0189208984375), Sample(0.0330810546875), Sample(-0.0582275390625),
     Sample(0.1015625), Sample(-0.2003173828125), Sample(0.77978515625),
     Sample(0.465087890625), Sample(-0.16650390625), Sample(0.089111328125),
     Sample(-0.0517578125), Sample(0.029296875), Sample(-0.0291748046875)},
    {Sample(-0.00830078125), Sample(0.014892578125), Sample(-0.026611328125),
     Sample(0.047607421875), Sample(-0.102294921875), Sample(0.97216796875),
     Sample(0.1373291015625), Sample(-0.0594482421875), Sample(0.033203125),
     Sample(-0.0196533203125), Sample(0.010986328125), Sample(0.001708984375)},
  }};
};

template<typename Sample> struct SOCPFIR {
  constexpr static size_t bufferSize = 12;
  constexpr static size_t intDelay = 5;

  constexpr static std::array<std::array<Sample, 12>, 4> coefficient{{
    {Sample(0.0003685641156728883), Sample(-0.002650918350234966),
     Sample(0.010431691474785761), Sample(-0.02960494116533279),
     Sample(0.06956040589153337), Sample(-0.16056336790983874),
     Sample(0.9619082375782309), Sample(0.19041213413772073),
     Sample(-0.052224367121096446), Sample(0.015475683352663345),
     Sample(-0.003574936706560431), Sample(0.00046731224485537984)},
    {Sample(0.0006010168055805836), Sample(-0.0042956015834363865),
     Sample(0.01674961587627823), Sample(-0.04685180005920451),
     Sample(0.10717670363659904), Sample(-0.23123864297504026),
     Sample(0.8077824959976656), Sample(0.4262349123407106), Sample(-0.09859571003641601),
     Sample(0.0279557646690599), Sample(-0.0063265909046811885),
     Sample(0.0008170930533285713)},
    {Sample(0.000612338707572119), Sample(-0.004350284800271867),
     Sample(0.016818981495704766), Sample(-0.04643237217237933),
     Sample(0.10377074972252216), Sample(-0.2121399428304628), Sample(0.564378803898816),
     Sample(0.669754196116537), Sample(-0.11798394837792067), Sample(0.03169486170413676),
     Sample(-0.007006550941503337), Sample(0.0008928958479879955)},
    {Sample(0.0003898604923927281), Sample(-0.0027537715341607303),
     Sample(0.010562059165944073), Sample(-0.02881152759651303),
     Sample(0.06308327336157231), Sample(-0.12339323584336521),
     Sample(0.2768518164085285), Sample(0.875685999063749), Sample(-0.08993915437525968),
     Sample(0.022579480010111178), Sample(-0.00485841074687993),
     Sample(0.0006099943345506027)},
  }};
};

template<typename Sample, uint8_t order> class FractionalDelayLagrange {
private:
  std::array<Sample, order> xd{};
  std::array<Sample, order> diff{};

public:
  void reset()
  {
    xd.fill(0.0);
    diff.fill(0.0);
  }

  void push(Sample input)
  {
    diff[0] = input - xd[0];
    for (size_t i = 1; i < order; ++i) diff[i] = diff[i - 1] - xd[i];

    xd[0] = input;
    for (size_t i = 1; i < order; ++i) xd[i] = diff[i - 1];
  }

  Sample at(Sample fraction)
  {
    Sample delta = fraction + (order - 1) / 2;
    Sample sig = 0.0;
    for (uint8_t i = order; i > 0; --i) sig = ((i - 1) - delta) / i * (diff[i - 1] + sig);
    return sig + xd[0];
  }

  Sample process(Sample input, Sample fraction)
  {
    push(input);
    return at(fraction);
  }
};

// SOS: Second order sections.
template<typename Sample, typename IIR> class SOSFilter {
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
      // clang-format off
      y0[i] =
        + IIR::co[i][0] * x0[i]
        + IIR::co[i][1] * x1[i]
        + IIR::co[i][2] * x2[i]
        - IIR::co[i][3] * y1[i]
        - IIR::co[i][4] * y2[i];
      // clang-format on
    }

    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }

  inline Sample output() { return y0[IIR::nSection - 1]; }

  Sample process(Sample input)
  {
    push(input);
    return output();
  }

  std::array<Sample, IIR::nSection> x0{};
  std::array<Sample, IIR::nSection> x1{};
  std::array<Sample, IIR::nSection> x2{};
  std::array<Sample, IIR::nSection> y0{};
  std::array<Sample, IIR::nSection> y1{};
  std::array<Sample, IIR::nSection> y2{};
};

template<typename Sample, typename FractionalDelayFIR> class TruePeakMeterFIR {
  std::array<Sample, FractionalDelayFIR::bufferSize> buf{};

public:
  void reset() { buf.fill(Sample(0)); }

  Sample process(Sample input)
  {
    for (size_t i = 0; i < buf.size() - 1; ++i) buf[i] = buf[i + 1];
    buf.back() = input;

    Sample max = std::fabs(buf[FractionalDelayFIR::intDelay]);
    for (const auto &phase : FractionalDelayFIR::coefficient) {
      Sample sum = 0;
      for (size_t i = 0; i < phase.size(); ++i) sum += buf[i] * phase[i];
      max = std::max(max, std::fabs(sum));
    }
    return max;
  }
};

template<typename Sample, uint8_t order> class TruePeakMeterLagrange {
private:
  std::array<Sample, (order - 1) / 2> buf{};
  uint8_t bufptr = 0;

public:
  FractionalDelayLagrange<Sample, order> delay;

  void reset()
  {
    buf.fill(0);
    bufptr = 0;
    delay.reset();
  }

  Sample process(Sample input)
  {
    delay.push(input);

    auto out0 = delay.at(Sample(1.0 / 5.0));
    auto out1 = delay.at(Sample(2.0 / 5.0));
    auto out2 = delay.at(Sample(3.0 / 5.0));
    auto out3 = delay.at(Sample(4.0 / 5.0));

    buf[bufptr] = input;
    bufptr = (bufptr + 1) % uint8_t(buf.size());

    auto output = std::max(std::fabs(out0), std::fabs(buf[bufptr]));
    output = std::max(std::fabs(out1), output);
    output = std::max(std::fabs(out2), output);
    return std::max(std::fabs(out3), output);
  }
};

template<typename Sample> class TruePeakMeterThiranAP {
private:
  std::array<Sample, 12> buf{};
  uint8_t bufptr = 0;

public:
  SOSFilter<Sample, ThiranAP0<Sample>> thiran0;
  SOSFilter<Sample, ThiranAP1<Sample>> thiran1;
  SOSFilter<Sample, ThiranAP2<Sample>> thiran2;
  SOSFilter<Sample, ThiranAP3<Sample>> thiran3;

  void reset()
  {
    buf.fill(0);
    bufptr = 0;

    thiran0.reset();
    thiran1.reset();
    thiran2.reset();
    thiran3.reset();
  }

  Sample process(Sample input)
  {
    auto out0 = thiran0.process(input);
    auto out1 = thiran1.process(input);
    auto out2 = thiran2.process(input);
    auto out3 = thiran3.process(input);

    buf[bufptr] = input;
    bufptr = (bufptr + 1) % uint8_t(buf.size());

    auto output = std::max(std::fabs(out0), std::fabs(buf[bufptr]));
    output = std::max(std::fabs(out1), output);
    output = std::max(std::fabs(out2), output);
    return std::max(std::fabs(out3), output);
  }
};

void testMeter()
{
  FractionalDelayLagrange<float, 15> lagrange;
  std::cout << "Lagurange: " << lagrange.process(1.0f, 0.0f) << "\n";

  TruePeakMeterFIR<float, BS1770FIR<float>> bs1770;
  std::cout << "bs1770: " << bs1770.process(1.0f) << "\n";

  TruePeakMeterFIR<float, SOCPFIR<float>> socp;
  std::cout << "socp: " << socp.process(1.0f) << "\n";

  SOSFilter<float, ThiranAP0<float>> thiran0;
  SOSFilter<float, ThiranAP1<float>> thiran1;
  SOSFilter<float, ThiranAP2<float>> thiran2;
  SOSFilter<float, ThiranAP3<float>> thiran3;
  std::cout << "thiran0: " << thiran0.process(1.0f) << "\n";
  std::cout << "thiran1: " << thiran1.process(1.0f) << "\n";
  std::cout << "thiran2: " << thiran2.process(1.0f) << "\n";
  std::cout << "thiran3: " << thiran3.process(1.0f) << "\n";
}

struct BenchmarkResult {
  constexpr static size_t nLoop = 1;
  std::string path;
  size_t samples;
  double elapsedMilliSeconds;
  float truepeak;
  float dbtp;
  bool inTolerance;
};

template<typename Meter>
nlohmann::json benchmark(SoundFile &snd, std::vector<double> tolerance)
{
  Meter meter;
  std::vector<float> wav(snd.data.size());

  double sumElapsed = 0.0;
  for (size_t n = 0; n < BenchmarkResult::nLoop; ++n) {
    meter.reset();
    auto start = std::chrono::steady_clock::now();
    for (size_t i = 0; i < snd.data.size(); ++i) {
      wav[i] = meter.process(snd.data[i]);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    sumElapsed += elapsed.count();
  }
  auto truepeak = *std::max_element(wav.begin(), wav.end());
  auto dbtp = toDecibel(truepeak);

  double upperBoundDb = tolerance[0] + tolerance[1];
  double lowerBoundDb = tolerance[0] + tolerance[2];
  bool inTolerance = (dbtp <= upperBoundDb) && (dbtp >= lowerBoundDb);

  std::cout << BenchmarkResult::nLoop << " * " << wav.size() << "[sample], " << sumElapsed
            << "[ms], " << dbtp << "[dB TP], " << inTolerance << "\n";

  return nlohmann::json{
    {"ElapsedMilliSeconds", sumElapsed},
    {"True-peak", truepeak},
    {"dBTP", dbtp},
    {"inTolerance", inTolerance},
  };
}

void testEBUTech3341()
{
  std::vector<std::vector<double>> tolerance{{
    {-6.0, +0.2, -0.4}, // 15
    {-6.0, +0.2, -0.4}, // 16
    {-6.0, +0.2, -0.4}, // 17
    {-6.0, +0.2, -0.4}, // 18
    {+3.0, +0.2, -0.4}, // 19
    {0.0, +0.2, -0.4},  // 20
    {0.0, +0.2, -0.4},  // 21
    {0.0, +0.2, -0.4},  // 22
    {0.0, +0.2, -0.4},  // 23
    {0.0, 0.0, 0.0},    // worstsinc
  }};

  std::vector<std::string> path{
    "data/ebu-loudness-test-setv05/seq-3341-15-24bit.wav.wav",
    "data/ebu-loudness-test-setv05/seq-3341-16-24bit.wav.wav",
    "data/ebu-loudness-test-setv05/seq-3341-17-24bit.wav.wav",
    "data/ebu-loudness-test-setv05/seq-3341-18-24bit.wav.wav",
    "data/ebu-loudness-test-setv05/seq-3341-19-24bit.wav.wav",
    "data/ebu-loudness-test-setv05/seq-3341-20-24bit.wav.wav",
    "data/ebu-loudness-test-setv05/seq-3341-21-24bit.wav.wav",
    "data/ebu-loudness-test-setv05/seq-3341-22-24bit.wav.wav",
    "data/ebu-loudness-test-setv05/seq-3341-23-24bit.wav.wav",
    "data/worstsinc/worst_48000Hz_01sec.wav",
  };

  std::string dataDir("../../../");

  nlohmann::json result = nlohmann::json::array({});
  for (size_t idx = 0; idx < path.size(); ++idx) {
    auto dataPath = dataDir + path[idx];
    SoundFile snd(dataPath.c_str());

    std::cout << "\n--- " << path[idx] << "\n";

    result.insert(
      result.end(),
      nlohmann::json{
        {"Path", path[idx]},
        {"Sample", snd.data.size()},
        {"Loop", BenchmarkResult::nLoop},
        {"ITU-R BS.1770 FIR",
         benchmark<TruePeakMeterFIR<float, BS1770FIR<float>>>(snd, tolerance[idx])},
        {"SOCP FIR",
         benchmark<TruePeakMeterFIR<float, SOCPFIR<float>>>(snd, tolerance[idx])},
        {"Thiran Allpass", benchmark<TruePeakMeterThiranAP<float>>(snd, tolerance[idx])},
        {"Lagurange Interpolation",
         benchmark<TruePeakMeterLagrange<float, 11>>(snd, tolerance[idx])},
      });
  }
  std::ofstream ofs;
  ofs.open("benchmark.json");
  ofs << result;
  ofs.close();
}

void measure()
{
  namespace fs = std::filesystem;

  nlohmann::json result;
  std::string dirPath("../../../");
  std::string dataDir("data");
  for (const auto &entry : fs::recursive_directory_iterator(dirPath + dataDir)) {
    auto path = entry.path();
    if (fs::is_directory(path)) continue;

    std::string unixPath;
    if (path.begin() != path.end()) {
      unixPath += path.begin()->string();
      for (auto it = ++(path.begin()); it != path.end(); ++it) {
        unixPath.append("/").append(it->string());
      }
    }

    SoundFile snd(unixPath.c_str());
    if (snd.data.size() == 0) {
      std::cout << "Skipped: " << unixPath << "\n";
      continue;
    }

    std::cout << "\n--- " << unixPath << "\n";
    std::vector<double> emptyTolerance{0.0, 0.0, 0.0};
    auto label = unixPath.erase(0, dirPath.size());
    result[label] = nlohmann::json{
      {"Path", unixPath},
      {"Sample", snd.data.size()},
      {"Loop", BenchmarkResult::nLoop},
      {"ITU-R BS.1770 FIR",
       benchmark<TruePeakMeterFIR<float, BS1770FIR<float>>>(snd, emptyTolerance)},
      {"SOCP FIR",
       benchmark<TruePeakMeterFIR<float, SOCPFIR<float>>>(snd, emptyTolerance)},
      {"Thiran Allpass", benchmark<TruePeakMeterThiranAP<float>>(snd, emptyTolerance)},
      {"Lagurange Interpolation",
       benchmark<TruePeakMeterLagrange<float, 11>>(snd, emptyTolerance)},
    };
  }

  auto now = std::chrono::system_clock::now();
  auto epoch = std::to_string(
    std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count());

  std::ofstream ofs;
  ofs.open(std::string("measure_cpp") + epoch + ".json");
  ofs << result;
  ofs.close();
}

int main()
{
  // testMeter();
  // testEBUTech3341();
  measure();
  return 0;
}
