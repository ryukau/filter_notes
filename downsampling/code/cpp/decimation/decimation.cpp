#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <sndfile.h>
#include <sstream>
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

/**
```python
import numpy
from scipy import signal
sos = signal.butter(16, 20500*2/44100/8, "low", output="sos")
```
*/
template<typename Sample> struct SosButterDecimation8 {
  constexpr static size_t oversample = 8;
  constexpr static size_t nSection = 8;

  constexpr static std::array<std::array<Sample, 5>, nSection> co{{
    {Sample(2.78583757e-13), Sample(5.57167513e-13), Sample(2.78583757e-13),
     Sample(-1.37840700e+00), Sample(4.75668194e-01)},
    {Sample(1.00000000e+00), Sample(2.00000000e+00), Sample(1.00000000e+00),
     Sample(-1.39243564e+00), Sample(4.90686697e-01)},
    {Sample(1.00000000e+00), Sample(2.00000000e+00), Sample(1.00000000e+00),
     Sample(-1.42080012e+00), Sample(5.21052591e-01)},
    {Sample(1.00000000e+00), Sample(2.00000000e+00), Sample(1.00000000e+00),
     Sample(-1.46409826e+00), Sample(5.67405872e-01)},
    {Sample(1.00000000e+00), Sample(2.00000000e+00), Sample(1.00000000e+00),
     Sample(-1.52317720e+00), Sample(6.30653459e-01)},
    {Sample(1.00000000e+00), Sample(2.00000000e+00), Sample(1.00000000e+00),
     Sample(-1.59904996e+00), Sample(7.11879848e-01)},
    {Sample(1.00000000e+00), Sample(2.00000000e+00), Sample(1.00000000e+00),
     Sample(-1.69274029e+00), Sample(8.12181022e-01)},
    {Sample(1.00000000e+00), Sample(2.00000000e+00), Sample(1.00000000e+00),
     Sample(-1.80501230e+00), Sample(9.32375009e-01)},
  }};
};

/**
```python
import numpy
from scipy import signal
sos = signal.ellip(16, 0.01, 120, 22000, "low", output="sos", fs=48000 * 16)
```
*/
template<typename Sample> struct SosEllipticDecimation8 {
  constexpr static size_t oversample = 8;
  constexpr static size_t nSection = 8;

  constexpr static std::array<std::array<Sample, 5>, nSection> co{{
    {Sample(2.9603570511348458e-06), Sample(1.9134308325212554e-06),
     Sample(2.960357051134845e-06), Sample(-1.7658654066203068),
     Sample(0.7820540053868873)},
    {Sample(1.0), Sample(-1.2291620715812779), Sample(1.0), Sample(-1.7791504726807224),
     Sample(0.811591673961548)},
    {Sample(1.0), Sample(-1.629608485127653), Sample(1.0), Sample(-1.799429713423469),
     Sample(0.8565233962004007)},
    {Sample(1.0), Sample(-1.755433170484052), Sample(1.0), Sample(-1.819719169754641),
     Sample(0.9010752193191813)},
    {Sample(1.0), Sample(-1.8071956550378283), Sample(1.0), Sample(-1.8363783422750304),
     Sample(0.936837780001277)},
    {Sample(1.0), Sample(-1.831500447999239), Sample(1.0), Sample(-1.8490580542244668),
     Sample(0.9625736181992885)},
    {Sample(1.0), Sample(-1.8432523070089726), Sample(1.0), Sample(-1.8591271322320881),
     Sample(0.9805644550845712)},
    {Sample(1.0), Sample(-1.84810578527457), Sample(1.0), Sample(-1.8686026300558227),
     Sample(0.9940001024564757)},
  }};
};

/**
```python
import numpy
from scipy import signal
sos = signal.ellip(16, 0.01, 120, 20500, "low", output="sos", fs=48000 * 16)
```
*/
template<typename Sample> struct SosEllipticDecimation16 {
  constexpr static size_t oversample = 16;
  constexpr static size_t nSection = 8;

  constexpr static std::array<std::array<Sample, 5>, nSection> co{{
    {Sample(1.325527960537483e-06), Sample(-1.0486296880714946e-06),
     Sample(1.3255279605374831e-06), Sample(-1.8869513625870566),
     Sample(0.8907702524266276)},
    {Sample(1.0), Sample(-1.7990799255799657), Sample(0.9999999999999999),
     Sample(-1.8984243919196817), Sample(0.90603953110456)},
    {Sample(1.0), Sample(-1.911565142173381), Sample(1.0), Sample(-1.9156790452684018),
     Sample(0.9289810487694654)},
    {Sample(1.0), Sample(-1.943108432116689), Sample(1.0), Sample(-1.9325800033745015),
     Sample(0.9513949185576018)},
    {Sample(1.0), Sample(-1.9556189334610117), Sample(1.0000000000000002),
     Sample(-1.9460569896516668), Sample(0.9691513322359767)},
    {Sample(1.0), Sample(-1.9614032628799634), Sample(1.0), Sample(-1.955820788205649),
     Sample(0.9818014494533305)},
    {Sample(1.0), Sample(-1.9641798599952223), Sample(1.0), Sample(-1.9628792094377048),
     Sample(0.9905808098567532)},
    {Sample(1.0), Sample(-1.965322766624071), Sample(0.9999999999999999),
     Sample(-1.968577366856729), Sample(0.9971005403967146)},
  }};
};

// SOS: Second order sections.
template<typename Sample, typename IIR> class SosFilter {
public:
  void reset()
  {
    x1.fill(0);
    x2.fill(0);
    y1.fill(0);
    y2.fill(0);
  }

  void push(Sample input)
  {
    for (size_t i = 0; i < IIR::nSection; ++i) {
      // clang-format off
      auto y0 =
        + IIR::co[i][0] * input
        + IIR::co[i][1] * x1[i]
        + IIR::co[i][2] * x2[i]
        - IIR::co[i][3] * y1[i]
        - IIR::co[i][4] * y2[i];
      // clang-format on

      x2[i] = x1[i];
      x1[i] = input;
      y2[i] = y1[i];
      y1[i] = y0;
      input = y0;
    }
  }

  inline Sample output() { return y1[IIR::nSection - 1]; }

  Sample process(const std::array<Sample, IIR::oversample> &input)
  {
    for (const auto &value : input) push(value);
    return output();
  }

  std::array<Sample, IIR::nSection> x1{};
  std::array<Sample, IIR::nSection> x2{};
  std::array<Sample, IIR::nSection> y1{};
  std::array<Sample, IIR::nSection> y2{};
};

// SOS: Second order sections.
// This implementation adds 1 sample delay for each section, but sounds almost as same as
// correct implementation (SosFilter). The length of latency caused by added delay is
// `nSection` samples. Slightly faster than no latency implementation.
template<typename Sample, typename IIR> class SosFilterLatency {
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

  Sample process(const std::array<Sample, IIR::oversample> &input)
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

// Direct form II.
template<typename Sample, typename IIR> class SosFilterDF2 {
public:
  void reset()
  {
    sig = 0;
    v0.fill(0);
    v1.fill(0);
    v2.fill(0);
  }

  void push(Sample input)
  {
    sig = input;
    for (size_t i = 0; i < IIR::nSection; ++i) {
      v0[i] = sig - IIR::co[i][3] * v1[i] - IIR::co[i][4] * v2[i];
      sig = IIR::co[i][0] * v0[i] + IIR::co[i][1] * v1[i] + IIR::co[i][2] * v2[i];
    }

    v2 = v1;
    v1 = v0;
  }

  inline Sample output() { return sig; }

  Sample process(const std::array<Sample, IIR::oversample> &input)
  {
    for (const auto &value : input) push(value);
    return output();
  }

  Sample sig = 0;
  std::array<Sample, IIR::nSection> v0{};
  std::array<Sample, IIR::nSection> v1{};
  std::array<Sample, IIR::nSection> v2{};
};

// Direct form II.
// This implementation adds 1 sample delay for each section, but sounds almost as same as
// correct implementation (SosFilterDF2). The length of latency caused by added delay is
// `nSection` samples. Slightly faster than no latency implementation.
template<typename Sample, typename IIR> class SosFilterDF2Latency {
public:
  void reset()
  {
    x0.fill(0);
    y0.fill(0);
    v0.fill(0);
    v1.fill(0);
    v2.fill(0);
  }

  void push(Sample input)
  {
    x0[0] = input;
    for (size_t i = 1; i < IIR::nSection; ++i) x0[i] = y0[i - 1];

    for (size_t i = 0; i < IIR::nSection; ++i) {
      v0[i] = x0[i] - IIR::co[i][3] * v1[i] - IIR::co[i][4] * v2[i];
      y0[i] = IIR::co[i][0] * v0[i] + IIR::co[i][1] * v1[i] + IIR::co[i][2] * v2[i];
    }

    v2 = v1;
    v1 = v0;
  }

  inline Sample output()
  {
    // return sig;
    return y0[IIR::nSection - 1];
  }

  Sample process(const std::array<Sample, IIR::oversample> &input)
  {
    for (const auto &value : input) push(value);
    return output();
  }

  std::array<Sample, IIR::nSection> x0{};
  std::array<Sample, IIR::nSection> y0{};
  std::array<Sample, IIR::nSection> v0{};
  std::array<Sample, IIR::nSection> v1{};
  std::array<Sample, IIR::nSection> v2{};
};

/*
Fir filter specialized to 8x oversampling.

±0.0052 dB ripple at pass band.
-64 dB attenuation at stop band.
For fs=48000Hz, aliasing overwraps to 19000Hz.

```python
import scipy.signal as signal
oversample = 8
signal.remez(128, [0, 0.375, 0.6, oversample], [1, 0], fs=oversample)
```
*/
template<typename Sample> struct FirRemez128Decimation8 {
  constexpr static size_t bufferSize = 16;
  constexpr static size_t nPhase = 8;

  constexpr static std::array<std::array<Sample, bufferSize>, nPhase> coefficient{{
    {Sample(-4.3367201090342294e-05), Sample(0.0004346243146011135),
     Sample(-0.0009251472353283694), Sample(0.0016612171550778898),
     Sample(-0.002677133553760478), Sample(0.00404890211257882),
     Sample(-0.006109727533258657), Sample(0.010979319633459658),
     Sample(0.12108731873320519), Sample(-0.004017909358618499),
     Sample(0.0005399248291031557), Sample(0.00042125979481843253),
     Sample(-0.0006423711752143175), Sample(0.000564226921945931),
     Sample(-0.0003828699690866557), Sample(0.00020312413171973427)},
    {Sample(-0.00044208529422296554), Sample(0.0006493789186035279),
     Sample(-0.001422335334297734), Sample(0.0026754675838967304),
     Sample(-0.004596497856809998), Sample(0.007585024502148825),
     Sample(-0.01297320726650012), Sample(0.02904408396554225),
     Sample(0.11511049947575919), Sample(-0.015095824264959638),
     Sample(0.006126939037219295), Sample(-0.002783579509504792),
     Sample(0.0012050070838781225), Sample(-0.00044648276821080427),
     Sample(0.00011763293979479671), Sample(-6.640131982004079e-06)},
    {Sample(-0.00025755503195241914), Sample(0.0008024597718016631),
     Sample(-0.0017785452728898907), Sample(0.0034230111708849664),
     Sample(-0.006070385040979504), Sample(0.010452060026711708),
     Sample(-0.019023606544203667), Sample(0.0489283284941896), Sample(0.10370363811315),
     Sample(-0.02186207622888707), Sample(0.010060638744526014),
     Sample(-0.005176024970725293), Sample(0.0026292526180395877),
     Sample(-0.0012390661848673606), Sample(0.0005106943331457618),
     Sample(-0.00016727396282846716)},
    {Sample(-0.0003017475020786858), Sample(0.0008500359484970119),
     Sample(-0.0019041568172862991), Sample(0.0037347545183922535),
     Sample(-0.006791848179837525), Sample(0.012090234525293302),
     Sample(-0.023170513523598623), Sample(0.06910081101343374),
     Sample(0.08789928618273786), Sample(-0.024389467061261057),
     Sample(0.012046987395531172), Sample(-0.006529575330676864),
     Sample(0.0034868084142677253), Sample(-0.0017331779258313819),
     Sample(0.0007579962630174413), Sample(-0.00026481024579566364)},
    {Sample(-0.00026481024579566364), Sample(0.0007579962630174413),
     Sample(-0.0017331779258313819), Sample(0.0034868084142677253),
     Sample(-0.006529575330676864), Sample(0.012046987395531172),
     Sample(-0.024389467061261057), Sample(0.08789928618273786),
     Sample(0.06910081101343374), Sample(-0.023170513523598623),
     Sample(0.012090234525293302), Sample(-0.006791848179837525),
     Sample(0.0037347545183922535), Sample(-0.0019041568172862991),
     Sample(0.0008500359484970119), Sample(-0.0003017475020786858)},
    {Sample(-0.00016727396282846716), Sample(0.0005106943331457618),
     Sample(-0.0012390661848673606), Sample(0.0026292526180395877),
     Sample(-0.005176024970725293), Sample(0.010060638744526014),
     Sample(-0.02186207622888707), Sample(0.10370363811315), Sample(0.0489283284941896),
     Sample(-0.019023606544203667), Sample(0.010452060026711708),
     Sample(-0.006070385040979504), Sample(0.0034230111708849664),
     Sample(-0.0017785452728898907), Sample(0.0008024597718016631),
     Sample(-0.00025755503195241914)},
    {Sample(-6.640131982004079e-06), Sample(0.00011763293979479671),
     Sample(-0.00044648276821080427), Sample(0.0012050070838781225),
     Sample(-0.002783579509504792), Sample(0.006126939037219295),
     Sample(-0.015095824264959638), Sample(0.11511049947575919),
     Sample(0.02904408396554225), Sample(-0.01297320726650012),
     Sample(0.007585024502148825), Sample(-0.004596497856809998),
     Sample(0.0026754675838967304), Sample(-0.001422335334297734),
     Sample(0.0006493789186035279), Sample(-0.00044208529422296554)},
    {Sample(0.00020312413171973427), Sample(-0.0003828699690866557),
     Sample(0.000564226921945931), Sample(-0.0006423711752143175),
     Sample(0.00042125979481843253), Sample(0.0005399248291031557),
     Sample(-0.004017909358618499), Sample(0.12108731873320519),
     Sample(0.010979319633459658), Sample(-0.006109727533258657),
     Sample(0.00404890211257882), Sample(-0.002677133553760478),
     Sample(0.0016612171550778898), Sample(-0.0009251472353283694),
     Sample(0.0004346243146011135), Sample(-4.3367201090342294e-05)},
  }};
};

/*
Fir filter specialized to 8x oversampling.

±0.0025 dB ripple at pass band.
-91 dB attenuation at stop band.
For fs=48000Hz, aliasing overwraps to 19000Hz.

```python
import scipy.signal as signal
oversample = 8
signal.remez(192, [0, 0.375, 0.6, oversample], [1, 0], fs=oversample)
```
*/
template<typename Sample> struct FirRemez192Decimation8 {
  constexpr static size_t bufferSize = 24;
  constexpr static size_t nPhase = 8;

  constexpr static std::array<std::array<Sample, bufferSize>, nPhase> coefficient{{
    {Sample(-2.1139354805935512e-05), Sample(5.347371195814565e-05),
     Sample(-0.0001415448461103393),  Sample(0.00030829717079647106),
     Sample(-0.000580667539010022),   Sample(0.0009885913753994082),
     Sample(-0.0015519689985549138),  Sample(0.0022900077235165363),
     Sample(-0.0032296376739828336),  Sample(0.004462800438655286),
     Sample(-0.006353256428726351),   Sample(0.011062021363425968),
     Sample(0.12111945173459678),     Sample(-0.004095637245109295),
     Sample(0.0005924057401233344),   Sample(0.00044754217123647627),
     Sample(-0.0007696998279471991),  Sample(0.0007815761010487296),
     Sample(-0.0006527016443240597),  Sample(0.00047788428784342713),
     Sample(-0.00030950623314961986), Sample(0.0001766722364536085),
     Sample(-8.495512809001416e-05),  Sample(3.318205426713193e-05)},
    {Sample(-1.7311366001983594e-05), Sample(7.022091088116382e-05),
     Sample(-0.00018710813538407),    Sample(0.0004147226606302399),
     Sample(-0.0008035437514774823),  Sample(0.0014178443825536837),
     Sample(-0.0023277165053251373),  Sample(0.0036297226119589463),
     Sample(-0.0054888746482107415),  Sample(0.00830883799285417),
     Sample(-0.013448184393528647),   Sample(0.029237825209169323),
     Sample(0.11517849894514598),     Sample(-0.015360406597174494),
     Sample(0.006497083959285075),    Sample(-0.003165466946011136),
     Sample(0.0015250623104165002),   Sample(-0.0006639794854765167),
     Sample(0.0002277689350273766),   Sample(-3.072343311579114e-05),
     Sample(-3.6145370651806164e-05), Sample(4.326228145123317e-05),
     Sample(-2.7968418213907746e-05), Sample(1.3354230369088297e-05)},
    {Sample(-2.0234722428936636e-05), Sample(7.882040718979042e-05),
     Sample(-0.0002105314022061763),  Sample(0.0004723098250846385),
     Sample(-0.0009323406279398252),  Sample(0.0016829793835441398),
     Sample(-0.0028383198611987726),  Sample(0.004566794646472994),
     Sample(-0.007167999753763226),   Sample(0.011365817967048296),
     Sample(-0.019637100479040544),   Sample(0.0491734740661636),
     Sample(0.10383224002494845),     Sample(-0.022307201803405022),
     Sample(0.010719310924237668),    Sample(-0.005924439955435388),
     Sample(0.0033518778532061547),   Sample(-0.0018500757368017946),
     Sample(0.00096700184746435),     Sample(-0.00046509263455289684),
     Sample(0.00019987354792748983),  Sample(-7.224476173388855e-05),
     Sample(2.086149646528279e-05),   Sample(-2.9270645975825263e-06)},
    {Sample(-1.9586499089068096e-05), Sample(7.504399735135108e-05),
     Sample(-0.00020222268386438486), Sample(0.00046136996856895984),
     Sample(-0.0009295075675721672),  Sample(0.0017161365611929771),
     Sample(-0.0029660267265837256),  Sample(0.004901976840164234),
     Sample(-0.007931518407270974),   Sample(0.013052295608969111),
     Sample(-0.02381756664215646),    Sample(0.06934111066420201),
     Sample(0.08809330442996116),     Sample(-0.024973136133000902),
     Sample(0.012916465866977297),    Sample(-0.0075451905802851),
     Sample(0.004506654813899305),    Sample(-0.002639227374819633),
     Sample(0.0014777182214186606),   Sample(-0.0007735066019763589),
     Sample(0.00037044987276235004),  Sample(-0.0001562418648010693),
     Sample(5.595830993682483e-05),   Sample(-1.4036430091615892e-05)},
    {Sample(-1.4036430091615892e-05), Sample(5.595830993682483e-05),
     Sample(-0.0001562418648010693),  Sample(0.00037044987276235004),
     Sample(-0.0007735066019763589),  Sample(0.0014777182214186606),
     Sample(-0.002639227374819633),   Sample(0.004506654813899305),
     Sample(-0.0075451905802851),     Sample(0.012916465866977297),
     Sample(-0.024973136133000902),   Sample(0.08809330442996116),
     Sample(0.06934111066420201),     Sample(-0.02381756664215646),
     Sample(0.013052295608969111),    Sample(-0.007931518407270974),
     Sample(0.004901976840164234),    Sample(-0.0029660267265837256),
     Sample(0.0017161365611929771),   Sample(-0.0009295075675721672),
     Sample(0.00046136996856895984),  Sample(-0.00020222268386438486),
     Sample(7.504399735135108e-05),   Sample(-1.9586499089068096e-05)},
    {Sample(-2.9270645975825263e-06), Sample(2.086149646528279e-05),
     Sample(-7.224476173388855e-05),  Sample(0.00019987354792748983),
     Sample(-0.00046509263455289684), Sample(0.00096700184746435),
     Sample(-0.0018500757368017946),  Sample(0.0033518778532061547),
     Sample(-0.005924439955435388),   Sample(0.010719310924237668),
     Sample(-0.022307201803405022),   Sample(0.10383224002494845),
     Sample(0.0491734740661636),      Sample(-0.019637100479040544),
     Sample(0.011365817967048296),    Sample(-0.007167999753763226),
     Sample(0.004566794646472994),    Sample(-0.0028383198611987726),
     Sample(0.0016829793835441398),   Sample(-0.0009323406279398252),
     Sample(0.0004723098250846385),   Sample(-0.0002105314022061763),
     Sample(7.882040718979042e-05),   Sample(-2.0234722428936636e-05)},
    {Sample(1.3354230369088297e-05), Sample(-2.7968418213907746e-05),
     Sample(4.326228145123317e-05),  Sample(-3.6145370651806164e-05),
     Sample(-3.072343311579114e-05), Sample(0.0002277689350273766),
     Sample(-0.0006639794854765167), Sample(0.0015250623104165002),
     Sample(-0.003165466946011136),  Sample(0.006497083959285075),
     Sample(-0.015360406597174494),  Sample(0.11517849894514598),
     Sample(0.029237825209169323),   Sample(-0.013448184393528647),
     Sample(0.00830883799285417),    Sample(-0.0054888746482107415),
     Sample(0.0036297226119589463),  Sample(-0.0023277165053251373),
     Sample(0.0014178443825536837),  Sample(-0.0008035437514774823),
     Sample(0.0004147226606302399),  Sample(-0.00018710813538407),
     Sample(7.022091088116382e-05),  Sample(-1.7311366001983594e-05)},
    {Sample(3.318205426713193e-05),  Sample(-8.495512809001416e-05),
     Sample(0.0001766722364536085),  Sample(-0.00030950623314961986),
     Sample(0.00047788428784342713), Sample(-0.0006527016443240597),
     Sample(0.0007815761010487296),  Sample(-0.0007696998279471991),
     Sample(0.00044754217123647627), Sample(0.0005924057401233344),
     Sample(-0.004095637245109295),  Sample(0.12111945173459678),
     Sample(0.011062021363425968),   Sample(-0.006353256428726351),
     Sample(0.004462800438655286),   Sample(-0.0032296376739828336),
     Sample(0.0022900077235165363),  Sample(-0.0015519689985549138),
     Sample(0.0009885913753994082),  Sample(-0.000580667539010022),
     Sample(0.00030829717079647106), Sample(-0.0001415448461103393),
     Sample(5.347371195814565e-05),  Sample(-2.1139354805935512e-05)},
  }};
};

/*
Fir filter specialized to 8x oversampling. Overwrapping.

±0.000014 dB ripple at pass band.
-116 dB attenuation at stop band.
For fs=48000Hz, aliasing overwraps to 19000Hz.

```python
import scipy.signal as signal
oversample = 8
signal.remez(256, [0, 0.375, 0.6, oversample], [1, 0], fs=oversample)
```
*/
template<typename Sample> struct FirRemez256SteepDecimation8 {
  constexpr static size_t bufferSize = 32;
  constexpr static size_t nPhase = 8;

  constexpr static std::array<std::array<Sample, bufferSize>, nPhase> coefficient{{
    {Sample(-1.2270515735399332e-06), Sample(4.978374107758859e-06),
     Sample(-1.692037416445694e-05),  Sample(4.461693315311093e-05),
     Sample(-9.963236300484039e-05),  Sample(0.00019733915552191524),
     Sample(-0.0003558261899578275),  Sample(0.0005941275220175334),
     Sample(-0.0009298972473056768),  Sample(0.0013774964702482737),
     Sample(-0.001948018138817861),   Sample(0.002654168479229834),
     Sample(-0.0035269837550773904),  Sample(0.004668248091400766),
     Sample(-0.006459360318620389),   Sample(0.011079405094818723),
     Sample(0.1211631534929385),      Sample(-0.004163649684579389),
     Sample(0.000646248393177665),    Sample(0.0004380730997247138),
     Sample(-0.0008218136254582062),  Sample(0.0008955216698469529),
     Sample(-0.0008155522625839855),  Sample(0.0006663536396049088),
     Sample(-0.0004993532881890077),  Sample(0.0003453033225176797),
     Sample(-0.0002199272981442123),  Sample(0.00012799940391747657),
     Sample(-6.706149793488577e-05),  Sample(3.082029451446905e-05),
     Sample(-1.1817702181024459e-05), Sample(3.432942861735908e-06)},
    {Sample(-1.0481047341488122e-06), Sample(6.084120811054507e-06),
     Sample(-2.0362265330418985e-05), Sample(5.389302851796453e-05),
     Sample(-0.00012200490669792629), Sample(0.00024648762854703494),
     Sample(-0.0004555744018196346),  Sample(0.0007833702780251432),
     Sample(-0.001269042561728215),   Sample(0.0019573262918314424),
     Sample(-0.002903593059158544),   Sample(0.004191199415042217),
     Sample(-0.005982994790662318),   Sample(0.008688329745651441),
     Sample(-0.01368262744553592),    Sample(0.02931818455364193),
     Sample(0.11523714134727606),     Sample(-0.015524079787076546),
     Sample(0.006719603069873568),    Sample(-0.0033997308646868474),
     Sample(0.0017308508628302282),   Sample(-0.000816056772945318),
     Sample(0.000316209035305098),    Sample(-6.146214353101596e-05),
     Sample(-4.90002396567423e-05),   Sample(7.977637253002407e-05),
     Sample(-7.229203152309606e-05),  Sample(5.154584985599571e-05),
     Sample(-3.0792606343665e-05),    Sample(1.5442043470289392e-05),
     Sample(-6.222837751623365e-06),  Sample(1.832460743063396e-06)},
    {Sample(-1.1896487321356817e-06), Sample(6.3281161650600315e-06),
     Sample(-2.096608390580413e-05),  Sample(5.584165579253361e-05),
     Sample(-0.00012814801023035996), Sample(0.0002634801703688661),
     Sample(-0.0004969208313631156),  Sample(0.0008738148403939522),
     Sample(-0.0014507051041812518),  Sample(0.0022985018452354026),
     Sample(-0.003512753893100307),   Sample(0.005243544598555015),
     Sample(-0.00778272400479136),    Sample(0.011855255036491867),
     Sample(-0.01995250856547176),    Sample(0.04928847352674567),
     Sample(0.1039154601560244),      Sample(-0.022560220787399346),
     Sample(0.011092674847938518),    Sample(-0.006359286174785589),
     Sample(0.003789811720337086),    Sample(-0.0022446287835192866),
     Sample(0.0012876570579126237),   Sample(-0.0007012003037570119),
     Sample(0.00035500612434411627),  Sample(-0.0001625813774172153),
     Sample(6.436232627245083e-05),   Sample(-1.992959559877075e-05),
     Sample(3.2585011857714164e-06),  Sample(1.0767850872025106e-06),
     Sample(-1.1125474540253126e-06), Sample(4.6282415117626426e-07)},
    {Sample(-1.0472095647280878e-06), Sample(5.3275443196532544e-06),
     Sample(-1.7786028288581525e-05), Sample(4.83491316503719e-05),
     Sample(-0.00011367596924483124), Sample(0.00023973161659623105),
     Sample(-0.0004638658455232798),  Sample(0.0008368921880120554),
     Sample(-0.0014257265302872132),  Sample(0.0023189906232642917),
     Sample(-0.003641661209345716),   Sample(0.005595228511457576),
     Sample(-0.008575113961467755),   Sample(0.013575354785558227),
     Sample(-0.024159901473236788),   Sample(0.06946318153423019),
     Sample(0.0882013862773577),      Sample(-0.025291087053468032),
     Sample(0.01339692957245307),     Sample(-0.008124592986248604),
     Sample(0.005115850595613574),    Sample(-0.0032171147744583907),
     Sample(0.0019780899887336814),   Sample(-0.0011719590755385947),
     Sample(0.0006608558903185911),   Sample(-0.0003503111083476137),
     Sample(0.00017210196601439494),  Sample(-7.695554137473143e-05),
     Sample(3.054795954641242e-05),   Sample(-1.0358423633036784e-05),
     Sample(2.8459368774448087e-06),  Sample(-5.135588807356153e-07)},
    {Sample(-5.135588807356153e-07),  Sample(2.8459368774448087e-06),
     Sample(-1.0358423633036784e-05), Sample(3.054795954641242e-05),
     Sample(-7.695554137473143e-05),  Sample(0.00017210196601439494),
     Sample(-0.0003503111083476137),  Sample(0.0006608558903185911),
     Sample(-0.0011719590755385947),  Sample(0.0019780899887336814),
     Sample(-0.0032171147744583907),  Sample(0.005115850595613574),
     Sample(-0.008124592986248604),   Sample(0.01339692957245307),
     Sample(-0.025291087053468032),   Sample(0.0882013862773577),
     Sample(0.06946318153423019),     Sample(-0.024159901473236788),
     Sample(0.013575354785558227),    Sample(-0.008575113961467755),
     Sample(0.005595228511457576),    Sample(-0.003641661209345716),
     Sample(0.0023189906232642917),   Sample(-0.0014257265302872132),
     Sample(0.0008368921880120554),   Sample(-0.0004638658455232798),
     Sample(0.00023973161659623105),  Sample(-0.00011367596924483124),
     Sample(4.83491316503719e-05),    Sample(-1.7786028288581525e-05),
     Sample(5.3275443196532544e-06),  Sample(-1.0472095647280878e-06)},
    {Sample(4.6282415117626426e-07), Sample(-1.1125474540253126e-06),
     Sample(1.0767850872025106e-06), Sample(3.2585011857714164e-06),
     Sample(-1.992959559877075e-05), Sample(6.436232627245083e-05),
     Sample(-0.0001625813774172153), Sample(0.00035500612434411627),
     Sample(-0.0007012003037570119), Sample(0.0012876570579126237),
     Sample(-0.0022446287835192866), Sample(0.003789811720337086),
     Sample(-0.006359286174785589),  Sample(0.011092674847938518),
     Sample(-0.022560220787399346),  Sample(0.1039154601560244),
     Sample(0.04928847352674567),    Sample(-0.01995250856547176),
     Sample(0.011855255036491867),   Sample(-0.00778272400479136),
     Sample(0.005243544598555015),   Sample(-0.003512753893100307),
     Sample(0.0022985018452354026),  Sample(-0.0014507051041812518),
     Sample(0.0008738148403939522),  Sample(-0.0004969208313631156),
     Sample(0.0002634801703688661),  Sample(-0.00012814801023035996),
     Sample(5.584165579253361e-05),  Sample(-2.096608390580413e-05),
     Sample(6.3281161650600315e-06), Sample(-1.1896487321356817e-06)},
    {Sample(1.832460743063396e-06),  Sample(-6.222837751623365e-06),
     Sample(1.5442043470289392e-05), Sample(-3.0792606343665e-05),
     Sample(5.154584985599571e-05),  Sample(-7.229203152309606e-05),
     Sample(7.977637253002407e-05),  Sample(-4.90002396567423e-05),
     Sample(-6.146214353101596e-05), Sample(0.000316209035305098),
     Sample(-0.000816056772945318),  Sample(0.0017308508628302282),
     Sample(-0.0033997308646868474), Sample(0.006719603069873568),
     Sample(-0.015524079787076546),  Sample(0.11523714134727606),
     Sample(0.02931818455364193),    Sample(-0.01368262744553592),
     Sample(0.008688329745651441),   Sample(-0.005982994790662318),
     Sample(0.004191199415042217),   Sample(-0.002903593059158544),
     Sample(0.0019573262918314424),  Sample(-0.001269042561728215),
     Sample(0.0007833702780251432),  Sample(-0.0004555744018196346),
     Sample(0.00024648762854703494), Sample(-0.00012200490669792629),
     Sample(5.389302851796453e-05),  Sample(-2.0362265330418985e-05),
     Sample(6.084120811054507e-06),  Sample(-1.0481047341488122e-06)},
    {Sample(3.432942861735908e-06),  Sample(-1.1817702181024459e-05),
     Sample(3.082029451446905e-05),  Sample(-6.706149793488577e-05),
     Sample(0.00012799940391747657), Sample(-0.0002199272981442123),
     Sample(0.0003453033225176797),  Sample(-0.0004993532881890077),
     Sample(0.0006663536396049088),  Sample(-0.0008155522625839855),
     Sample(0.0008955216698469529),  Sample(-0.0008218136254582062),
     Sample(0.0004380730997247138),  Sample(0.000646248393177665),
     Sample(-0.004163649684579389),  Sample(0.1211631534929385),
     Sample(0.011079405094818723),   Sample(-0.006459360318620389),
     Sample(0.004668248091400766),   Sample(-0.0035269837550773904),
     Sample(0.002654168479229834),   Sample(-0.001948018138817861),
     Sample(0.0013774964702482737),  Sample(-0.0009298972473056768),
     Sample(0.0005941275220175334),  Sample(-0.0003558261899578275),
     Sample(0.00019733915552191524), Sample(-9.963236300484039e-05),
     Sample(4.461693315311093e-05),  Sample(-1.692037416445694e-05),
     Sample(4.978374107758859e-06),  Sample(-1.2270515735399332e-06)},
  }};
};

/*
Fir filter specialized to 8x oversampling. Non-overwrapping.

+-0.00021 dB ripple at pass band.
-92 dB attenuation at stop band.
No aliasing. Constant -92 dB noise floor.

```python
import scipy.signal as signal
oversample = 8
signal.remez(256, [0, 0.325, 0.5, oversample], [1, 0], fs=oversample)
```
*/
template<typename Sample> struct FirRemez256GentleDecimation8 {
  constexpr static size_t bufferSize = 32;
  constexpr static size_t nPhase = 8;

  constexpr static std::array<std::array<Sample, bufferSize>, nPhase> coefficient{{
    {Sample(-1.1993081333552215e-05), Sample(1.964351220155972e-05),
     Sample(-5.8542034232794804e-05), Sample(0.00010729638207563003),
     Sample(-0.00010861906408281117), Sample(-3.867915381346443e-05),
     Sample(0.0004322331823420826),   Sample(-0.0010693301612759568),
     Sample(0.0017399991433914232),   Sample(-0.0019825991186250935),
     Sample(0.0011686283811654897),   Sample(0.001274210241157238),
     Sample(-0.005583154628339127),   Sample(0.011492044293419247),
     Sample(-0.018433541682628757),   Sample(0.027315874351768445),
     Sample(0.10270431736879226),     Sample(0.013879731900139584),
     Sample(-0.01452736311674651),    Sample(0.011166138134262784),
     Sample(-0.00674755211588896),    Sample(0.002769935016919955),
     Sample(-4.009490587425957e-05),  Sample(-0.001274898178295907),
     Sample(0.0014817554685805524),   Sample(-0.0010903987836278143),
     Sample(0.0005610052774460024),   Sample(-0.00016371848490421748),
     Sample(-3.021698523771811e-05),  Sample(7.299171143099649e-05),
     Sample(-4.9361141768451935e-05), Sample(1.928471889670624e-05)},
    {Sample(3.2249589769729163e-07),  Sample(1.7250259436960175e-05),
     Sample(-6.21865122856873e-05),   Sample(0.0001348810486536517),
     Sample(-0.00018690253193160978), Sample(0.00010897334152472398),
     Sample(0.000237006555520949),    Sample(-0.0009240014134695712),
     Sample(0.001829639950691649),    Sample(-0.0025392077747600893),
     Sample(0.002363175865517446),    Sample(-0.0005054495164705803),
     Sample(-0.003682261303252013),   Sample(0.010571294285325182),
     Sample(-0.02084932000189474),    Sample(0.041675754428785475),
     Sample(0.09909143706061187),     Sample(0.002054354856634085),
     Sample(-0.009687967298185867),   Sample(0.009752192980934491),
     Sample(-0.0071251831835546856),  Sample(0.003854818212075314),
     Sample(-0.0011385633017672824),  Sample(-0.000502278062690789),
     Sample(0.0010977093770892592),   Sample(-0.0009979787446873758),
     Sample(0.0006175383360506621),   Sample(-0.0002559487752727118),
     Sample(4.000172807136199e-05),   Sample(3.672733065133769e-05),
     Sample(-3.665115799920039e-05),  Sample(1.6870248227605994e-05)},
    {Sample(2.2503681738789036e-06),  Sample(1.1649917134935601e-05),
     Sample(-5.8513126534126e-05),    Sample(0.00015097289101137599),
     Sample(-0.00025549218821049695), Sample(0.00026487066856409976),
     Sample(-1.0320382994734447e-05), Sample(-0.0006568375452585839),
     Sample(0.0017226463034250854),   Sample(-0.00286712707082035),
     Sample(0.0034089105189780726),   Sample(-0.002398837725595235),
     Sample(-0.0011802716309675257),  Sample(0.008359540601009104),
     Sample(-0.02128054679171583),    Sample(0.056165953198831345),
     Sample(0.09209852583342776),     Sample(-0.007625313358707562),
     Sample(-0.004485531288872734),   Sample(0.007497632466445405),
     Sample(-0.006748480704862311),   Sample(0.004454423163837809),
     Sample(-0.0020244163573003857),  Sample(0.0002502048197780836),
     Sample(0.0006395888736872797),   Sample(-0.0008135421619874107),
     Sample(0.0006042820181490843),   Sample(-0.00030991397078378364),
     Sample(9.581097375344377e-05),   Sample(2.7367287878378366e-06),
     Sample(-2.2556287813925137e-05), Sample(1.3220180164909415e-05)},
    {Sample(5.300322984843467e-06),   Sample(2.784425972711471e-06),
     Sample(-4.64128664581005e-05),   Sample(0.00015143395461486247),
     Sample(-0.00030449440307228257), Sample(0.00041167081465447847),
     Sample(-0.00028746546287128664), Sample(-0.00028463844665736407),
     Sample(0.0014105295845201326),   Sample(-0.002906953178611909),
     Sample(0.0041751036972818185),   Sample(-0.004207170906490688),
     Sample(0.0017087101235156495),   Sample(0.004943184372771973),
     Sample(-0.01934395014863162),    Sample(0.069942604852471),
     Sample(0.08217310732294353),     Sample(-0.01480836904577092),
     Sample(0.000544552019992989),    Sample(0.004707184186562935),
     Sample(-0.005722272775313291),   Sample(0.004549876432132227),
     Sample(-0.002625974909866911),   Sample(0.0009073076371669507),
     Sample(0.00016176855322777128),  Sample(-0.0005659165925500945),
     Sample(0.0005307435464635439),   Sample(-0.00032491953302187285),
     Sample(0.000133383985341516),    Sample(-2.569079832105044e-05),
     Sample(-8.933338715263757e-06),  Sample(9.130442103556708e-06)},
    {Sample(9.130442103556708e-06),   Sample(-8.933338715263757e-06),
     Sample(-2.569079832105044e-05),  Sample(0.000133383985341516),
     Sample(-0.00032491953302187285), Sample(0.0005307435464635439),
     Sample(-0.0005659165925500945),  Sample(0.00016176855322777128),
     Sample(0.0009073076371669507),   Sample(-0.002625974909866911),
     Sample(0.004549876432132227),    Sample(-0.005722272775313291),
     Sample(0.004707184186562935),    Sample(0.000544552019992989),
     Sample(-0.01480836904577092),    Sample(0.08217310732294353),
     Sample(0.069942604852471),       Sample(-0.01934395014863162),
     Sample(0.004943184372771973),    Sample(0.0017087101235156495),
     Sample(-0.004207170906490688),   Sample(0.0041751036972818185),
     Sample(-0.002906953178611909),   Sample(0.0014105295845201326),
     Sample(-0.00028463844665736407), Sample(-0.00028746546287128664),
     Sample(0.00041167081465447847),  Sample(-0.00030449440307228257),
     Sample(0.00015143395461486247),  Sample(-4.64128664581005e-05),
     Sample(2.784425972711471e-06),   Sample(5.300322984843467e-06)},
    {Sample(1.3220180164909415e-05),  Sample(-2.2556287813925137e-05),
     Sample(2.7367287878378366e-06),  Sample(9.581097375344377e-05),
     Sample(-0.00030991397078378364), Sample(0.0006042820181490843),
     Sample(-0.0008135421619874107),  Sample(0.0006395888736872797),
     Sample(0.0002502048197780836),   Sample(-0.0020244163573003857),
     Sample(0.004454423163837809),    Sample(-0.006748480704862311),
     Sample(0.007497632466445405),    Sample(-0.004485531288872734),
     Sample(-0.007625313358707562),   Sample(0.09209852583342776),
     Sample(0.056165953198831345),    Sample(-0.02128054679171583),
     Sample(0.008359540601009104),    Sample(-0.0011802716309675257),
     Sample(-0.002398837725595235),   Sample(0.0034089105189780726),
     Sample(-0.00286712707082035),    Sample(0.0017226463034250854),
     Sample(-0.0006568375452585839),  Sample(-1.0320382994734447e-05),
     Sample(0.00026487066856409976),  Sample(-0.00025549218821049695),
     Sample(0.00015097289101137599),  Sample(-5.8513126534126e-05),
     Sample(1.1649917134935601e-05),  Sample(2.2503681738789036e-06)},
    {Sample(1.6870248227605994e-05), Sample(-3.665115799920039e-05),
     Sample(3.672733065133769e-05),  Sample(4.000172807136199e-05),
     Sample(-0.0002559487752727118), Sample(0.0006175383360506621),
     Sample(-0.0009979787446873758), Sample(0.0010977093770892592),
     Sample(-0.000502278062690789),  Sample(-0.0011385633017672824),
     Sample(0.003854818212075314),   Sample(-0.0071251831835546856),
     Sample(0.009752192980934491),   Sample(-0.009687967298185867),
     Sample(0.002054354856634085),   Sample(0.09909143706061187),
     Sample(0.041675754428785475),   Sample(-0.02084932000189474),
     Sample(0.010571294285325182),   Sample(-0.003682261303252013),
     Sample(-0.0005054495164705803), Sample(0.002363175865517446),
     Sample(-0.0025392077747600893), Sample(0.001829639950691649),
     Sample(-0.0009240014134695712), Sample(0.000237006555520949),
     Sample(0.00010897334152472398), Sample(-0.00018690253193160978),
     Sample(0.0001348810486536517),  Sample(-6.21865122856873e-05),
     Sample(1.7250259436960175e-05), Sample(3.2249589769729163e-07)},
    {Sample(1.928471889670624e-05),   Sample(-4.9361141768451935e-05),
     Sample(7.299171143099649e-05),   Sample(-3.021698523771811e-05),
     Sample(-0.00016371848490421748), Sample(0.0005610052774460024),
     Sample(-0.0010903987836278143),  Sample(0.0014817554685805524),
     Sample(-0.001274898178295907),   Sample(-4.009490587425957e-05),
     Sample(0.002769935016919955),    Sample(-0.00674755211588896),
     Sample(0.011166138134262784),    Sample(-0.01452736311674651),
     Sample(0.013879731900139584),    Sample(0.10270431736879226),
     Sample(0.027315874351768445),    Sample(-0.018433541682628757),
     Sample(0.011492044293419247),    Sample(-0.005583154628339127),
     Sample(0.001274210241157238),    Sample(0.0011686283811654897),
     Sample(-0.0019825991186250935),  Sample(0.0017399991433914232),
     Sample(-0.0010693301612759568),  Sample(0.0004322331823420826),
     Sample(-3.867915381346443e-05),  Sample(-0.00010861906408281117),
     Sample(0.00010729638207563003),  Sample(-5.8542034232794804e-05),
     Sample(1.964351220155972e-05),   Sample(-1.1993081333552215e-05)},
  }};
};

template<typename Sample> struct FirTest {
  constexpr static size_t bufferSize = 4;
  constexpr static size_t nPhase = 2;

  constexpr static std::array<std::array<Sample, bufferSize>, nPhase> coefficient{{
    {Sample(0), Sample(1), Sample(2), Sample(3)},
    {Sample(4), Sample(5), Sample(6), Sample(7)},
  }};
};

template<typename Sample, typename FIR> class FirFilterRot {
  std::array<std::array<Sample, FIR::bufferSize>, FIR::nPhase> buf;
  std::
    array<std::array<std::array<Sample, FIR::bufferSize>, FIR::bufferSize>, FIR::nPhase>
      coefficient;

  uint8_t bufptr = 0;

public:
  FirFilterRot()
  {
    for (uint8_t ph = 0; ph < FIR::nPhase; ++ph) {
      const auto &src = FIR::coefficient[ph];
      std::array<Sample, FIR::bufferSize> fir;
      std::reverse_copy(src.begin(), src.end(), fir.begin());

      for (uint8_t rot = 0; rot < FIR::bufferSize; ++rot) {
        std::rotate(fir.rbegin(), fir.rbegin() + 1, fir.rend());
        coefficient[ph][rot] = fir;
      }
    }

    // // debug: FirFilterRot で以下が出力される。
    // //
    // // 0, 3, 2, 1,
    // // 1, 0, 3, 2,
    // // 2, 1, 0, 3,
    // // 3, 2, 1, 0,
    // //
    // // 4, 7, 6, 5,
    // // 5, 4, 7, 6,
    // // 6, 5, 4, 7,
    // // 7, 6, 5, 4,
    // for (const auto &ph : co) {
    //   for (const auto &rot : ph) {
    //     for (const auto &val : rot) std::cout << val << ", ";
    //     std::cout << "\n";
    //   }
    //   std::cout << "\n";
    // }

    reset();
  }

  void reset()
  {
    for (auto &bf : buf) bf.fill(Sample(0));
  }

  Sample process(std::array<Sample, FIR::nPhase> &input)
  {
    Sample sum = 0;

    for (size_t ph = 0; ph < FIR::nPhase; ++ph) {
      auto &bf = buf[ph];
      bf[bufptr] = input[ph];

      const auto &fir = coefficient[ph][bufptr];
      for (size_t i = 0; i < fir.size(); ++i) sum += bf[i] * fir[i];
    }

    ++bufptr;
    bufptr = (bufptr + 1) % FIR::bufferSize;

    return sum;
  }
};

template<typename Sample, typename FIR> class FirFilter {
  std::array<std::array<Sample, FIR::bufferSize>, FIR::nPhase> buf;

public:
  FirFilter() { reset(); }

  void reset()
  {
    for (auto &bf : buf) bf.fill(Sample(0));
  }

  Sample process(std::array<Sample, FIR::nPhase> &input)
  {
    Sample sum = 0;

    for (size_t ph = 0; ph < FIR::nPhase; ++ph) {
      auto &bf = buf[ph];
      for (size_t i = 0; i < bf.size() - 1; ++i) bf[i] = bf[i + 1];
      bf.back() = input[ph];

      const auto &fir = FIR::coefficient[ph];
      for (size_t i = 0; i < fir.size(); ++i) sum += bf[i] * fir[i];
    }

    return sum;
  }
};

template<typename Sample>
std::vector<Sample> generateSaw(size_t samples, Sample samplerate, Sample frequency)
{
  std::vector<Sample> buf;
  buf.reserve(samples);

  Sample delta = frequency / samplerate;
  Sample phase = 0;
  for (size_t i = 0; i < samples; ++i) {
    buf.push_back(2.0f * phase - 1.0f);
    phase += delta;
    phase -= std::floor(phase);
  }

  return buf;
}

template<typename Sample, size_t oversample> struct Bypass {
  Sample process(std::array<Sample, oversample> &input) { return input[oversample - 1]; }
};

constexpr size_t samplerate = 48000;
constexpr size_t nLoop = 128;

template<size_t oversample, typename Lowpass> void bench(std::string name)
{
  // float freq = 440.0f * std::pow(2.0f, (84.0f - 69.0f) / 12.0f);
  // auto signal
  //   = generateSaw<float>(oversample * samplerate, oversample * samplerate, freq);

  std::string path = "sample/saw" + std::to_string(oversample) + ".wav";
  SoundFile snd(path);
  auto &signal = snd.data;

  std::vector<float> wav(signal.size() / oversample);

  Lowpass lp;
  std::array<float, oversample> buffer;

  double sumElapsed = 0.0;
  for (size_t loop = 0; loop < nLoop; ++loop) {
    for (size_t idx = 0; idx < wav.size(); ++idx) {
      auto start = std::chrono::steady_clock::now();

      for (size_t os = 0; os < oversample; ++os)
        buffer[os] = signal[idx * oversample + os];

      wav[idx] = lp.process(buffer);

      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<double, std::milli> elapsed = end - start;
      sumElapsed += elapsed.count();
    }
  }
  std::cout << name << ": " << signal.size() << "[sample], " << sumElapsed << "[ms], "
            << "\n";

  name = "snd/out_" + std::to_string(oversample) + "x_" + name + ".wav";
  writeWave(name, wav, samplerate);
}

int main()
{
  std::cout << "\nDry Run ----------------\n\n";
  bench<16, Bypass<float, 16>>("Naive");

  std::cout << "\nBenchmark ----------------\n\n";

  // bench<8, SosFilter<float, SosButterDecimation8<float>>>("Butter8");
  // bench<8, SosFilter<float, SosEllipticDecimation8<float>>>("Elliptic8");
  // bench<8, SosFilterDF2<float, SosButterDecimation8<float>>>("Butter8_DF2");
  // bench<8, SosFilterDF2<float, SosEllipticDecimation8<float>>>("Elliptic8_DF2");
  // bench<8, FirFilter<float, FirRemez128Decimation8<float>>>("Remez128_8");
  // bench<8, FirFilter<float, FirRemez192Decimation8<float>>>("Remez192_8");
  // bench<8, FirFilter<float, FirRemez256SteepDecimation8<float>>>("Remez256Steep_8");
  // bench<8, FirFilter<float, FirRemez256GentleDecimation8<float>>>("Remez256Gentle_8");
  // bench<8, FirFilterRot<float, FirRemez128Decimation8<float>>>("RotRemez128_8");
  // bench<8, FirFilterRot<float, FirRemez192Decimation8<float>>>("RotRemez192_8");
  // bench<8, FirFilterRot<float, FirRemez256GentleDecimation8<float>>>(
  //   "RotRemez256Gentle_8");
  bench<16, SosFilter<float, SosEllipticDecimation16<float>>>("Elliptic16");
  bench<16, SosFilterLatency<float, SosEllipticDecimation16<float>>>(
    "Elliptic16_Latency");
  bench<16, SosFilterDF2<float, SosEllipticDecimation16<float>>>("Elliptic16_DF2");
  bench<16, SosFilterDF2Latency<float, SosEllipticDecimation16<float>>>(
    "Elliptic16_DF2_Latency");
  return 0;
}
