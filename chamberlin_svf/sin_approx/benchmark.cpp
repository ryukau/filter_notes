#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

namespace approx {

// Fast approximation of 2 * sin(pi * x) for x in [0.0f, 0.5f]
inline float sin_pi_x2(float x) noexcept {
  const float t = x * x;
  float res = 0.15512077399924159f;
  res = std::fma(res, t, -1.1964842528437377f);
  res = std::fma(res, t, 5.100139452670464f);
  res = std::fma(res, t, -10.335419369597881f);
  res = std::fma(res, t, 6.283185273790786f);
  return x * res;
}

// Fast approximation of 2 * sin(pi * x) for x in [0.0, 0.5]
inline double sin_pi_x2(double x) noexcept {
  const double t = x * x;
  double res = -0.000042245905984219777;
  res = std::fma(res, t, 0.0009319552340645613);
  res = std::fma(res, t, -0.014740722143372288);
  res = std::fma(res, t, 0.16429175651310609);
  res = std::fma(res, t, -1.1985290575514195);
  res = std::fma(res, t, 5.100328079719553);
  res = std::fma(res, t, -10.335425560099505);
  res = std::fma(res, t, 6.283185307179586);
  return x * res;
}

} // namespace approx

// Helper to prevent compiler optimization on variables
template<typename T> void do_not_optimize(T const& val) {
  // asm volatile("" : : "g"(&val) : "memory");
  std::cout << val << "\n";
}

int main() {
  const int N = 10000000;
  std::vector<float> input_f(N);
  std::vector<double> input_d(N);

  std::mt19937 gen(42);
  std::uniform_real_distribution<double> dis(0.0, 0.5);
  for (int i = 0; i < N; ++i) {
    double val = dis(gen);
    input_f[i] = static_cast<float>(val);
    input_d[i] = val;
  }

  // Warmup
  double dummy = 0;
  for (int i = 0; i < 100000; ++i) { dummy += approx::sin_pi_x2(input_f[i % N]); }
  do_not_optimize(dummy);

  // --- FLOAT BENCHMARKS ---
  {
    auto start = std::chrono::high_resolution_clock::now();
    float sum = 0.0f;
    for (int i = 0; i < N; ++i) { sum += 2.0f * std::sin(3.141592653589793f * input_f[i]); }
    do_not_optimize(sum);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> diff = end - start;
    std::cout << "Float std::sin (2*sin(pi*x)): " << diff.count() << " ms\n";
  }

  {
    auto start = std::chrono::high_resolution_clock::now();
    float sum = 0.0f;
    for (int i = 0; i < N; ++i) { sum += approx::sin_pi_x2(input_f[i]); }
    do_not_optimize(sum);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> diff = end - start;
    std::cout << "Float approx::sin_pi_x2:   " << diff.count() << " ms\n";
  }

  // --- DOUBLE BENCHMARKS ---
  {
    auto start = std::chrono::high_resolution_clock::now();
    double sum = 0.0;
    for (int i = 0; i < N; ++i) { sum += 2.0 * std::sin(3.141592653589793 * input_d[i]); }
    do_not_optimize(sum);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> diff = end - start;
    std::cout << "Double std::sin (2*sin(pi*x)): " << diff.count() << " ms\n";
  }

  {
    auto start = std::chrono::high_resolution_clock::now();
    double sum = 0.0;
    for (int i = 0; i < N; ++i) { sum += approx::sin_pi_x2(input_d[i]); }
    do_not_optimize(sum);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> diff = end - start;
    std::cout << "Double approx::sin_pi_x2:   " << diff.count() << " ms\n";
  }

  return 0;
}
