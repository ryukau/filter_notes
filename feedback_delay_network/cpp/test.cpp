#include "lib/pcg-cpp/pcg_random.hpp"
#include <array>
#include <chrono>
#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <string>

/**
At least on older version of `cl.exe`, the size of `matrix` may be truncated if `length`
is too long.
*/
template<typename Sample, size_t length> struct FeedbackMatrix {
  std::array<std::array<Sample, length>, length> matrix{};

  enum class MatrixType {
    orthogonal,
    specialOrthogonal,
    householder,
    circulant,
    upperTriangular,
    lowerTriangular,
    schroeder,
    absorbent,
    hadamard,
    conference,
  };

  template<size_t dim>
  void randomOrthogonal(unsigned seed, std::array<std::array<Sample, dim>, dim> &H)
  {
    pcg64 rng{};
    rng.seed(seed);
    std::normal_distribution<Sample> dist{}; // mean 0, stddev 1.

    H.fill({});
    for (size_t i = 0; i < dim; ++i) H[i][i] = Sample(1);

    std::array<Sample, dim> x;
    for (size_t n = 0; n < dim; ++n) {
      auto xRange = dim - n;
      for (size_t i = 0; i < xRange; ++i) x[i] = dist(rng);

      Sample norm2 = 0;
      for (size_t i = 0; i < xRange; ++i) norm2 += x[i] * x[i];

      Sample x0 = x[0];

      Sample D = x0 >= 0 ? Sample(1) : Sample(-1);
      x[0] += D * std::sqrt(norm2);

      Sample denom = std::sqrt((norm2 - x0 * x0 + x[0] * x[0]) / Sample(2));
      for (size_t i = 0; i < xRange; ++i) x[i] /= denom;

      for (size_t row = 0; row < dim; ++row) {
        Sample dotH = 0;
        for (size_t col = 0; col < xRange; ++col) dotH += H[col][row] * x[col];
        for (size_t col = 0; col < xRange; ++col) {
          H[col][row] = D * (H[col][row] - dotH * x[col]);
        }
      }
    }
  }

  /**
  Randomize `matrix` as special orthogonal matrix. The algorithm is ported from
  `scipy.stats.special_ortho_group` in SciPy v1.8.0.
  */
  template<size_t dim>
  void randomSpecialOrthogonal(unsigned seed, std::array<std::array<Sample, dim>, dim> &H)
  {
    pcg64 rng{};
    rng.seed(seed);
    std::normal_distribution<Sample> dist{}; // mean 0, stddev 1.

    H.fill({});
    for (size_t i = 0; i < dim; ++i) H[i][i] = Sample(1);

    std::array<Sample, dim> x;
    std::array<Sample, dim> D;
    for (size_t n = 0; n < dim; ++n) {
      auto xRange = dim - n;
      for (size_t i = 0; i < xRange; ++i) x[i] = dist(rng);

      Sample norm2 = 0;
      for (size_t i = 0; i < xRange; ++i) norm2 += x[i] * x[i];

      Sample x0 = x[0];

      D[n] = x0 >= 0 ? Sample(1) : Sample(-1);
      x[0] += D[n] * std::sqrt(norm2);

      Sample denom = std::sqrt((norm2 - x0 * x0 + x[0] * x[0]) / Sample(2));
      for (size_t i = 0; i < xRange; ++i) x[i] /= denom;

      for (size_t row = 0; row < dim; ++row) {
        Sample dotH = 0;
        for (size_t col = 0; col < xRange; ++col) dotH += H[col][row] * x[col];
        for (size_t col = 0; col < xRange; ++col) H[col][row] -= dotH * x[col];
      }
    }

    size_t back = dim - 1;
    D[back] = (back & 0b1) == 0 ? Sample(1) : Sample(-1);
    for (size_t i = 0; i < back; ++i) D[back] *= D[i];

    for (size_t row = 0; row < dim; ++row) {
      for (size_t col = 0; col < dim; ++col) H[col][row] *= D[row];
    }
  }

  /**
  Reference: https://nhigham.com/2020/09/15/what-is-a-householder-matrix/
  */
  template<size_t dim>
  void randomHouseholder(unsigned seed, std::array<std::array<Sample, dim>, dim> &matrix)
  {
    pcg64 rng{};
    rng.seed(seed);
    std::uniform_real_distribution<Sample> dist{Sample(0), Sample(1)};

    std::array<Sample, dim> vec{};
    for (size_t i = 0; i < dim; ++i) vec[i] = dist(rng);

    Sample denom = 0;
    for (size_t i = 0; i < dim; ++i) denom += vec[i] * vec[i];

    // Return identity matrix if `vec` is all 0.
    if (denom <= std::numeric_limits<Sample>::epsilon()) {
      for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
          matrix[i][j] = i == j ? Sample(1) : Sample(0);
        }
      }
      return;
    }

    auto scale = Sample(-2) / denom;

    for (size_t i = 0; i < dim; ++i) {
      // Diagonal elements.
      matrix[i][i] = Sample(1) + scale * vec[i] * vec[i];

      // Non-diagonal elements.
      for (size_t j = i + 1; j < dim; ++j) {
        auto value = scale * vec[i] * vec[j];
        matrix[i][j] = value;
        matrix[j][i] = value;
      }
    }
  }

  template<size_t dim>
  void randomOrthogonalCirculant(
    unsigned seed, std::array<std::array<Sample, dim>, dim> &matrix)
  {
    pcg64 rng{};
    rng.seed(seed);
    std::uniform_real_distribution<Sample> dist{Sample(0), Sample(1)};

    std::array<Sample, length> source;
    Sample sum = 0;
    do {
      sum = 0;
      for (auto &value : source) {
        value = dist(rng);
        sum += value;
      }
    } while (sum == 0); // Avoid 0 division.

    Sample scale = Sample(2) / sum;

    std::array<Sample, length> squared;
    for (size_t i = 0; i < length; ++i) squared[i] = std::sqrt(source[i]);

    for (size_t row = 0; row < length; ++row) {
      for (size_t col = 0; col < length; ++col) {
        matrix[row][col] = row == col ? scale * source[row] - Sample(1)
                                      : scale * squared[row] * squared[col];
      }
    }
  }

  template<size_t dim>
  void
  randomUpperTriangular(unsigned seed, std::array<std::array<Sample, dim>, dim> &matrix)
  {
    pcg64 rng{};
    rng.seed(seed);
    // TODO: Add option for [-1, 0] and [0, 1].
    std::uniform_real_distribution<Sample> dist{Sample(-1), Sample(0)};

    matrix.fill({});

    for (size_t row = 0; row < length; ++row) {
      for (size_t col = row; col < length; ++col) matrix[row][col] = dist(rng);
    }
    for (size_t col = 0; col < length; ++col) {
      Sample sum = 0;
      for (size_t row = 0; row < col + 1; ++row) sum += matrix[row][col];
      Sample scale = Sample(2) / sum;
      matrix[col][col] = scale * matrix[col][col] - Sample(1);
      for (size_t row = 0; row < col; ++row) matrix[row][col] *= scale;
    }
  }

  template<size_t dim>
  void
  randomLowerTriangular(unsigned seed, std::array<std::array<Sample, dim>, dim> &matrix)
  {
    pcg64 rng{};
    rng.seed(seed);
    std::uniform_real_distribution<Sample> dist{Sample(-1), Sample(0)};

    matrix.fill({});

    for (size_t row = 0; row < length; ++row) {
      for (size_t col = 0; col < row + 1; ++col) matrix[row][col] = dist(rng);
    }
    for (size_t col = 0; col < length; ++col) {
      Sample sum = 0;
      for (size_t row = col; row < length; ++row) sum += matrix[row][col];
      Sample scale = Sample(2) / sum;
      matrix[col][col] = scale * matrix[col][col] - Sample(1);
      for (size_t row = col + 1; row < length; ++row) matrix[row][col] *= scale;
    }
  }

  template<size_t dim>
  void randomSchroeder(unsigned seed, std::array<std::array<Sample, dim>, dim> &matrix)
  {
    static_assert(
      length >= 2, "FeedbackDelayNetwork::randomSchroeder(): length must be >= 2.");

    pcg64 rng{};
    rng.seed(seed);
    std::uniform_real_distribution<Sample> dist{Sample(0), Sample(1)};

    matrix.fill({});

    for (size_t idx = 0; idx < length; ++idx) matrix[idx][idx] = dist(rng);

    auto &&parallelGain = matrix[length - 2][length - 2];
    for (size_t col = 0; col < length - 2; ++col) {
      matrix[length - 2][col] = Sample(1);
      matrix[length - 1][col] = -parallelGain;
    }
    matrix[length - 1][length - 2] = Sample(1) - parallelGain * parallelGain;
  }

  /**
  Construct following matrix:

  ```
  [[-A * G  , A ],
   [ I - G^2, G ]]
  ```

  - I is identity matrix.
  - G is diagonal matrix represents all-pass gain. diag(g1, g2, ...).
  - A is orthogonal matrix.
  */
  template<size_t dim>
  void randomAbsorbent(unsigned seed, std::array<std::array<Sample, dim>, dim> &matrix)
  {
    static_assert(
      length >= 2, "FeedbackDelayNetwork::randomAbsorbent(): length must be >= 2.");
    static_assert(
      length % 2 == 0, "FeedbackDelayNetwork::randomAbsorbent(): length must be even.");

    pcg64 rng{};
    rng.seed(seed);
    std::uniform_real_distribution<Sample> dist{Sample(-1), Sample(0)};
    std::uniform_int_distribution<unsigned> seeder{
      0, std::numeric_limits<unsigned>::max()};

    constexpr size_t half = length / 2;

    matrix.fill({});

    std::array<std::array<Sample, half>, half> A;
    randomOrthogonal(seeder(rng), A);

    for (size_t col = 0; col < half; ++col) {
      auto gain = dist(rng);
      matrix[half + col][half + col] = gain;             // Fill lower right.
      matrix[half + col][col] = Sample(1) - gain * gain; // Fill lower left.
      for (size_t row = 0; row < half; ++row) {
        matrix[row][half + col] = A[row][col];  // Fill top right.
        matrix[row][col] = -A[row][col] * gain; // Fill top left
      }
    }
  }

  /** Sylvester's construction of Hadamard matrix. */
  template<size_t dim>
  void constructHadamardSylvester(std::array<std::array<Sample, dim>, dim> &mat)
  {
    // This static_assert condition is obtained from: https://stackoverflow.com/a/19399478
    static_assert(
      dim && ((dim & (dim - 1)) == 0),
      "FeedbackDelayNetwork::constructHadamardSylvester(): dim must be power of 2.");

    mat[0][0] = Sample(1) / std::sqrt(Sample(dim));

    size_t start = 1;
    size_t end = 2;
    while (start < dim) {
      for (size_t row = start; row < end; ++row) {
        for (size_t col = start; col < end; ++col) {
          auto value = mat[row - start][col - start];
          mat[row - start][col] = value; // Upper right.
          mat[row][col - start] = value; // Lower left.
          mat[row][col] = -value;        // Lower right.
        }
      }
      start *= 2;
      end *= 2;
    }
  }

  /**
  Construct a kind of conference matrix which is used in Paley's construction of Hadamard
  matrix.

  `dim` must follow the conditions below:
  - `dim mod 4 == 2`.
  - `dim - 1` is sum of 2 squared integer.

  This implementation use the sequence from https://oeis.org/A000952 to determine the size
  of matrix. It's possible to construct this kind of conference matrix greater than size
  of 62, but they are out of scope of FDN64Reverb.
  */
  template<size_t dim>
  void constructConference(std::array<std::array<Sample, dim>, dim> &mat)
  {
    constexpr std::array<size_t, 13> candidates{
      62, 54, 50, 46, 42, 38, 30, 26, 18, 14, 10, 6, 2,
    };

    auto found = std::find_if(
      candidates.begin(), candidates.end(), [](size_t size) { return size <= dim; });
    if (found == candidates.end()) return; // mat is too small.

    size_t dimension = *found;
    size_t modulo = dimension - 1;

    std::set<size_t> quadraticResidue;
    for (size_t i = 1; i < modulo; ++i) quadraticResidue.emplace((i * i) % modulo);
    quadraticResidue.erase(0); // Just in case.

    Sample value = Sample(1) / std::sqrt(Sample(modulo));
    std::vector<Sample> symbol; // Legendre symbol of quadratic residue.
    symbol.reserve(modulo);
    symbol.push_back(0);
    for (size_t i = 1; i < modulo; ++i) {
      symbol.push_back(quadraticResidue.count(i) ? value : -value);
    }

    mat.fill({});
    mat[0][0] = 0;
    for (size_t i = 1; i < dimension; ++i) {
      mat[0][i] = value;
      mat[i][0] = value;
    }
    for (size_t row = 1; row < dimension; ++row) {
      std::copy(symbol.begin(), symbol.end(), mat[row].begin() + 1);
      std::rotate(symbol.rbegin(), symbol.rbegin() + 1, symbol.rend());
    }
  }

  void randomizeMatrix(MatrixType matrixType, unsigned seed = 0)
  {
    if (matrixType == MatrixType::specialOrthogonal) {
      randomSpecialOrthogonal(seed, matrix);
    } else if (matrixType == MatrixType::householder) {
      randomHouseholder(seed, matrix);
    } else if (matrixType == MatrixType::circulant) {
      randomOrthogonalCirculant(seed, matrix);
    } else if (matrixType == MatrixType::upperTriangular) {
      randomUpperTriangular(seed, matrix);
    } else if (matrixType == MatrixType::lowerTriangular) {
      randomLowerTriangular(seed, matrix);
    } else if (matrixType == MatrixType::schroeder) {
      randomSchroeder(seed, matrix);
    } else if (matrixType == MatrixType::absorbent) {
      randomAbsorbent(seed, matrix);
    } else if (matrixType == MatrixType::hadamard) {
      constructHadamardSylvester(matrix);
    } else if (matrixType == MatrixType::conference) {
      constructConference(matrix);
    } else { // matrixType == MatrixType::orthogonal, or default.
      randomOrthogonal(seed, matrix);
    }
  }

  void write(std::string name, MatrixType matrixType, unsigned seed = 0)
  {
    auto start = std::chrono::steady_clock::now();
    randomizeMatrix(matrixType, seed);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << std::format("{:24}, seed={:}: {} ms\n", name, seed, elapsed.count());

    std::string text{"["};
    for (size_t i = 0; i < length; ++i) {
      text += "[";
      for (size_t j = 0; j < length; ++j) text += std::format("{},", matrix[i][j]);
      text.pop_back();
      text += "],";
    }
    text.pop_back();
    text += "]\n";

    std::ofstream os{"json/" + name + ".json"};
    os << text;
  }
};

int main()
{
  using FBMat = FeedbackMatrix<float, 16>;
  FBMat matrix;
  unsigned seed = 123456;

  std::cout << "--- Warm up\n";
  matrix.write("absorbent", FBMat::MatrixType::absorbent, seed);
  matrix.write("absorbent", FBMat::MatrixType::absorbent, seed);

  std::cout << "\n--- Benchmark\n";
  matrix.write("orthogonal", FBMat::MatrixType::orthogonal, seed);
  matrix.write("specialOrthogonal", FBMat::MatrixType::specialOrthogonal, seed);
  matrix.write("householder", FBMat::MatrixType::householder, seed);
  matrix.write("circulant", FBMat::MatrixType::circulant, seed);
  matrix.write("upperTriangular", FBMat::MatrixType::upperTriangular, seed);
  matrix.write("lowerTriangular", FBMat::MatrixType::lowerTriangular, seed);
  matrix.write("schroeder", FBMat::MatrixType::schroeder, seed);
  matrix.write("absorbent", FBMat::MatrixType::absorbent, seed);
  matrix.write("hadamard", FBMat::MatrixType::hadamard, seed);
  matrix.write("conference", FBMat::MatrixType::conference, seed);

  return 0;
}
