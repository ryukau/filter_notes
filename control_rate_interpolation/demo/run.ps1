#
# Note: Using libsndfile through vcpkg.
#

if (!(test-path build)) {
  cmake -S . -B build `
    -DCMAKE_TOOLCHAIN_FILE="/src/vcpkg/scripts/buildsystems/vcpkg.cmake" `
    -DCMAKE_BUILD_TYPE=Release
}

if (!(test-path snd)) {
  mkdir snd
}

cmake --build build --config release

./build/Release/naive.exe
./build/Release/linterp.exe
./build/Release/smoother.exe
./build/Release/emafilter.exe
./build/Release/rate_limiter.exe

Write-Output "--- bench_linterp"
./build/Release/bench_linterp.exe

Write-Output "--- bench_smoother"
./build/Release/bench_smoother.exe
