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
./build/Release/test.exe
python plot.py
python test.py
