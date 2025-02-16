if (!(test-path build)) {
  cmake -S . -B build `
    -DCMAKE_TOOLCHAIN_FILE="/src/vcpkg/scripts/buildsystems/vcpkg.cmake" `
    -DCMAKE_BUILD_TYPE=Release
}

if (!(test-path snd)) {
  mkdir snd
  python synth.py
}

if (!(test-path img)) {
  mkdir img
}

cmake --build build --config release
# ./build/Release/adaa.exe
