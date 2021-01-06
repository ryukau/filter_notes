if (!(test-path build)) {
  cmake -S . -B build `
    -DCMAKE_TOOLCHAIN_FILE="/src/vcpkg/scripts/buildsystems/vcpkg.cmake" `
    -DCMAKE_BUILD_TYPE=Release
}

# python synth.py

cmake --build build --config release
./build/Release/decimation.exe
python plot.py
