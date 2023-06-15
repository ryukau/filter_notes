To build, run following command. CMake and C++ 20 compiler are required.

```bash
cd biquad_filter_comparison/cpp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config release
./build/Release/benchmark
```

On Windows, above commands are already written in `build.ps1`.
