To build, run following command. CMake and C++ 20 compiler are required.

```bash
cd matched_iir_filter/demo
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config release
./build/Release/benchmark

python3 matchediir.py
```

On Windows, above commands are already written in `build.ps1`.

## `matchediir.py`
`testResponses()` can be used to plot filter responses without dealing with C++. To use:

- Uncomment `testResponses()` in main.
- Comment out `compareResponsesToCpp()` in main.

Optionally, comment out particular plot in `testResponses()`.
