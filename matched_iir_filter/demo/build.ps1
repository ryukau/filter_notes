if (!(test-path build)) {
  cmake -S . -B build `
    -DCMAKE_BUILD_TYPE=Release
}

cmake --build build --config release
./build/Release/test.exe
python matchediir.py
