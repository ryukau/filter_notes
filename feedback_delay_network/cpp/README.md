## Usage
This instruction uses PowerShell. Following is required.

- Python 3
  - SciPy, Numpy, matplotlib
- cl.exe (MSVC C++ compiler)

Make directory `json` and `img`.

```powershell
mkdir json
mkdir img
```

Build `test.cpp` and run. C++20 is required.

```powershell
cl.exe /std:c++latest /EHsc /O2 /fp:fast test.cpp
.\test.exe
```

`test.exe` writes randomized matrices into `json` directory.

Run `test.py` to plot results.

```powershell
python test.py
```

`test.py` writes some plots of randomized matrices to `img`.

---

`randomunitary.py` is a Python 3 script I used for prototyping.
