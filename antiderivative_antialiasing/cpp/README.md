# Antiderivative Antialiasing (ADAA)
A C++ implementation of ADAA for various nonlinearities.

To build on Windows, run `build.ps1`. libsndfile is required via vcpkg.

Following libaries were ported and placed on `lib`.

- [cephes](https://netlib.org/cephes/)
- [GitHub - Expander/polylogarithm: Implementation of polylogarithms in C/C++/Fortran](https://github.com/Expander/polylogarithm)

## License
### Cephes
Below is the readme of cephes obtained in 2025-02-17. It mentions the license.

```
   Some software in this archive may be from the book _Methods and
Programs for Mathematical Functions_ (Prentice-Hall or Simon & Schuster
International, 1989) or from the Cephes Mathematical Library, a
commercial product. In either event, it is copyrighted by the author.
What you see here may be used freely but it comes with no support or
guarantee.

   The two known misprints in the book are repaired here in the
source listings for the gamma function and the incomplete beta
integral.


   Stephen L. Moshier
   moshier@na-net.ornl.gov
```

### Polylogarithm
Below is the license of polylogarithm obtained in 2025-02-17.

```
MIT License

Copyright (c) 2021 Alexander Voigt

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
