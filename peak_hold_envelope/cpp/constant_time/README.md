Evaluation of the constant-time peak-hold described in following link.

- [Constant-time peak-hold  :  Blog  :  Signalsmith Audio](https://signalsmith-audio.co.uk/writing/2022/constant-time-peak-hold/)

A reference implementaion is available on the link below:

- https://github.com/Signalsmith-Audio/dsp/blob/5c7d1f3eb375b4862b682d310ccfbaafb6a4477e/envelopes.h#L362

There's a bug in the reference implementation (2025-02-25). In `PeakHold::resize()`, there's a following line:

```c++
while (bufferLength < maxLength) bufferLength *= 2;
```

It should be:

```c++
while (bufferLength <= maxLength) bufferLength *= 2;
```

The comparison is changed from `<` to `<=`. Otherwise, the peak-hold introduces spikes when `maxLength` is exactly 2^n.

Below is the license of `SignalsmithAudio::PeakHold`.

```
MIT License

Copyright (c) 2021 Geraint Luff / Signalsmith Audio Ltd.

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
