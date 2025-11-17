import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import soundfile
import math
import json
import tqdm
from pathlib import Path
from collections import deque


def lagrange3Interp(y0: float, y1: float, y2: float, y3: float, t: float):
    u = 1 + t
    d0 = y0 - y1
    d1 = d0 - (y1 - y2)
    d2 = d1 - ((y1 - y2) - (y2 - y3))
    return y0 - u * (d0 + (1 - u) / 2 * (d1 + (2 - u) / 3 * d2))


class SosFilter:
    def __init__(self, nSection):
        self.x1 = np.zeros(nSection)
        self.x2 = np.zeros(nSection)
        self.y1 = np.zeros(nSection)
        self.y2 = np.zeros(nSection)

    def process(self, x0: float, sos) -> float:
        for i, co in enumerate(sos):
            y0 = (
                co[0] * x0
                + co[1] * self.x1[i]
                + co[2] * self.x2[i]
                - co[4] * self.y1[i]
                - co[5] * self.y2[i]
            )

            self.x2[i] = self.x1[i]
            self.x1[i] = x0
            self.y2[i] = self.y1[i]
            self.y1[i] = y0

            x0 = y0
        return x0


class DelayAaIir:
    iirOrder = 16

    def __init__(self, maxTimeSample: int):
        self.prevTime: float = 0
        self.holdingPitch: float = 1
        self.lowpass = SosFilter(self.iirOrder // 2)
        self.wptr: int = 0
        self.buf: np.ndarray = np.zeros((maxTimeSample + 1))

        self.pitch: deque = deque([0] * (maxTimeSample + 1))

        self.minTimeSample: int = 1
        self.maxTimeSample: int = maxTimeSample

    def process(self, input: float, timeInSample: float):
        size: int = len(self.buf)

        clamped: float = np.clip(timeInSample, self.minTimeSample, self.maxTimeSample)

        inputPitch = self.prevTime - clamped + 1
        self.pitch.append(inputPitch)
        self.prevTime = clamped
        carry = 1
        while True:
            self.pitch[0] -= carry

            carry = -self.pitch[0]
            if carry <= 0:
                if len(self.pitch) == 1:
                    self.holdingPitch = inputPitch
                break

            self.pitch.popleft()
            if len(self.pitch) > 0:
                self.holdingPitch = self.pitch[0]
            else:
                self.holdingPitch = 1
                break

        cutoff: float = 0.5 if self.holdingPitch <= 1 else 2.0 ** (-self.holdingPitch)
        cutoff = min(cutoff, 0.45)
        sos = signal.butter(self.iirOrder, cutoff, output="sos", fs=1)
        lp = self.lowpass.process(input, sos)

        # Write to buffer.
        timeInt: int = int(clamped)
        fraction: float = clamped - timeInt

        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = lp

        # Read from buffer.
        rptr0: int = self.wptr - timeInt
        rptr1: int = rptr0 - 1
        if rptr0 < 0:
            rptr0 += size
        if rptr1 < 0:
            rptr1 += size

        return self.buf[rptr0] + fraction * (self.buf[rptr1] - self.buf[rptr0])

    def debug(self, input: float, timeInSample: float):
        clamped: float = np.clip(timeInSample, self.minTimeSample, self.maxTimeSample)

        inputPitch = self.prevTime - clamped + 1
        self.pitch.append(inputPitch)
        self.prevTime = clamped
        carry = 1
        while True:
            self.pitch[0] -= carry

            carry = -self.pitch[0]
            if carry <= 0:
                if len(self.pitch) == 1:
                    self.holdingPitch = inputPitch
                break

            self.pitch.popleft()
            if len(self.pitch) > 0:
                self.holdingPitch = self.pitch[0]
            else:
                self.holdingPitch = 1
                break

        cutoff: float = 0.5 if self.holdingPitch <= 1 else 2.0 ** (-self.holdingPitch)
        cutoff = min(cutoff, 0.45)
        sos = signal.butter(self.iirOrder, cutoff, output="sos", fs=1)
        input = self.lowpass.process(input, sos)

        return [inputPitch, self.holdingPitch, cutoff, input]


class Delay:
    maxTap = 256
    maxCutoffOctave: int = int(np.round(np.log2(maxTap)))
    recursiveOscFractionMrgin = 1024 * np.finfo(np.float64).eps

    def __init__(self, maxTimeSample: int):
        self.prevTime: float = 0
        self.wptr: int = 0
        self.buf: np.ndarray = np.zeros(
            max(self.maxTap, maxTimeSample + self.maxTap // 2 + 1)
        )

        self.minTimeSample: int = self.maxTap // 2 - (self.maxTap + 1) % 2
        self.maxTimeSample: int = maxTimeSample

    def processLinear(self, input: float, timeInSample: float):
        size: int = len(self.buf)
        clamped: float = np.clip(timeInSample, 0, size - 2)
        timeInt: int = int(clamped)
        rFraction: float = clamped - timeInt

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        rptr0: int = self.wptr - timeInt
        rptr1: int = rptr0 - 1
        if rptr0 < 0:
            rptr0 += size
        if rptr1 < 0:
            rptr1 += size

        return self.buf[rptr0] + rFraction * (self.buf[rptr1] - self.buf[rptr0])

    def processCubic(self, input: float, timeInSample: float):
        size: int = len(self.buf)
        clamped: float = np.clip(timeInSample - 1, 1, size - 4)
        timeInt: int = int(clamped)
        rFraction: float = clamped - timeInt

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        rptr0: int = self.wptr - timeInt
        rptr1: int = rptr0 - 1
        rptr2: int = rptr0 - 2
        rptr3: int = rptr0 - 3
        if rptr0 < 0:
            rptr0 += size
        if rptr1 < 0:
            rptr1 += size
        if rptr2 < 0:
            rptr2 += size
        if rptr3 < 0:
            rptr3 += size
        return lagrange3Interp(
            self.buf[rptr0],
            self.buf[rptr1],
            self.buf[rptr2],
            self.buf[rptr3],
            rFraction,
        )

    def processSincRectNaive(self, input: float, timeInSample: float):
        """
        Lowpass cutoff is fixed to Nyquist frequency. It causes aliasing.
        """
        size: int = len(self.buf)

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        clamped: float = np.clip(timeInSample, self.minTimeSample, self.maxTimeSample)
        timeInt: int = int(clamped)
        fraction: float = clamped - timeInt
        fir = lowpassFirReversed(self.maxTap, 0.5, fraction)

        rptr: int = self.wptr - timeInt - self.minTimeSample - 1
        if rptr < 0:
            rptr += size

        sig = 0
        for idx in range(self.maxTap):
            sig += fir[idx] * self.buf[rptr]
            rptr += 1
            if rptr >= size:
                rptr = 0
        return sig

    def processSincRectFixedTap(self, input: float, timeInSample: float):
        """
        Lowpass cutoff is adaptively changed.

        Longer FIR taps sounds better, but short delay time can't be used. When sampling rate is 48000 Hz:

        - `self.maxTap = 16` sounds a bit rough.
        - `self.maxTap = 64` sounds okay.
        - `self.maxTap = 256` sounds clean.
        """
        size: int = len(self.buf)

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        clamped: float = np.clip(timeInSample, self.minTimeSample, self.maxTimeSample)
        timeInt: int = int(clamped)
        fraction: float = clamped - timeInt

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)
        fir = lowpassFirReversed(self.maxTap, cutoff, fraction)

        rptr: int = self.wptr - timeInt - self.minTimeSample - 1
        if rptr < 0:
            rptr += size

        sig = 0
        for idx in range(self.maxTap):
            sig += fir[idx] * self.buf[rptr]
            rptr += 1
            if rptr >= size:
                rptr = 0
        return sig

    def processSincRectFullExp2(self, input: float, timeInSample: float):
        """2^n taps only."""
        size: int = len(self.buf)

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        localTap = np.clip(1 << max(0, math.frexp(timeInSample)[1]), 2, self.maxTap)
        halfTap: int = localTap // 2
        clamped: float = np.clip(timeInSample, halfTap - 1, self.maxTimeSample)
        timeInt: int = int(clamped)
        fraction: float = clamped - timeInt

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)
        if timeInSample <= 0:
            return input * modifiedSinc(0, cutoff)
        fir = lowpassFirReversed(localTap, cutoff, fraction)

        rptr: int = self.wptr - timeInt - halfTap
        if rptr < 0:
            rptr += size

        sig = 0
        for idx in range(localTap):
            sig += fir[idx] * self.buf[rptr]
            rptr += 1
            if rptr >= size:
                rptr = 0
        return sig

    def processSincRectFullEven(self, input: float, timeInSample: float):
        """
        Even length FIR. Assuming that `self.maxTap` is even number.
        """
        size: int = len(self.buf)

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        localTap = np.clip(2 * int(timeInSample), 2, self.maxTap)
        halfTap: int = localTap // 2
        clamped: float = np.clip(timeInSample, halfTap - 1, self.maxTimeSample)
        timeInt: int = int(clamped)

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)
        if timeInSample <= 0:
            return input * modifiedSinc(0, cutoff)
        fir = lowpassFirReversed(localTap, cutoff, clamped - timeInt)

        rptr: int = self.wptr - timeInt - halfTap
        if rptr < 0:
            rptr += size

        sig = 0
        for idx in range(localTap):
            sig += fir[idx] * self.buf[rptr]
            rptr += 1
            if rptr >= size:
                rptr = 0
        return sig

    def processSincRectFullEven2(
        self, input: float, timeInSample: float, slewRate: float = 8
    ):
        """
        Slew rate limiting of time.
        """
        size: int = len(self.buf)

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        td = timeInSample - self.prevTime
        if td > slewRate:
            timeInSample = self.prevTime + slewRate
        elif td < -slewRate:
            timeInSample = self.prevTime - slewRate

        localTap = np.clip(2 * int(timeInSample), 2, self.maxTap)
        halfTap: int = localTap // 2
        clamped: float = np.clip(timeInSample, halfTap - 1, self.maxTimeSample)
        timeInt: int = int(clamped)

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)
        if timeInSample <= 0:
            return input * modifiedSinc(0, cutoff)
        fir = lowpassFirReversed(localTap, cutoff, clamped - timeInt)

        rptr: int = self.wptr - timeInt - halfTap
        if rptr < 0:
            rptr += size

        sig = 0
        for idx in range(localTap):
            sig += fir[idx] * self.buf[rptr]
            rptr += 1
            if rptr >= size:
                rptr = 0
        return sig

    def processSincRectFullEvenNormalized(self, input: float, timeInSample: float):
        """
        The amplitude of FIR filter is normalized to: `sum(fir) == 1`.
        """
        size: int = len(self.buf)

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        localTap = np.clip(2 * int(timeInSample), 2, self.maxTap)
        halfTap: int = localTap // 2
        clamped: float = np.clip(timeInSample, halfTap - 1, self.maxTimeSample)
        timeInt: int = int(clamped)
        fraction: float = clamped - timeInt

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)
        if timeInSample <= 0:
            return input * modifiedSinc(0, cutoff)
        fir = lowpassFirReversed(localTap, cutoff, fraction)
        if cutoff != 0:
            fir /= np.sum(fir)

        rptr: int = self.wptr - timeInt - halfTap
        if rptr < 0:
            rptr += size

        sig = 0
        for idx in range(localTap):
            sig += fir[idx] * self.buf[rptr]
            rptr += 1
            if rptr >= size:
                rptr = 0
        return sig

    def processSincRectFullEvenBiquadSine(self, input: float, timeInSample: float):
        """
        FIR filter is constructed using recursive sine algorithm (biquad oscillator).
        """
        size: int = len(self.buf)

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        localTap = np.clip(2 * int(timeInSample), 2, self.maxTap)
        halfTap: int = localTap // 2
        clamped: float = np.clip(timeInSample, halfTap - 1, self.maxTimeSample)
        timeInt: int = int(clamped)

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)
        if timeInSample <= 0:
            return input * modifiedSinc(0, cutoff)

        rptr: int = self.wptr - timeInt - halfTap
        if rptr < 0:
            rptr += size

        # Setup recursive sine oscillator.
        fraction = clamped - timeInt
        if fraction <= self.recursiveOscFractionMrgin:
            fraction = 0
        elif fraction >= 1 - self.recursiveOscFractionMrgin:
            fraction = 1
        mid = fraction - halfTap
        omega = 2 * np.pi * cutoff
        phi = mid * omega
        k = 2 * np.cos(omega)
        u1 = np.sin(phi - omega)
        u2 = np.sin(phi - 2 * omega)

        sig = 0
        for idx in range(localTap):
            u0 = k * u1 - u2
            u2 = u1
            u1 = u0

            x = idx + mid
            sinc = 2 * cutoff if x == 0 else u0 / (np.pi * x)

            sig += sinc * self.buf[rptr]
            rptr += 1
            if rptr >= size:
                rptr = 0
        return sig

    def processLanczosBiquadSine(self, input: float, timeInSample: float):
        """
        FIR filter is constructed using recursive sine algorithm (biquad oscillator).
        """
        size: int = len(self.buf)

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        localTap = np.clip(2 * int(timeInSample), 2, self.maxTap)
        halfTap: int = localTap // 2
        clamped: float = np.clip(timeInSample, halfTap - 1, self.maxTimeSample)
        timeInt: int = int(clamped)

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)
        if timeInSample <= 0:
            return input * modifiedSinc(0, cutoff)

        rptr: int = self.wptr - timeInt - halfTap
        if rptr < 0:
            rptr += size

        # Setup recursive sine oscillator.
        fraction = clamped - timeInt
        if fraction <= self.recursiveOscFractionMrgin:
            fraction = 0
        elif fraction >= 1 - self.recursiveOscFractionMrgin:
            fraction = 1
        mid = fraction - halfTap

        o1_omega = 2 * np.pi * cutoff
        o1_phi = mid * o1_omega
        o1_u1 = np.sin(o1_phi - o1_omega)
        o1_u2 = np.sin(o1_phi - 2 * o1_omega)
        o1_k = 2 * np.cos(o1_omega)

        a = max(1, np.sqrt(halfTap))
        o2_omega = 2 * np.pi * cutoff / a
        o2_phi = mid * o2_omega
        o2_u1 = np.sin(o2_phi - o2_omega)
        o2_u2 = np.sin(o2_phi - 2 * o2_omega)
        o2_k = 2 * np.cos(o2_omega)

        sig = 0
        A = 1 if cutoff == 0 else a / (2 * cutoff * np.pi * np.pi)
        for idx in range(localTap):
            o1_u0 = o1_k * o1_u1 - o1_u2
            o1_u2 = o1_u1
            o1_u1 = o1_u0

            o2_u0 = o2_k * o2_u1 - o2_u2
            o2_u2 = o2_u1
            o2_u1 = o2_u0

            x = idx + mid
            sinc = 2 * cutoff if x == 0 else A * o1_u0 * o2_u0 / (x * x)

            sig += sinc * self.buf[rptr]
            rptr += 1
            if rptr >= size:
                rptr = 0
        return sig

    def processSincRectFullEvenStableQuadSine(self, input: float, timeInSample: float):
        """
        FIR filter is constructed using recursive sine algorithm (stable quadrature oscillator).
        """
        size: int = len(self.buf)

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        localTap = np.clip(2 * int(timeInSample), 2, self.maxTap)
        halfTap: int = localTap // 2
        clamped: float = np.clip(timeInSample, halfTap - 1, self.maxTimeSample)
        timeInt: int = int(clamped)

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)
        if timeInSample <= 0:
            return input * modifiedSinc(0, cutoff)

        rptr: int = self.wptr - timeInt - halfTap
        if rptr < 0:
            rptr += size

        # Setup recursive sine oscillator.
        fraction = clamped - timeInt
        if fraction <= self.recursiveOscFractionMrgin:
            fraction = 0
        elif fraction >= 1 - self.recursiveOscFractionMrgin:
            fraction = 1
        mid = fraction - halfTap
        omega = 2 * np.pi * cutoff
        phi = mid * omega
        k1 = np.tan(omega / 2)
        k2 = np.sin(omega)
        u = np.cos(phi - omega)
        v = np.sin(phi - omega)

        sig = 0
        # gain = 0
        for idx in range(localTap):
            w = u - k1 * v
            v += k2 * w
            u = w - k1 * v

            x = idx + mid
            sinc = 2 * cutoff if x == 0 else v / (np.pi * x)

            # gain += sinc
            sig += sinc * self.buf[rptr]
            rptr += 1
            if rptr >= size:
                rptr = 0
        # return sig if gain == 0 else sig / gain
        return sig

    def processSincRectFullOdd(self, input: float, timeInSample: float):
        """
        Odd length FIR. Assuming that `self.maxTap` is even number.
        """
        size: int = len(self.buf)

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        if timeInSample < 1:
            return self.processSincRectFullEven(input, timeInSample)
        localTap = np.clip(2 * int(timeInSample) + 1, 3, self.maxTap + 1)
        halfTap: int = localTap // 2
        clamped: float = np.clip(timeInSample, halfTap, self.maxTimeSample)
        timeInt: int = int(clamped)
        fraction: float = clamped - timeInt

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)
        fir = lowpassFirReversed(localTap, cutoff, fraction)

        rptr: int = self.wptr - timeInt - halfTap - 1
        if rptr < 0:
            rptr += size

        sig = 0
        for idx in range(localTap):
            sig += fir[idx] * self.buf[rptr]
            rptr += 1
            if rptr >= size:
                rptr = 0
        return sig

    def processSincRectFull(self, input: float, timeInSample: float):
        size: int = len(self.buf)

        # Write to buffer.
        self.wptr += 1
        if self.wptr >= size:
            self.wptr = 0
        self.buf[self.wptr] = input

        # Read from buffer.
        localTap = np.clip(int(2 * timeInSample), 2, self.maxTap)
        isEven: int = (localTap + 1) & 1
        halfTap: int = localTap // 2
        clamped: float = np.clip(timeInSample, halfTap - isEven, self.maxTimeSample)
        timeInt: int = int(clamped)
        fraction: float = clamped - timeInt

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)
        fir = lowpassFirReversed(localTap, cutoff, fraction)

        rptr: int = self.wptr - timeInt - halfTap - (isEven ^ 1)
        if rptr < 0:
            rptr += size

        sig = 0
        for idx in range(localTap):
            sig += fir[idx] * self.buf[rptr]
            rptr += 1
            if rptr >= size:
                rptr = 0
        return sig

    def debugCutoff(self, timeInSample: float):
        localTap = np.clip(2 * int(timeInSample), 2, self.maxTap)
        halfTap: int = localTap // 2
        clamped: float = np.clip(timeInSample, halfTap - 1, self.maxTimeSample)

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)

        timeInt: int = int(clamped)
        fraction = clamped - timeInt
        mid = fraction - halfTap

        return [int(localTap), float(cutoff), float(fraction)]

    def debugDiffFir(self, timeInSample: float):
        localTap = np.clip(2 * int(timeInSample), 2, self.maxTap)
        halfTap: int = localTap // 2
        clamped: float = np.clip(timeInSample, halfTap - 1, self.maxTimeSample)
        timeInt: int = int(clamped)

        timeDiff: float = abs(self.prevTime - clamped + 1)
        self.prevTime = clamped
        cutoff: float = 0.5 if timeDiff <= 1 else 2.0 ** (-timeDiff)
        if timeInSample <= 0:
            return 2 * cutoff
        fir = lowpassFirReversed(localTap, cutoff, clamped - timeInt)

        # Setup recursive sine oscillator.
        mid = clamped - timeInt - halfTap
        omega = 2 * np.pi * cutoff
        phi = mid * omega
        k1 = np.tan(omega / 2)
        k2 = np.sin(omega)
        u = np.cos(phi - omega)
        v = np.sin(phi - omega)

        diff = 0
        for idx in range(localTap):
            w = u - k1 * v
            v += k2 * w
            u = w - k1 * v

            x = idx + mid
            sinc = 2 * cutoff if x == 0 else v / (np.pi * x)

            diff += abs(sinc - fir[idx])
        return diff


def plotResponse(firList, cutoffList, nameList, worN=8192, fs=48000, title=""):
    fig, ax = plt.subplots(4, 1)
    if len(title) >= 1:
        fig.suptitle(title)

    # worN = np.hstack([[0], np.geomspace(0.1, fs / 2)])

    gdMedians = []
    cmap = plt.get_cmap("plasma")
    for idx, fir in enumerate(firList):
        freq, resp = signal.freqz(fir, 1, worN=worN, fs=fs)
        freq, delay = signal.group_delay((fir, 1), w=worN, fs=fs)
        gain = 20 * np.log10(np.abs(resp))
        phase = np.unwrap(np.angle(resp))

        cut = max(1, int(len(delay) * cutoffList[idx] / fs))
        # print(f"delay: {np.mean(delay[:cut])}")
        gdMedians.append(np.median(delay[:cut]))

        color = cmap(idx / len(firList))
        ax[0].axvline(cutoffList[idx], alpha=0.33, color="black", ls="--")
        ax[0].plot(freq, gain, alpha=0.66, lw=1, color=color)
        ax[1].plot(freq, phase, alpha=0.66, lw=1, color=color)
        ax[2].plot(freq, delay, alpha=0.66, lw=1, color=color)
        ax[3].plot(fir, alpha=0.66, lw=1, color=color, label=f"{nameList[idx]}")

    ax[0].set_ylabel("Gain [dB]")
    ax[0].set_ylim((-40, 6))
    ax[1].set_ylabel("Phase [rad/sample]")
    ax[2].set_ylabel("Delay [sample]")
    gdMid = np.median(gdMedians)
    gdRange = np.max(
        [
            1.5 * (gdMid - np.min(gdMedians)),
            1.5 * (np.max(gdMedians) - gdMid),
            1,
        ]
    )
    ax[2].set_ylim([gdMid - gdRange, gdMid + gdRange])
    ax[3].set_ylabel("FIR Amplitude")
    ax[3].legend(ncol=2)

    for axis in ax[:3]:
        axis.set_xlim([100, 25000])
        axis.set_xscale("log")
        axis.axvline(fs / 2, color="black", ls="--")
        # axis.axvline(cutoffHz, color="black", ls="--")

    for axis in ax:
        axis.grid(color="#f0f0f0", which="both")
        # axis.legend(ncol=2)
        # axis.set_xscale("log")
    fig.set_size_inches((8, 8))
    fig.tight_layout()
    plt.show()


def modifiedSinc(x, cutoff):
    if x == 0:
        return 2 * cutoff
    return np.sin(np.pi * 2 * cutoff * x) / (np.pi * x)


def lowpassFir(length: int, cutoff: float, fractionSample: float):
    mid = length // 2 if length % 2 == 1 else length // 2 - 0.5

    fir = np.zeros(length)
    for i in range(length):
        x = i - mid
        fir[i] = modifiedSinc(x - fractionSample, cutoff)
    return fir


def lowpassFirReversed(length: int, cutoff: float, fractionSample: float):
    mid = length // 2 + length % 2
    mid -= fractionSample

    fir = np.zeros(length)
    for i in range(length):
        x = i - mid
        fir[i] = modifiedSinc(x, cutoff)
    return fir


def lowpassFirStableQuad(length: int, cutoff: float, fractionSample: float):
    mid = length // 2 + 1 if length % 2 == 1 else length // 2

    margin = 256 * np.finfo(np.float64).eps
    if fractionSample <= margin:
        fractionSample = 0
    elif fractionSample >= 1 - margin:
        fractionSample = 1
    mid -= fractionSample

    omega = 2 * np.pi * cutoff
    phi = -mid * omega
    k1 = np.tan(omega / 2)
    k2 = np.sin(omega)
    u = np.cos(phi - omega)
    v = np.sin(phi - omega)

    fir = np.zeros(length)
    for i in range(length):
        w = u - k1 * v
        v += k2 * w
        u = w - k1 * v

        x = i - mid
        fir[i] = 2 * cutoff if x == 0 else v / (np.pi * x)

        # print(f"{i:8}: {v:+20.16f} {np.sin(np.pi * 2 * cutoff * x):+20.16f}")
    # print("quad", fir.tolist())
    return fir


def lowpassFirReinsch(length: int, cutoff: float, fractionSample: float):
    mid = length // 2 + 1 if length % 2 == 1 else length // 2

    margin = 256 * np.finfo(np.float64).eps
    if fractionSample <= margin:
        fractionSample = 0
    elif fractionSample >= 1 - margin:
        fractionSample = 1
    mid -= fractionSample

    omega = 2 * np.pi * cutoff
    phi = -mid * omega
    A = 2 * np.sin(omega / 2)
    u = np.sin(phi - omega)
    v = A * np.cos(phi - omega / 2)
    k = A * A

    fir = np.zeros(length)
    for i in range(length):
        u = u + v
        v = v - k * u

        x = i - mid
        if x == 0:
            fir[i] = 2 * cutoff
        else:
            fir[i] = u / (np.pi * x)
    return fir


def lowpassFirBiquad(length: int, cutoff: float, fractionSample: float):
    mid = length // 2 + 1 if length % 2 == 1 else length // 2

    margin = 256 * np.finfo(np.float64).eps
    if fractionSample <= margin:
        fractionSample = 0
    elif fractionSample >= 1 - margin:
        fractionSample = 1
    mid -= fractionSample

    omega = 2 * np.pi * cutoff
    phi = -mid * omega
    k = 2 * np.cos(omega)
    u1 = np.sin(phi - omega)
    u2 = np.sin(phi - 2 * omega)

    fir = np.zeros(length)
    for i in range(length):
        u0 = k * u1 - u2
        u2 = u1
        u1 = u0

        x = i - mid
        if x == 0:
            fir[i] = 2 * cutoff
        else:
            fir[i] = u0 / (np.pi * x)
    return fir


def lowpassLanczosReversed(
    length: int, cutoff: float, fractionSample: float, lanczos_a: float = 3
):
    mid = length // 2 + 1 if length % 2 == 1 else length // 2
    mid -= fractionSample

    fir = np.zeros(length)
    for i in range(length):
        x = i - mid
        lanczos = modifiedSinc(x / lanczos_a, cutoff) / cutoff / 2
        fir[i] = modifiedSinc(x, cutoff) * lanczos
    return fir


def lowpassWelchReversed(length: int, cutoff: float, fractionSample: float):
    mid = length // 2 + 1 if length % 2 == 1 else length // 2
    mid -= fractionSample

    fir = np.zeros(length)
    M = max(1, (length - 1) / 2)
    for i in range(length):
        A = (i - M) / M
        window = 1 - A * A
        fir[i] = modifiedSinc(i - mid, cutoff) * window
    return fir


def lowpassWindowedReversed(length: int, cutoff: float, fractionSample: float):
    mid = length // 2 + 1 if length % 2 == 1 else length // 2
    mid -= fractionSample

    fir = np.zeros(length)
    for i in range(length):
        fir[i] = modifiedSinc(i - mid, cutoff)
    return fir
    # return fir * signal.get_window("blackman", length)
    # return fir * signal.get_window("blackmanharris", length)
    # return fir * signal.get_window(("chebwin", 120), length)
    # return fir * signal.get_window("cosine", length)
    # return fir * signal.get_window(("dpss", 4), length)
    # return fir * signal.get_window(("exponential", None, 16), length)
    # return fir * signal.get_window(("gaussian", length / 9), length)
    # return fir * signal.get_window("hann", length)
    # return fir * signal.get_window("hamming", length)
    # return fir * signal.get_window(("kaiser", 4 * np.pi), length)
    # return fir * signal.get_window("nuttall", length)
    # return fir * signal.get_window("triang", length)


def lowpassLanczosBiquad(
    length: int, cutoff: float, fractionSample: float, lanczos_width: int = 8
):
    mid = length // 2 + 1 if length % 2 == 1 else length // 2

    margin = 256 * np.finfo(np.float64).eps
    if fractionSample <= margin:
        fractionSample = 0
    elif fractionSample >= 1 - margin:
        fractionSample = 1
    mid -= fractionSample

    o1_omega = 2 * np.pi * cutoff
    o1_phi = -mid * o1_omega
    o1_u1 = np.sin(o1_phi - o1_omega)
    o1_u2 = np.sin(o1_phi - 2 * o1_omega)
    o1_k = 2 * np.cos(o1_omega)

    o2_omega = 2 * np.pi * cutoff / float(lanczos_width)
    o2_phi = -mid * o2_omega
    o2_u1 = np.sin(o2_phi - o2_omega)
    o2_u2 = np.sin(o2_phi - 2 * o2_omega)
    o2_k = 2 * np.cos(o2_omega)

    fir = np.zeros(length)
    A = lanczos_width / (2 * cutoff * np.pi * np.pi)
    for i in range(length):
        o1_u0 = o1_k * o1_u1 - o1_u2
        o1_u2 = o1_u1
        o1_u1 = o1_u0

        o2_u0 = o2_k * o2_u1 - o2_u2
        o2_u2 = o2_u1
        o2_u1 = o2_u0

        x = i - mid
        if x == 0:
            fir[i] = 2 * cutoff
        else:
            fir[i] = A * o1_u0 * o2_u0 / (x * x)
    return fir


def testFir():
    sampleRate = 48000
    tap = 256
    fraction = 0.0

    nOctave = int(np.round(np.log2(tap)))

    fir = []
    cutoff = []
    name = []
    for idx in range(2, nOctave + 1):
        cutoff.append(48000 / 2**idx)
        name.append(str(cutoff[-1]))
        fir.append(lowpassWindowedReversed(tap, cutoff[-1] / sampleRate, fraction))
    plotResponse(fir, cutoff, name, fs=sampleRate)


def compareFir():
    # Labels originally represented times in the signal that these bad cases were found.
    # They don't have any particular meaning in this function.
    #
    # {label: [localTap, cutoff, fraction], ...}
    badCases = {
        2000: [126, 0.5, 0.9999999999999929],
        10000: [126, 0.4862462530029181, 0.9999999999999929],
        12000: [32, 0.4959896743541015, 7.105427357601002e-15],
        14000: [8, 0.49912897434665304, 2.6645352591003757e-15],
        22000: [6, 0.5, 0.9999999999999951],
        24000: [30, 0.5, 0.9999999999999893],
        34000: [126, 0.4862462530029589, 0.9999999999999574],
        36000: [32, 0.49598967435411495, 1.7763568394002505e-14],
        38000: [8, 0.4991289743466527, 8.881784197001252e-15],
        46000: [8, 0.5, 8.881784197001252e-16],
    }

    for label, param in badCases.items():
        sampleRate = 48000
        # tap = 126
        # cutoff = 0.5
        # fraction = 0.9999999999999929
        # param = [tap, cutoff, fraction]
        # param = badCases[46000]

        fir = []
        name = ["ref.", "quad."]
        fir.append(lowpassFirReversed(*param))
        fir.append(lowpassFirBiquad(*param))
        print(f"diff on {label}: {np.sum(np.abs(fir[0]-fir[1]))}")
        plotResponse(
            fir, [param[1] * sampleRate] * 2, name, fs=sampleRate, title=str(label)
        )
    plt.show()


def generateSawtooth(
    sampleRate: float, frequencyHz: float, durationSecond: float
) -> np.ndarray:
    """Integer sample period only. It makes easier to check the aliasing on spectrogram."""
    period: int = int(sampleRate / frequencyHz)
    durationSample = int(sampleRate * durationSecond)
    if period == 0:
        return np.zeros(durationSample)
    sig = np.arange(durationSample, dtype=np.int64) % period
    return sig.astype(np.float64) / float(period) * 2 - 1


def plotSpectrogram(sampleRateHz: float, sig: np.ndarray):
    frameSize = 1024
    win = signal.get_window("hann", frameSize)
    sft = signal.ShortTimeFFT(
        win, len(win) // 128, fs=sampleRateHz, scale_to="magnitude"
    )
    mag = sft.stft(sig)

    plt.figure(figsize=(10, 5))
    im1 = plt.imshow(
        np.log(abs(mag) + 1e-7),
        origin="lower",
        aspect="auto",
        extent=sft.extent(len(sig)),
        cmap="magma",
    )
    plt.xlabel("Time [s]")
    plt.ylabel("Frequency [Hz]")
    # plt.ylim([10, 20000])
    # plt.yscale("log")
    # plt.colorbar(im1, label="Magnitude $|S_x(t, f)|$")
    plt.tight_layout()


def showSpectrogram(sampleRateHz: float, sig: np.ndarray):
    plotSpectrogram(sampleRateHz, sig)
    plt.show()
    plt.close()


def saveSpectrogram(path: Path, sampleRateHz: float, sig: np.ndarray):
    plotSpectrogram(sampleRateHz, sig)
    plt.savefig(path)
    plt.close()


def testGroupDelay():
    sampleRate = 48000
    worN = 2048
    print(f"{"cutoff [Hz]":24} | {"median delay [sample]":24}")
    for idx in range(2, 9):
        cutoff = 48000 / 2**idx
        fir = lowpassFirReversed(16, cutoff / sampleRate, -0.5)

        _, gd = signal.group_delay((fir, 1), w=worN, fs=sampleRate)

        cut = int(len(gd) * cutoff / sampleRate)
        delay = np.mean(gd[:cut])
        print(f"{cutoff:24} | {delay:24}")


def testDelayTime():
    impulse = np.zeros(32)
    impulse[0] = 1

    maxDelayTimeSample = 1024
    delayTimeSample = 2.5
    for method in [m for m in dir(Delay) if "process" in m]:
        delay = Delay(maxDelayTimeSample)
        func = getattr(delay, method)
        delayed = np.zeros_like(impulse)
        for i, v in enumerate(impulse):
            delayed[i] = func(v, delayTimeSample)

        plt.figure()
        plt.title(method)
        plt.plot(impulse, alpha=0.75, lw=1, color="red", label="impulse")
        plt.plot(delayed, alpha=0.75, lw=1, color="blue", label="delayed")
        plt.xlim([-1, 2 + max(2, delayTimeSample)])
        plt.legend()
        plt.grid()
        plt.savefig(f"img/testDelayTime_{method}.png")
    # plt.show()


def testMaxDelayTime():
    impulse = np.zeros(32)
    impulse[0] = 1

    data = {}
    maxDelayTimeSample: float = 8
    delayTimeSample: float = 8
    for method in [m for m in dir(Delay) if "process" in m]:
        delay = Delay(maxDelayTimeSample)
        func = getattr(delay, method)
        delayed = np.zeros_like(impulse)
        for i, v in enumerate(impulse):
            delayed[i] = func(v, delayTimeSample)
        data[method] = list(delayed)

    with open("textMaxDelayTime.json", "w", encoding="utf-8") as fp:
        json.dump(data, fp)


def testDelay(inputPath: Path | None = None):
    def loadSound(path: Path | None):
        if path is None or not path.exists():
            sampleRate = 48000
            return (sampleRate, generateSawtooth(sampleRate, 4000, 1))
        data, fs = soundfile.read(path, always_2d=True)
        return (fs, data.T[0])

    sampleRate, source = loadSound(inputPath)
    saveSpectrogram(Path("img/testDelay_source.png"), sampleRate, source)

    # delayTime = np.linspace(8001, 1, len(source))
    # delayTime = 8001 - np.geomspace(1, 8000, len(source))
    # delayTime = np.geomspace(8000, 1, len(source))
    # delayTime = np.full_like(source, 8000)
    # delayTime = np.hstack(
    #     [
    #         np.linspace(0, 512, len(source) // 4),
    #         np.linspace(512, 0, len(source) // 4),
    #         np.linspace(0, 512, len(source) // 4),
    #         np.linspace(512, 0, len(source) // 4),
    #     ]
    # )
    # delayTime = 200 * 2 ** (8 * np.sin(4 * np.pi * np.arange(len(source)) / sampleRate))
    delayTime = 16 * 2 ** (4 * np.sin(4 * np.pi * np.arange(len(source)) / sampleRate))

    outputGain = 0.25
    maxDelayTimeSample = 65536
    # for method in [m for m in dir(Delay) if "process" in m]:
    for method in tqdm.tqdm(
        [
            "processLinear",
            "processSincRectFullEven",
            "processLanczosBiquadSine",
        ]
    ):
        delay = Delay(maxDelayTimeSample)
        func = getattr(delay, method)
        delayed = np.zeros_like(source)
        for i, v in enumerate(source):
            delayed[i] = func(v, delayTime[i])

        saveSpectrogram(Path(f"img/testDelay_{method}.png"), sampleRate, delayed)

        soundfile.write("snd/source.wav", outputGain * source, sampleRate, "FLOAT")
        soundfile.write(f"snd/{method}.wav", outputGain * delayed, sampleRate, "FLOAT")


def testDelaySlewRate():
    sampleRate = 48000

    sineFrequencyHz = 1111
    durationSecond = 1

    freq: float = sineFrequencyHz / sampleRate
    durationSample = int(sampleRate * durationSecond)
    source = np.sin(np.linspace(0, 2 * np.pi * freq * durationSample, durationSample))

    # source = generateSawtooth(sampleRate, 4000, 1)

    soundfile.write(
        "snd/testDelaySlewRate_source.wav", 0.25 * source, sampleRate, "FLOAT"
    )

    for slewRate in tqdm.tqdm(range(1, 9)):
        mod = 3  # Delay time modulation.
        minDelayTimeSample = 512
        maxDelayTimeSample = minDelayTimeSample + 2**mod

        delay = Delay(maxDelayTimeSample)
        baseTime = minDelayTimeSample

        y1 = 0
        sig = np.zeros_like(source)
        for i, v in enumerate(source):
            y1 = delay.processSincRectFullEven2(v, baseTime * 2 ** (mod * y1), slewRate)
            sig[i] = y1

        saveSpectrogram(Path(f"img/testDelaySlewRate_{slewRate}.png"), sampleRate, sig)

        soundfile.write(
            f"snd/testDelaySlewRate_{slewRate}.wav", 0.25 * sig, sampleRate, "FLOAT"
        )


def compareDelay():
    methodTarget = "processSincRectFullEven"
    methodActual = "processSincRectFullEven2"
    # methodActual = "processSincRectFullEvenBiquadSine"
    # methodActual = "processSincRectFullEvenStableQuadSine"

    sampleRate = 48000
    maxDelayTime = sampleRate
    source = generateSawtooth(sampleRate, 4000, 1)

    # delayTime = 256 * 2 ** (8 * np.sin(4 * np.pi * np.arange(len(source)) / sampleRate))
    delayTime = 16 * 2 ** (4 * np.sin(4 * np.pi * np.arange(len(source)) / sampleRate))

    delayTarget = Delay(maxDelayTime)
    funcTarget = getattr(delayTarget, methodTarget)
    sigTarget = np.zeros_like(source)

    delayActual = Delay(maxDelayTime)
    funcActual = getattr(delayActual, methodActual)
    sigActual = np.zeros_like(source)

    diffTracker = Delay(maxDelayTime)
    diff = np.zeros_like(source)
    for i, v in enumerate(source):
        sigTarget[i] = funcTarget(v, delayTime[i])
        sigActual[i] = funcActual(v, delayTime[i])
        diff[i] = diffTracker.debugDiffFir(delayTime[i])

    frameSize = 1024
    win = signal.get_window("hann", frameSize)
    sft = signal.ShortTimeFFT(win, len(win) // 128, fs=sampleRate, scale_to="magnitude")

    extent = list(sft.extent(len(sigTarget)))
    extent[0] *= sampleRate
    extent[1] *= sampleRate

    fig, ax = plt.subplots(3, 1, height_ratios=[3, 3, 1])
    fig.set_size_inches((12, 8))

    mag0 = np.log(abs(sft.stft(sigTarget)) + 1e-7)
    ax[0].imshow(mag0, origin="lower", aspect="auto", extent=extent, cmap="magma")
    ax[0].set_ylabel("Frequency [Hz]")

    mag1 = np.log(abs(sft.stft(sigActual)) + 1e-7)
    ax[1].imshow(mag1, origin="lower", aspect="auto", extent=extent, cmap="magma")
    ax[1].set_ylabel("Frequency [Hz]")

    ax[2].plot(diff)
    ax[2].set_xlabel("Time [sample]")
    ax[2].set_ylabel("diff")
    ax[2].set_xlim([extent[0], extent[1]])
    # ax[2].set_ylim([0.45, 0.51])

    for row in ax:
        row.grid(alpha=0.5)
    fig.tight_layout()
    plt.show()


def printFirGlitch():
    sampleRate = 48000
    maxDelayTime = sampleRate
    duration: int = sampleRate

    delayTime = 16 * 2 ** (4 * np.sin(4 * np.pi * np.arange(duration) / sampleRate))

    cutoffTracker = Delay(maxDelayTime)
    cutoff = np.zeros_like(delayTime)
    debugTracker = Delay(maxDelayTime)
    firParam = []
    for i in range(sampleRate):
        cutoff[i] = cutoffTracker.debugDiffFir(delayTime[i])
        firParam.append(debugTracker.debugCutoff(delayTime[i]))

    indices = np.argwhere(cutoff > 0.001)
    for index in indices:
        idx = index[0]
        print(f"{idx}: {firParam[idx]},")


def testSosFilter():
    sampleRate = 48000
    order = 8
    cutoff = 1000

    impulse = np.zeros(sampleRate)
    impulse[0] = 1

    sos = signal.butter(order, cutoff, output="sos", fs=sampleRate)
    target = signal.sosfilt(sos, impulse)

    sf = SosFilter(order // 2)
    actual = np.zeros_like(impulse)
    for i in range(len(impulse)):
        actual[i] = sf.process(impulse[i], sos)

    plotResponse(
        [target, actual],
        [cutoff] * 2,
        ["target", "actual"],
        fs=sampleRate,
        title="SosFilter Test",
    )


def testDelayAaIir(inputPath: Path | None = None):
    def loadSound(path: Path | None):
        if path is None or not path.exists():
            sampleRate = 48000
            return (sampleRate, generateSawtooth(sampleRate, 4000, 1))
        data, fs = soundfile.read(path, always_2d=True)
        return (fs, data.T[0])

    sampleRate, source = loadSound(inputPath)
    saveSpectrogram(Path("img/testDelayAaIir_source.png"), sampleRate, source)

    delayTime = np.linspace(8001, 1, len(source))
    # delayTime = 1001 - np.geomspace(1, 1000, len(source))
    # delayTime = np.geomspace(8000, 1, len(source))
    # delayTime = np.full_like(source, 8000)
    # delayTime = np.hstack(
    #     [
    #         np.linspace(0, 512, len(source) // 4),
    #         np.linspace(512, 0, len(source) // 4),
    #         np.linspace(0, 512, len(source) // 4),
    #         np.linspace(512, 0, len(source) // 4),
    #     ]
    # )
    # delayTime = 200 * 2 ** (8 * np.sin(4 * np.pi * np.arange(len(source)) / sampleRate))
    # delayTime = 16 * 2 ** (4 * np.sin(4 * np.pi * np.arange(len(source)) / sampleRate))

    # Test if antialising is effective.
    assert not np.any(np.diff(delayTime) < -1)

    outputGain = 0.25
    maxDelayTimeSample = 65536
    delay = DelayAaIir(maxDelayTimeSample)
    delayed = np.zeros_like(source)
    for i, v in enumerate(source):
        delayed[i] = delay.process(v, delayTime[i])

    saveSpectrogram(Path(f"img/testDelayAaIir_output.png"), sampleRate, delayed)

    soundfile.write("snd/source.wav", outputGain * source, sampleRate, "FLOAT")
    soundfile.write(
        f"snd/testDelayAaIir.wav", outputGain * delayed, sampleRate, "FLOAT"
    )


def debugDelayAaIir():
    sampleRate = 48000
    source = generateSawtooth(int(0.25 * sampleRate), 4000, 1)

    delayTime = np.linspace(8001, 1, len(source))
    # delayTime = np.hstack(
    #     [
    #         np.linspace(0, 512, len(source) // 4),
    #         np.linspace(512, 0, len(source) // 4),
    #         np.linspace(0, 512, len(source) // 4),
    #         np.linspace(512, 0, len(source) // 4),
    #     ]
    # )
    # delayTime = 128 * 2 ** (4 * np.sin(4 * np.pi * np.arange(len(source)) / sampleRate))

    # Test if antialising is effective.
    assert not np.any(np.diff(delayTime) < -1)

    # diff = np.diff(delayTime)
    # plt.plot(diff)
    # plt.grid()
    # plt.show()
    # exit()

    maxDelayTimeSample = 65536
    delay = DelayAaIir(maxDelayTimeSample)
    inputPitch = np.zeros_like(source)
    holdingPitch = np.zeros_like(source)
    cutoff = np.zeros_like(source)
    lowpassed = np.zeros_like(source)
    for i, v in enumerate(source):
        data = delay.debug(v, delayTime[i])
        inputPitch[i] = data[0]
        holdingPitch[i] = data[1]
        cutoff[i] = data[2]
        lowpassed[i] = data[3]

    fig, ax = plt.subplots(4, 1)
    fig.set_size_inches((8, 8))

    ax[0].plot(delayTime)
    ax[0].set_ylabel("Delay Time [sample]")

    ax[1].plot(inputPitch)
    ax[1].set_ylabel("Input Pitch")
    ax[1].set_ylim([0, 2])

    ax[2].plot(holdingPitch)
    ax[2].set_ylabel("Holding Pitch")

    ax[3].plot(cutoff)
    ax[3].set_ylabel("Cutoff [rad/2π]")

    for axis in ax:
        axis.grid()

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    Path("img").mkdir(parents=True, exist_ok=True)
    Path("snd").mkdir(parents=True, exist_ok=True)

    testFir()
    # compareFir()
    # testGroupDelay()
    # testDelayTime()
    # testMaxDelayTime()
    # testDelay(Path("snd/yey.wav"))
    # testDelay()
    # testDelaySlewRate()
    # compareDelay()
    # printFirGlitch()
    # testSosFilter()
    # testDelayAaIir()
    # debugDelayAaIir()
