import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy.special as special
import soundfile
import mpmath
from pathlib import Path


def generateSine(durationSamples: int, frequencyNormalized: float) -> np.ndarray:
    phase = 0
    sig = np.empty(durationSamples)
    for i in range(len(sig)):
        sig[i] = np.sin(2 * np.pi * phase)
        phase += frequencyNormalized
        phase -= np.floor(phase)
    return sig


def adaa1(x, f0, f1):
    tolerance = 1 / 2**24
    y = np.zeros_like(x)
    x1 = 0
    s1 = 0
    for n in range(len(x)):
        x0 = x[n]
        s0 = f1(x0)
        if (x1 == 0 and s1 == 0) or np.abs(x0 - x1) < tolerance:  # fallback
            y[n] = f0((x0 + x1) / 2)
        else:
            y[n] = (s0 - s1) / (x0 - x1)
        s1 = s0
        x1 = x0
    return y


def adaa2(x, f0, f1, f2):
    tolerance = 1 / 2**24
    y = np.zeros_like(x)

    x1 = 0
    x2 = 0
    s1 = 0
    for n in range(len(x)):
        x0 = x[n]

        f2_x1 = f2(x1)
        s0 = (
            f1((x0 + x1) / 2)
            if np.abs(x0 - x1) < tolerance
            else (f2(x0) - f2_x1) / (x0 - x1)
        )

        if x1 == 0 and x2 == 0:
            y[n] = f0((x0 + 2 * x1 + x2) / 4)
        elif np.abs(x0 - x2) < tolerance:
            x_bar = (x0 + x2) / 2
            delta = x_bar - x1
            if np.abs(delta) < tolerance:
                y[n] = f0((x_bar + x1) / 2)
            else:
                y[n] = (2 / delta) * (f1(x_bar) + (f2_x1 - f2(x_bar)) / delta)
        else:
            y[n] = 2 * (s0 - s1) / (x0 - x2)
        s1 = s0
        x2 = x1
        x1 = x0
    return y


def hardclipJ0(x):
    return np.clip(x, -1, 1)


def hardclipJ1(x):
    absed = np.abs(x)
    if absed < 1:
        return x * x / 2
    return absed - 1 / 2


def hardclipJ2(x):
    if np.abs(x) < 1:
        return x * x * x / 6
    return (x * x / 2 + 1 / 6) * np.sign(x) - (x / 2)


def hardclipAa0(x, gain=10):
    return hardclipJ0(gain * x)


def hardclipAa1(x, gain=10):
    return adaa1(gain * x, hardclipJ0, hardclipJ1)


def hardclipAa2(x, gain=10):
    return adaa2(gain * x, hardclipJ0, hardclipJ1, hardclipJ2)


def halfrectJ0(x):
    return np.maximum(0, x)


def halfrectJ1(x):
    if x < 0:
        return 0
    return x * x / 2


def halfrectJ2(x):
    if x < 0:
        return 0
    return x * x * x / 6


def halfrectAa0(x, gain=10):
    return (1 / gain) * halfrectJ0(gain * x)


def halfrectAa1(x, gain=10):
    return (1 / gain) * adaa1(gain * x, halfrectJ0, halfrectJ1)


def halfrectAa2(x, gain=10):
    return (1 / gain) * adaa2(gain * x, halfrectJ0, halfrectJ1, halfrectJ2)


def powerJ0(x, β=2.345):
    return np.sign(x) * np.power(np.abs(x), β)


def powerJ1(x, β=2.345):
    if β == -1:
        return np.log(np.abs(x))
    return np.power(np.abs(x), β + 1) / (β + 1)


def powerJ2(x, β=2.345):
    if β == -1:
        return x * (np.log(np.abs(x)) - 1)
    return np.sign(x) * np.power(np.abs(x), β + 2) / ((β + 3) * β + 2)


def powerAa0(x, gain=10, β=2.345):
    return (1 / np.power(gain, β)) * powerJ0(gain * x)


def powerAa1(x, gain=10, β=2.345):
    return (1 / np.power(gain, β)) * adaa1(gain * x, powerJ0, powerJ1)


def powerAa2(x, gain=10, β=2.345):
    return (1 / np.power(gain, β)) * adaa2(gain * x, powerJ0, powerJ1, powerJ2)


def softclip2J0(x, h=1, ratio=0.5):
    """`h` is threshold of clipping."""
    absed = np.abs(x)

    a1 = h * ratio
    if absed <= a1:
        return x

    a2 = 2 * h - a1
    if absed >= a2:
        return np.sign(x) * h

    C1 = a2 - absed
    return np.sign(x) * (h + 0.25 * C1 * C1 / (a1 - h))


def softclip2J1(x, h=1, ratio=0.5):
    absed = np.abs(x)

    a1 = h * ratio
    if absed <= a1:
        return absed * absed / 2

    a2 = 2 * h - a1
    C0 = a1 - a2
    if absed >= a2:
        return a1 * (a1 / 2 - h) + h * absed + C0 * C0 * C0 / (h - a1) / 12

    C1 = absed - a2
    return a1 * (a1 / 2 - h) + h * absed + (C0 * C0 * C0 - C1 * C1 * C1) / (h - a1) / 12


def softclip2J2(x, h=1, ratio=0.5):
    absed = np.abs(x)

    a1 = h * ratio
    if absed <= a1:
        return x * x * x / 6

    a2 = 2 * h - a1
    C0 = a1 - a2
    C1 = absed - a1
    if absed >= a2:
        return np.sign(x) * (
            a1 * a1 * (3 * absed - 2 * a1) / 6
            + C1 * C1 * h / 2
            + (C0 * C0 * C0 * (4 * absed - 3 * a1 - a2)) / (h - a1) / 48
        )

    C2 = absed + a1
    return np.sign(x) * (
        a1 * a1 * (3 * absed - 2 * a1) / 6
        + C1
        * C1
        * (h / 2 - (C2 * C2 + 2 * C0 * C0 - 4 * a2 * (C0 + absed)) / (h - a1) / 48)
    )


def softclip2Aa0(x, gain=10):
    y = np.zeros_like(x)
    for i, sig in enumerate(x):
        y[i] = softclip2J0(gain * sig)
    return y


def softclip2Aa1(x, gain=10):
    return adaa1(gain * x, softclip2J0, softclip2J1)


def softclip2Aa2(x, gain=10):
    return adaa2(gain * x, softclip2J0, softclip2J1, softclip2J2)


def softclipNJ0(x0, C=1, R=0.5, beta=2, S=0.1):
    absed = np.abs(x0)

    rc = C * R
    if absed <= rc:
        return x0

    xc = rc + beta * (C - rc)
    A = (rc - C) / (xc - rc) ** beta
    xs = xc - (-S / (A * beta)) ** (1 / (beta - 1))
    if absed < xs:
        return np.sign(x0) * (A * (xc - absed) ** beta + C)
    return np.sign(x0) * (A * (xc - xs) ** beta + C + S * (absed - xs))


def softclipNJ1(x0, C=1, R=0.5, beta=2, S=0.1):
    absed = np.abs(x0)

    rc = C * R
    if absed <= rc:
        return x0 * x0 / 2

    xc = rc + beta * (C - rc)
    Q0 = xc - rc
    A = (rc - C) / Q0**beta
    xs = xc - (-S / (A * beta)) ** (1 / (beta - 1))
    b1 = 1 + beta
    if absed < xs:
        return A * (Q0**b1 - (xc - absed) ** b1) / b1 + rc * rc / 2 + C * (absed - rc)
    return (
        A * Q0**b1 / b1
        + C * Q0
        + S * (absed * absed - xc * xc) / 2
        + rc * rc / 2
        + (absed - xc) * (A * (xc - xs) ** beta + C - S * xs)
    )


def softclipNJ2(x0, C=1, R=0.5, beta=2, S=0.1):
    absed = np.abs(x0)

    rc = C * R
    if absed <= rc:
        return x0 * x0 * x0 / 6

    xc = rc + beta * (C - rc)
    Q0 = xc - rc
    A = (rc - C) / Q0**beta
    xs = xc - (-S / (A * beta)) ** (1 / (beta - 1))
    b1 = 1 + beta
    b2 = 2 + beta
    if absed < xs:
        Q1 = absed - rc
        return np.sign(x0) * (
            A * (((xc - absed) ** b2 - Q0**b2) / (b1 * b2) + Q0**b1 * Q1 / b1)
            + rc * rc * (absed / 2 - rc / 3)
            + C * Q1 * Q1 / 2
        )
    Q2 = xc - xs
    return np.sign(x0) * (
        A * Q0**b2 * (1 - 1 / b2) / b1
        + C * Q0 * Q0 / 2
        + S * (absed * absed * absed - xc * xc * xc) / 6
        + rc * rc * (xc / 2 - rc / 3)
        + (absed - xc)
        * (
            A * (Q0**b1 / b1 - xc * Q2**beta)
            - C * rc
            + S * xc * (xs - xc / 2)
            + rc * rc / 2
            + (absed + xc) / 2 * (A * Q2**beta + C - S * xs)
        )
    )


def softclipNAa0(x, gain=10):
    y = np.zeros_like(x)
    for i, sig in enumerate(x):
        y[i] = softclipNJ0(gain * sig, 1, 0.5)
    return y


def softclipNAa1(x, gain=10):
    return adaa1(gain * x, softclipNJ0, softclipNJ1)


def softclipNAa2(x, gain=10):
    return adaa2(gain * x, softclipNJ0, softclipNJ1, softclipNJ2)


def tanhJ0(x):
    return np.tanh(x)


def tanhJ1(x):
    return np.log(np.cosh(x))


def tanhJ2(x):
    e2x = np.exp(2 * x)
    return x * (np.log(np.cosh(x) / (e2x + 1)) + x / 2) - special.spence(e2x) / 2


def tanhAa0(x, gain=10):
    return tanhJ0(gain * x)


def tanhAa1(x, gain=10):
    return adaa1(gain * x, tanhJ0, tanhJ1)


def tanhAa2(x, gain=10):
    return adaa2(gain * x, tanhJ0, tanhJ1, tanhJ2)


def atanJ0(x):
    return 2 / np.pi * np.arctan(x)


def atanJ1(x):
    return 2 / np.pi * (x * np.arctan(x) - np.log1p(x * x) / 2)


def atanJ2(x):
    return (x - x * np.log1p(x * x) + (x * x - 1) * np.arctan(x)) / np.pi


def atanAa0(x, gain=10):
    return atanJ0(gain * x)


def atanAa1(x, gain=10):
    return adaa1(gain * x, atanJ0, atanJ1)


def atanAa2(x, gain=10):
    return adaa2(gain * x, atanJ0, atanJ1, atanJ2)


def algebraicJ0(x):
    return x / (np.abs(x) + 1)


def algebraicJ1(x):
    z = np.abs(x)
    return z - np.log1p(z)


def algebraicJ2(x):
    z = np.abs(x)
    w = np.log1p(z)
    return np.sign(x) * (z * (1 + z / 2 - w) - w)


def algebraicAa0(x, gain=10):
    return algebraicJ0(gain * x)


def algebraicAa1(x, gain=10):
    return adaa1(gain * x, algebraicJ0, algebraicJ1)


def algebraicAa2(x, gain=10):
    return adaa2(gain * x, algebraicJ0, algebraicJ1, algebraicJ2)


def softplusJ0(x):
    return np.log(np.exp(x) + 1)


def softplusJ1(x):
    """
    `special.spence(0)` is a constant to compensate the difference between SciPy's `spence` and dilogarithm (Li_2). If Li_2 can be used, the term `special.spence(0)` should be removed.
    """
    return special.spence(0) - special.spence(np.exp(x))


def softplusJ2(x):
    """
    Slow because of arbitrary precision.
    It may cause glitch because of casting into `float`.
    """
    return -float(mpmath.polylog(3, -np.exp(x)))


def softplusAa0(x, gain=10):
    return 1 / gain * softplusJ0(gain * x)


def softplusAa1(x, gain=10):
    return 1 / gain * adaa1(gain * x, softplusJ0, softplusJ1)


def softplusAa2(x, gain=10):
    return 1 / gain * adaa2(gain * x, softplusJ0, softplusJ1, softplusJ2)


def swishJ0(x, β=2):
    if x == 0:
        return 1 / 2
    return x / (np.exp(-x * β) + 1)


def swishJ1(x, β=2):
    """
    `special.spence(0)` is a constant to compensate the difference between SciPy's `spence` and dilogarithm (Li_2). If Li_2 can be used, the term `special.spence(0)` should be removed.
    """
    exb = np.exp(x * β)
    return (x * β * np.log1p(exb) + special.spence(exb) - special.spence(0)) / (β * β)

    # # Reference implementation with mpmath.
    # return float(x * β * np.log1p(exb) + mpmath.polylog(2, -exb)) / (β * β)


def swishJ2(x, β=2):
    """
    Slow because of arbitrary precision.
    It may cause glitch because of casting into `float`.
    """
    exb = np.exp(x * β)
    return float(2 * mpmath.polylog(3, -exb) - x * β * mpmath.polylog(2, -exb)) / (
        β * β * β
    )


def swishAa0(x, gain=10):
    y = np.zeros_like(x)
    for i, sig in enumerate(x):
        y[i] = swishJ0(gain * sig)
    return y / gain


def swishAa1(x, gain=10):
    return 1 / gain * adaa1(gain * x, swishJ0, swishJ1)


def swishAa2(x, gain=10):
    return 1 / gain * adaa2(gain * x, swishJ0, swishJ1, swishJ2)


def exppolyJ0(x, β=2.5):
    absed = np.abs(x)
    return np.sign(x) * absed**β * np.exp(-absed)


def upperIncGamma(a, x):
    """
    `scipy.special.gammaincc` is regularized, which means it includes `1 / Γ(a)`.
    """
    return special.gamma(a) * special.gammaincc(a, x)


def exppolyJ1(x, β=2.5):
    b1 = β + 1
    C0 = upperIncGamma(b1, 0)  # Constant of integration.
    return C0 - upperIncGamma(b1, np.abs(x))

    # # Slightly efficient implementation. Leaving this for porting to other languages.
    # b1 = β + 1
    # return special.gamma(b1) * (
    #     special.gammaincc(b1, 0) - special.gammaincc(b1, np.abs(x))
    # )


def exppolyJ2(x, β=2.5):
    z = np.abs(x)
    b1 = β + 1
    b2 = β + 2
    C0 = upperIncGamma(b1, 0)
    C1 = upperIncGamma(b2, 0)  # Another constant of integration.
    return np.sign(x) * (upperIncGamma(b2, z) - C1 + z * (C0 - upperIncGamma(b1, z)))

    # # Slightly efficient implementation. Leaving this for porting to other languages.
    # z = np.abs(x)
    # b1 = β + 1
    # b2 = β + 2
    # reg1 = special.gamma(b1)
    # reg2 = special.gamma(b2)
    # return np.sign(x) * (
    #     reg2 * (special.gammaincc(b2, z) - special.gammaincc(b2, 0))
    #     + z * reg1 * (special.gammaincc(b1, 0) - special.gammaincc(b1, z))
    # )


def exppolyAa0(x, gain=10):
    return exppolyJ0(gain * x)


def exppolyAa1(x, gain=10):
    return adaa1(gain * x, exppolyJ0, exppolyJ1)


def exppolyAa2(x, gain=10):
    return adaa2(gain * x, exppolyJ0, exppolyJ1, exppolyJ2)


def cosdecayJ0(x):
    if x == 0:
        return 0
    return (1 - np.cos(x)) / x


def cosdecayJ1(x):
    if x == 0:
        return 0
    z = np.abs(x)
    _, ci = special.sici(z)
    return np.log(z) - ci


def cosdecayJ2(x):
    if x == 0:
        return 0
    z = np.abs(x)
    _, ci = special.sici(z)
    return np.sign(x) * (np.sin(z) + z * (np.log(z) - ci - 1))


def cosdecayAa0(x, gain=10):
    y = np.zeros_like(x)
    for i, sig in enumerate(x):
        y[i] = cosdecayJ0(gain * sig)
    return y


def cosdecayAa1(x, gain=10):
    return adaa1(gain * x, cosdecayJ0, cosdecayJ1)


def cosdecayAa2(x, gain=10):
    return adaa2(gain * x, cosdecayJ0, cosdecayJ1, cosdecayJ2)


def sinalgexpJ0(x):
    return np.sign(x) * np.sin(np.pi * (1 - np.exp(-np.abs(x))))


def sinalgexpJ1(x):
    C0, _ = special.sici(np.pi)
    si, _ = special.sici(np.pi * np.exp(-np.abs(x)))
    return C0 - si


def sinalgexpAa0(x, gain=10):
    return sinalgexpJ0(gain * x)


def sinalgexpAa1(x, gain=10):
    return adaa1(gain * x, sinalgexpJ0, sinalgexpJ1)


def sinatanexpJ0(x):
    return np.sign(x) * np.sin(2 * np.arctan(np.expm1(np.abs(x))))


def sinatanexpJ1(x):
    z = np.abs(x)
    em1x = np.expm1(z)
    return np.log(em1x * em1x + 1) / 2 + np.arctan(em1x) - z


def sinatanexpJ1Complex(x):
    z = np.abs(x)
    ex = np.exp(z)
    em1x = np.expm1(z)
    questionable = 1 + 1j  # No idea why this works. Shouldn't this be exp(i*π/4)?
    value = (
        questionable
        * (
            np.log(em1x * em1x + 1)
            - 2 * z
            - (1 - 1j) * np.arctan(1 - ex)
            + (1 + 1j) * np.arctan(ex / (2 - ex))
        )
        / 2
    )
    return value.real


def sinatanexpJ2Complex(x):
    z = np.abs(x)
    ex = np.exp(z)
    em1x = np.expm1(z)
    log_ex = np.log(ex)
    value = (
        -z * z
        + log_ex * np.log(em1x * em1x + 1)
        - (1 - 1j) * z * np.arctan(1 - ex)
        + (1 + 1j) * z * np.arctan(ex / (2 - ex))
        - (log_ex - 1j) * z * np.log(1 - (0.5 - 0.5j) * ex)
        - (log_ex + 1j) * z * np.log(1 - (0.5 + 0.5j) * ex)
        - (1 - 1j) * special.spence(-(0.5 - 0.5j) * ex)
        - (1 + 1j) * special.spence(-(0.5 + 0.5j) * ex)
    ) / 2
    return np.sign(x) * value.real


def sinatanexpAa0(x, gain=10):
    return sinatanexpJ0(gain * x)


def sinatanexpAa1(x, gain=10):
    return adaa1(gain * x, sinatanexpJ0, sinatanexpJ1)


def sinatanexpAa2(x, gain=10):
    return adaa2(gain * x, sinatanexpJ0, sinatanexpJ1, sinatanexpJ2Complex)


def sinrungeJ0(x):
    return np.sin(x) / (25 * x * x + 1)


def sinrungeJ1(x):
    si_neg, ci_neg = special.sici(1j / 5 - x)
    si_pos, ci_pos = special.sici(1j / 5 + x)
    value = 0.1 * (
        np.sinh(1 / 5) * (ci_neg + ci_pos) + 1j * np.cosh(1 / 5) * (si_neg + si_pos)
    )
    return value.real


def sinrungeJ2(x):
    si_neg, ci_neg = special.sici(1j / 5 - x)
    si_pos, ci_pos = special.sici(1j / 5 + x)

    # 100 * np.exp(1 / 5) = 0.0122140275816017
    # -1 + np.exp(2/5) = 0.49182469764127035
    # +1 + np.exp(2/5) = 2.4918246976412703
    value = (
        +(-1 + np.exp(2 / 5)) * ((-1j + 5 * x) * ci_neg + (+1j + 5 * x) * ci_pos)
        + (+1 + np.exp(2 / 5)) * ((+1 + 5j * x) * si_neg + (-1 + 5j * x) * si_pos)
    ) / (100 * np.exp(1 / 5))
    return value.real


def sinrungeAa0(x, gain=10):
    return 10 * sinrungeJ0(gain * x)


def sinrungeAa1(x, gain=10):
    return 10 * adaa1(gain * x, sinrungeJ0, sinrungeJ1)


def sinrungeAa2(x, gain=10):
    return 10 * adaa2(gain * x, sinrungeJ0, sinrungeJ1, sinrungeJ2)


def log1pJ0(x):
    return np.sign(x) * np.log1p(np.abs(x))


def log1pJ1(x):
    z = np.abs(x)
    return (z + 1) * np.log1p(z) - z


def log1pJ2(x):
    z = np.abs(x)
    z1 = z + 1
    return np.sign(x) * (2 * z1 * z1 * np.log1p(z) - (3 * z + 2) * z) / 4


def log1pAa0(x, gain=10):
    return 1 / np.log1p(gain) * log1pJ0(gain * x)


def log1pAa1(x, gain=10):
    return 1 / np.log1p(gain) * adaa1(gain * x, log1pJ0, log1pJ1)


def log1pAa2(x, gain=10):
    return 1 / np.log1p(gain) * adaa2(gain * x, log1pJ0, log1pJ1, log1pJ2)


# def sinrungeJ2(x):
#     if x == 0:
#         return 0
#     z = np.abs(x)
#     _, ci = special.sici(z)
#     return np.sign(x) * (np.sin(z) + z * (np.log(z) - ci - 1))


def render(frequencyHz=1661, durationSecond=1, gain=10, sampleRate=48000):
    mpmath.mp.prec = 256

    clipFunctions = [
        hardclipAa0,
        hardclipAa1,
        hardclipAa2,
        halfrectAa0,
        halfrectAa1,
        halfrectAa2,
        powerAa0,
        powerAa1,
        powerAa2,
        softclip2Aa0,
        softclip2Aa1,
        softclip2Aa2,
        softclipNAa0,
        softclipNAa1,
        softclipNAa2,
        tanhAa0,
        tanhAa1,
        tanhAa2,
        atanAa0,
        atanAa1,
        atanAa2,
        algebraicAa0,
        algebraicAa1,
        algebraicAa2,
        softplusAa0,
        softplusAa1,
        softplusAa2,  # Slow
        swishAa0,
        swishAa1,
        swishAa2,  # Slow
        exppolyAa0,
        exppolyAa1,
        exppolyAa2,
        cosdecayAa0,
        cosdecayAa1,
        cosdecayAa2,
        sinalgexpAa0,
        sinalgexpAa1,
        sinatanexpAa0,
        sinatanexpAa1,
        sinrungeAa0,
        sinrungeAa1,
        sinrungeAa2,
        log1pAa0,
        log1pAa1,
        log1pAa2,
    ]

    sndDir = Path("snd")
    sndDir.mkdir(parents=True, exist_ok=True)

    for clipFn in clipFunctions:
        sig = generateSine(int(sampleRate * durationSecond), frequencyHz / sampleRate)
        # sig *= signal.windows.tukey(len(sig), 0.1 / durationSecond)
        sig = clipFn(sig, gain)

        spc = np.fft.rfft(sig, norm="forward")
        absed = np.abs(spc)
        mag = 20 * np.log10(absed)
        freq = np.fft.rfftfreq(len(sig), 1 / sampleRate)

        soundfile.write(
            sndDir / Path(f"{clipFn.__name__}_{frequencyHz}Hz_gain{gain}.wav"),
            0.1 * sig,
            sampleRate,
            "FLOAT",
        )

        fig, ax = plt.subplots(2, 1)
        fig.set_size_inches(8, 6)
        fig.suptitle(f"{clipFn.__name__}")

        ax[0].plot(sig)
        ax[0].set_xlabel("Time [sample]")
        ax[0].set_ylabel("Amplitude")
        ax[0].set_xlim([0, 3 * int(sampleRate / frequencyHz)])

        ax[1].plot(freq, mag)
        # ax[1].set_xscale("log")
        ax[1].set_xlabel("Frequency [Hz]")
        ax[1].set_ylabel("Gain [dB]")
        ax[1].set_xlim([10, sampleRate // 2])
        ax[1].set_ylim([-120, 10])

        for axis in ax:
            axis.grid()
        plt.tight_layout()
        plt.show()


def testClipFunction():
    def plot(clipFuncs):
        def process(x, func):
            y = np.zeros_like(x)
            for i, sig in enumerate(x):
                y[i] = func(sig)
            return y

        x = np.linspace(-2, 2, 2048)
        for idx, fn in enumerate(clipFuncs):
            plt.plot(x, process(x, fn), label=fn.__name__)
        plt.grid()
        plt.legend()
        plt.show()

    plot([hardclipJ0, hardclipJ1, hardclipJ2])
    plot([halfrectJ0, halfrectJ1, halfrectJ2])
    plot([powerJ0, powerJ1, powerJ2])
    plot([softclip2J0, softclip2J1, softclip2J2])
    plot([softclipNJ0, softclipNJ1, softclipNJ2])
    plot([tanhJ0, tanhJ1, tanhJ2])
    plot([atanJ0, atanJ1, atanJ2])
    plot([algebraicJ0, algebraicJ1, algebraicJ2])
    plot([softplusJ0, softplusJ1, softplusJ2])
    plot([swishJ0, swishJ1, swishJ2])
    plot([exppolyJ0, exppolyJ1, exppolyJ2])
    plot([cosdecayJ0, cosdecayJ1, cosdecayJ2])
    plot([sinalgexpJ0, sinalgexpJ1])
    plot([sinatanexpJ0, sinatanexpJ1])
    plot([sinrungeJ0, sinrungeJ1, sinrungeJ2])
    plot([log1pJ0, log1pJ1, log1pJ2])


if __name__ == "__main__":
    render()
    testClipFunction()
