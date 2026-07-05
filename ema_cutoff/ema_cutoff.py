import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from mpmath import mp
from numpy.typing import ArrayLike, NDArray, DTypeLike as Dtl
from typing import Any, Union
import math

mp.dps = 50


def round_to_f64(x):
    x = mp.mpf(x)
    if mp.isnan(x):
        return np.float64(np.nan)
    if mp.isinf(x):
        return np.float64(np.inf) if x > 0 else np.float64(-np.inf)
    if x == 0:
        return np.float64(0.0)

    # Extract sign
    sign = 1
    if x < 0:
        sign = -1
        x = -x

    e = int(mp.floor(mp.log(x, 2)))

    # Overflow check: round-to-nearest-even cutoff to infinity is 2**1024 - 2**970
    if e >= 1024 or (e == 1023 and x >= mp.mpf(2) ** 1024 - mp.mpf(2) ** 970):
        return np.float64(np.inf) if sign > 0 else np.float64(-np.inf)

    if e < -1022:
        # Subnormal
        r = x / mp.mpf(2) ** -1074
        r_floor = mp.floor(r)
        diff = r - r_floor
        if diff < 0.5:
            r_rounded = r_floor
        elif diff > 0.5:
            r_rounded = r_floor + 1
        else:  # exactly 0.5
            if r_floor % 2 == 0:
                r_rounded = r_floor
            else:
                r_rounded = r_floor + 1

        bits = int(r_rounded)
        res = np.uint64(bits).view(np.float64)
        return res if sign > 0 else -res
    else:
        # Normal
        shift = e - 52
        r = x / (mp.mpf(2) ** shift)

        r_floor = mp.floor(r)
        diff = r - r_floor
        if diff < 0.5:
            r_rounded = r_floor
        elif diff > 0.5:
            r_rounded = r_floor + 1
        else:  # exactly 0.5
            if r_floor % 2 == 0:
                r_rounded = r_floor
            else:
                r_rounded = r_floor + 1

        if r_rounded >= 2**53:
            r_rounded = 2**52
            e += 1
            if e >= 1024:
                return np.float64(np.inf) if sign > 0 else np.float64(-np.inf)

        biased_exp = e + 1023
        mantissa = int(r_rounded) - 2**52

        # Build bits with sign bit (bit 63)
        sign_bit = (1 << 63) if sign < 0 else 0
        bits = sign_bit | (biased_exp << 52) | mantissa
        return np.uint64(bits).view(np.float64)


def round_to_f32(x):
    x = mp.mpf(x)
    if mp.isnan(x):
        return np.float32(np.nan)
    if mp.isinf(x):
        return np.float32(np.inf) if x > 0 else np.float32(-np.inf)
    if x == 0:
        return np.float32(0.0)

    # Extract sign
    sign = 1
    if x < 0:
        sign = -1
        x = -x

    e = int(mp.floor(mp.log(x, 2)))

    # Overflow check: round-to-nearest-even cutoff to infinity is 2**128 - 2**103
    if e >= 128 or (e == 127 and x >= mp.mpf(2) ** 128 - mp.mpf(2) ** 103):
        return np.float32(np.inf) if sign > 0 else np.float32(-np.inf)

    if e < -126:
        # Subnormal
        r = x / mp.mpf(2) ** -149
        r_floor = mp.floor(r)
        diff = r - r_floor
        if diff < 0.5:
            r_rounded = r_floor
        elif diff > 0.5:
            r_rounded = r_floor + 1
        else:
            if r_floor % 2 == 0:
                r_rounded = r_floor
            else:
                r_rounded = r_floor + 1

        bits = int(r_rounded)
        res = np.uint32(bits).view(np.float32)
        return res if sign > 0 else -res
    else:
        # Normal
        shift = e - 23
        r = x / (mp.mpf(2) ** shift)

        r_floor = mp.floor(r)
        diff = r - r_floor
        if diff < 0.5:
            r_rounded = r_floor
        elif diff > 0.5:
            r_rounded = r_floor + 1
        else:
            if r_floor % 2 == 0:
                r_rounded = r_floor
            else:
                r_rounded = r_floor + 1

        if r_rounded >= 2**24:
            r_rounded = 2**23
            e += 1
            if e >= 128:
                return np.float32(np.inf) if sign > 0 else np.float32(-np.inf)

        biased_exp = e + 127
        mantissa = int(r_rounded) - 2**23

        # Build bits with sign bit (bit 31)
        sign_bit = (1 << 31) if sign < 0 else 0
        bits = sign_bit | (biased_exp << 23) | mantissa
        return np.uint32(bits).view(np.float32)


def to_float(x, dtype: Dtl = np.float64) -> float:
    if dtype == np.float32:
        return float(round_to_f32(x))
    return float(round_to_f64(x))


def cutoffToAlphaReference(cutoffNormalized: float, dtype: Dtl = np.float64) -> mp.mpf:
    with mp.workdps(720):
        cn = mp.mpf(float(dtype(cutoffNormalized)))
        y = 1 - mp.cos(2 * mp.pi * cn)
        result = mp.sqrt((y + 2) * y) - y
        return result


def cutoffToAlphaNaive(cutoffNormalized: float, dtype: Dtl = np.float64) -> float:
    y = 1 - np.cos(2 * np.pi * cutoffNormalized, dtype=dtype)
    result = np.sqrt((y + 2) * y, dtype=dtype) - y
    return float(result)


def cutoffToAlphaA(cutoffNormalized: float, dtype: Dtl = np.float64) -> float:
    """
    `cutoffNormalized < 1e-37` is required to avoid division by 0.
    The branch is commented out for the comparison.
    Accuracy degrades when `cutoffNormalized < 1e-38`.
    """
    # if cutoffNormalized < 1e-37:
    #     return float(0)
    sn = np.sin(np.pi * cutoffNormalized, dtype=dtype)
    y = 2 * sn * sn
    result = 2 * y / (np.sqrt(y * y + 2 * y, dtype=dtype) + y)
    return float(result)


def cutoffToAlphaB(cutoffNormalized: float, dtype: Dtl = np.float64) -> float:
    sn = np.sin(np.pi * cutoffNormalized, dtype=dtype)
    y = 2 * sn * sn
    result = 2 / (np.sqrt(1 + 2 / y, dtype=dtype) + 1)
    return float(result)


def cutoffToAlphaC(cutoffNormalized: float, dtype: Dtl = np.float64) -> float:
    """
    Normal numbers are accurate.
    """
    sn = np.sin(np.pi * cutoffNormalized, dtype=dtype)
    result = 2 * sn / (np.sqrt(sn * sn + 1, dtype=dtype) + sn)
    return float(result)


def cutoffToAlphaC_fma(cutoffNormalized: float, dtype: Dtl = np.float64) -> float:
    sn = np.sin(np.pi * cutoffNormalized, dtype=dtype)
    sn2p1 = math.fma(sn, sn, 1.0)
    result = 2 * sn / (np.sqrt(sn2p1, dtype=dtype) + sn)
    return float(result)


def cutoffToAlphaD(cutoffNormalized: float, dtype: Dtl = np.float64) -> float:
    """
    Subnormals are handled.
    """
    pi = dtype(np.pi)
    twopi = dtype(2.0 * pi)
    cn = dtype(cutoffNormalized)

    omega = twopi * cn
    if omega < np.finfo(dtype).eps:
        return float(omega)

    sn = np.sin(pi * cn, dtype=dtype)
    result = 2 * sn / (np.sqrt(sn * sn + 1, dtype=dtype) + sn)
    return float(result)


def ulp_error(ref, approx, dtype=np.float64):
    approx = dtype(approx)

    if approx == ref:
        return 0.0

    if math.isnan(approx) or mp.isnan(ref):
        return math.nan

    ref_flt = to_float(ref, dtype)
    if math.isinf(approx) or math.isinf(ref):
        return math.inf

    eps = abs(np.spacing(ref_flt, dtype=dtype))
    ulp = abs(approx - ref) / eps
    return to_float(ulp, dtype)


def relative_error(ref, approx, signed_zero=True, dtype=np.float64):
    approx = dtype(approx)

    if math.isnan(approx) and mp.isnan(ref):
        return 0.0
    if math.isinf(approx) or mp.isinf(ref):
        return 0.0 if to_float(ref, dtype) == approx else 1.0

    if ref == 0.0:
        if approx != 0.0:
            return 1.0
        if signed_zero and math.copysign(1.0, ref) != math.copysign(1.0, approx):
            return 1.0  # Mismatching signs on 0.
        return 0.0

    diff = ref - approx
    error = abs(diff / ref)
    return to_float(error, dtype)


def ulp_error_array(ref, approx, dtype=np.float64):
    return [ulp_error(p, q, dtype) for p, q in zip(ref, approx)]


def relative_error_array(ref, approx, signed_zero=True, dtype=np.float64):
    return [relative_error(p, q, signed_zero, dtype) for p, q in zip(ref, approx)]


def plotError(x=None, data=None, savefig=""):
    if x is None:
        # x = np.geomspace(1e-7, 0.5, 8193)
        x = np.geomspace(1e-53, 1e-7, 8193)
        # x = np.geomspace(np.finfo(np.float64).smallest_subnormal, 1e-53, 8193)
    ref_f32 = [cutoffToAlphaReference(fc, np.float32) for fc in x]
    ref_f64 = [cutoffToAlphaReference(fc, np.float64) for fc in x]

    print(f"--- Value at f_c = 0\nref:  {cutoffToAlphaReference(0)}")

    if data is None:
        data = {
            # "n32": lambda f: cutoffToAlphaNaive(f, np.float32),
            # "n64": lambda f: cutoffToAlphaNaive(f, np.float64),
            # "a32": lambda f: cutoffToAlphaA(f, np.float32),
            # "a64": lambda f: cutoffToAlphaA(f, np.float64),
            # "b32": lambda f: cutoffToAlphaB(f, np.float32),
            # "b64": lambda f: cutoffToAlphaB(f, np.float64),
            # "c32": lambda f: cutoffToAlphaC(f, np.float32),
            # "c64": lambda f: cutoffToAlphaC(f, np.float64),
            "d32": lambda f: cutoffToAlphaD(f, np.float32),
            "d64": lambda f: cutoffToAlphaD(f, np.float64),
        }

    plt.figure(figsize=(8, 4))
    cmap = plt.get_cmap("plasma")
    for idx, (key, func) in enumerate(data.items()):
        dtype = np.float32 if "32" in key else np.float64
        ref = ref_f32 if dtype == np.float32 else ref_f64
        print(f"{key}:  {func(0)}")

        y = [func(fc) for fc in x]
        relativeError = relative_error_array(ref, y, dtype)
        plt.plot(
            x, relativeError, color=cmap(idx / len(data)), alpha=0.5, lw=1, label=key
        )

    plt.axhline(
        np.finfo(np.float32).eps,
        color="black",
        lw=1,
        ls="--",
        label="f32 eps.",
    )
    plt.axhline(
        np.finfo(np.float64).eps,
        color="black",
        lw=1,
        ls="-.",
        label="f64 eps.",
    )

    sn32 = np.finfo(np.float32).smallest_normal
    if min(x) <= sn32 and max(x) >= sn32:
        plt.axvline(sn32, color="yellow", lw=1, ls="--", label="f32 tiny")
    ss32 = np.finfo(np.float32).smallest_subnormal
    if min(x) <= ss32 and max(x) >= ss32:
        plt.axvline(ss32, color="orange", lw=1, ls="-.", label="f32 min.")
    sn64 = np.finfo(np.float64).smallest_normal
    if min(x) <= sn64:
        plt.axvline(sn64, color="yellow", lw=1, ls="--", label="f64 tiny")
    ss64 = np.finfo(np.float64).smallest_subnormal
    if min(x) <= ss64:
        plt.axvline(ss64, color="orange", lw=1, ls="-.", label="f64 min.")

    plt.title("Relative Errors of EMA Filter Cutoff to α")
    plt.xlabel("Normalized Cutoff [rad/2π]")
    plt.ylabel("Relative Error")
    plt.xscale("log")
    plt.yscale("log")
    # plt.ylim([1e-17, 1])
    plt.yticks(np.geomspace(1e-16, 1, 9))
    plt.legend()
    plt.grid(color="#f4f4f4")
    plt.tight_layout()
    if len(savefig) > 0:
        plt.savefig(f"{savefig}.svg")
    plt.show()


def listSubnormalsFrom0(n: int, start: float = 0.0, dtype: Dtl = np.float64):
    def generate(n: int, start: float = 0.0, dtype: Dtl = np.float64):
        dtype = np.dtype(dtype)
        v = dtype.type(start)
        target = dtype.type(np.inf)
        for _ in range(n):
            v = np.nextafter(v, target)
            yield v

    return np.array(list(generate(n, start, dtype)), dtype=dtype)


def plotUlp(target_fn=cutoffToAlphaD, savefig=""):
    def plotAxis(x, axis, dtype, fn, annotation):
        y = np.array([fn(fc, dtype) for fc in x])
        ref = [cutoffToAlphaReference(fc, dtype) for fc in x]

        bits = "32" if dtype == np.float32 else "64"
        axis.set_title(f"f{bits} Subnormals - {annotation}")

        scalar = 2**30 if dtype == np.float32 else 2**300
        errors = ulp_error_array(ref, y, dtype)
        axis.plot(x * scalar, errors, color="black", alpha=0.5, lw=1)
        axis.axhline(0.5, color="red", lw=1, ls="--", label="0.5 ULP")
        axis.set_xlabel("Normalized Cutoff [rad/2π]")
        axis.set_ylabel("ULP Error")
        axis.grid(color="#f4f4f4")

        def subnormal_formatter(x, pos):
            if pos is not None and pos % 2 == 0:
                return ""
            return f"{x / scalar:.1e}"

        axis.xaxis.set_major_formatter(FuncFormatter(subnormal_formatter))
        # axis.set_xscale("log")
        # axis.set_yscale("log")

    def x_subrnomal_linspace(dtype, size=2048):
        sn = np.finfo(dtype).smallest_normal
        ss = np.finfo(dtype).smallest_subnormal
        return np.linspace(ss, sn, size, dtype=dtype)

    x_lin_f32 = x_subrnomal_linspace(np.float32)
    x_lin_f64 = x_subrnomal_linspace(np.float64)
    x_from0_f32 = listSubnormalsFrom0(2048, dtype=np.float32)
    x_from0_f64 = listSubnormalsFrom0(2048, dtype=np.float64)

    fig, ax = plt.subplots(2, 2)
    fig.set_size_inches(8, 6)
    fig.suptitle(f"Subnormal ULP Errors on {target_fn.__name__} (Red line is 0.5 ULP)")
    plotAxis(x_lin_f32, ax[0][0], np.float32, target_fn, "linspace")
    plotAxis(x_lin_f64, ax[0][1], np.float64, target_fn, "linspace")
    plotAxis(x_from0_f32, ax[1][0], np.float32, target_fn, "from 0")
    plotAxis(x_from0_f64, ax[1][1], np.float64, target_fn, "from 0")
    fig.tight_layout()
    if len(savefig) >= 1:
        plt.savefig(f"{savefig}.svg")
    plt.show()


def testSingle():
    dtype = np.float32
    # x = 1e-323
    # x = 2.2826844664704e-310
    # x = 2.39138372677854e-310
    # x = 2.1163745981989623e-308
    x = 1.1675765e-41  # Known smallest f32 that exceeds 0.5 ULP.
    # x = 5.48450035e-316  # Known smallest f64 that exceeds 0.5 ULP.

    r = cutoffToAlphaReference(x, dtype)
    a = cutoffToAlphaD(x, dtype)
    ulp = ulp_error(r, a, dtype)

    decimal_digits = 9 if dtype == np.float32 else 17
    print(f"{r:.{decimal_digits}e}")
    print(f"{a:.{decimal_digits}e}")
    print(ulp)


def sweep():
    dtype = np.float64
    sn = np.finfo(dtype).smallest_normal
    ss = np.finfo(dtype).smallest_subnormal

    x = np.linspace(ss, sn, 2048)
    # x = listSubnormalsFrom0(100000, dtype=dtype)

    for v in x:
        v = float(v)
        r = cutoffToAlphaReference(v)
        a = cutoffToAlphaD(v)
        ulp = ulp_error(r, a, dtype)
        if ulp > 0:
            rel = relative_error(r, a, dtype)
            if rel >= np.finfo(dtype).eps:
                print(v, r, a, ulp, rel)
    # print(v)


def renderPlots():
    x1 = np.geomspace(1e-7, 0.5, 8193)
    x2 = np.geomspace(1e-53, 1e-7, 8193)
    x3 = np.geomspace(np.finfo(np.float64).smallest_subnormal, 1e-53, 8193)

    naive_fn = {
        "n32": lambda f: cutoffToAlphaNaive(f, np.float32),
        "n64": lambda f: cutoffToAlphaNaive(f, np.float64),
    }
    plotError(x1, naive_fn, "img/naive_relative_errors")

    c_fn = {
        "c32": lambda f: cutoffToAlphaC(f, np.float32),
        "c64": lambda f: cutoffToAlphaC(f, np.float64),
    }
    plotError(x1, c_fn, "img/c_relative_errors_small")
    plotError(x2, c_fn, "img/c_relative_errors_small_f32")
    plotError(x3, c_fn, "img/c_relative_errors_small_f64")

    d_fn = {
        "d32": lambda f: cutoffToAlphaD(f, np.float32),
        "d64": lambda f: cutoffToAlphaD(f, np.float64),
    }
    plotError(x1, d_fn, "img/d_relative_errors_small")
    plotError(x2, d_fn, "img/d_relative_errors_small_f32")
    plotError(x3, d_fn, "img/d_relative_errors_small_f64")

    plotUlp(cutoffToAlphaC, "img/c_ulp_errors_subnormal")
    plotUlp(cutoffToAlphaD, "img/d_ulp_errors_subnormal")


if __name__ == "__main__":
    # renderPlots()
    # plotError()
    # plotUlp()
    # sweep()
    testSingle()
