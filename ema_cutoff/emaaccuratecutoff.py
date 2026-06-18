import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp
from numpy.typing import ArrayLike, NDArray, DTypeLike as Dtl
from typing import Any, Union
import math


def round_to_f64(x):
    if x == 0:
        return np.float64(0.0)

    e = int(mp.floor(mp.log(x, 2)))

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
        return np.uint64(bits).view(np.float64)
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
            shift = e - 52

        biased_exp = e + 1023
        mantissa = int(r_rounded) - 2**52
        bits = (biased_exp << 52) | mantissa
        return np.uint64(bits).view(np.float64)


def round_to_f32(x):
    if x == 0:
        return np.float32(0.0)

    e = int(mp.floor(mp.log(x, 2)))

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
        return np.uint32(bits).view(np.float32)
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

        biased_exp = e + 127
        mantissa = int(r_rounded) - 2**23
        bits = (biased_exp << 23) | mantissa
        return np.uint32(bits).view(np.float32)


def to_float(x, dtype: Dtl = np.float64) -> float:
    if dtype == np.float32:
        return float(round_to_f32(x))
    return float(round_to_f64(x))


def cutoffToAlphaReference(cutoffNormalized: float, dtype: Dtl = np.float64) -> float:
    with mp.workdps(720):
        cn = mp.mpf(float(dtype(cutoffNormalized)))
        y = 1 - mp.cos(2 * mp.pi * cn)
        result = mp.sqrt((y + 2) * y) - y
        return float(dtype(result))
        # return to_float(result)


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


def relative_error_array(
    ref: ArrayLike, approx: ArrayLike, dtype: Dtl = np.float64
) -> NDArray[Any]:
    ref_arr = np.asarray(ref, dtype=dtype)
    approx_arr = np.asarray(approx, dtype=dtype)
    with np.errstate(divide="ignore", invalid="ignore"):
        error = np.abs((ref_arr - approx_arr) / ref_arr, dtype=dtype)
    return np.nan_to_num(error, nan=0.0)


def ulp_error_array(
    ref: ArrayLike, approx: ArrayLike, dtype: Dtl = np.float64
) -> Union[NDArray[Any], Any]:
    approx = np.asanyarray(approx, dtype=dtype)
    ref = np.asanyarray(ref, dtype=dtype)

    eps = np.spacing(ref, dtype=dtype)
    with np.errstate(divide="ignore", invalid="ignore"):
        error = np.abs(approx - ref, dtype=dtype) / eps

    error = np.where(approx == ref, 0.0, error)
    error = np.where(np.isinf(ref) & (approx != ref), np.inf, error)

    return error.item() if error.ndim == 0 else error


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


def plotUlp():
    def plotAxis(axis, dtype, fn, prefix):
        sn = np.finfo(dtype).smallest_normal
        ss = np.finfo(dtype).smallest_subnormal

        # x = np.linspace(ss, sn, 2048, dtype=dtype)
        x = listSubnormalsFrom0(2048, dtype=dtype)
        y = np.array([fn(fc, dtype) for fc in x])
        ref = np.array([cutoffToAlphaReference(fc, dtype) for fc in x])

        bits = "32" if dtype == np.float32 else "64"
        axis.set_title(f"f{bits} Subnormals ({prefix}{bits})")

        scalar = 2**30 if dtype == np.float32 else 2**300
        relativeError = ulp_error_array(ref, y, dtype)
        # relativeError = relative_error_array(ref, y, dtype)
        axis.plot(x * scalar, relativeError, color="black", alpha=0.5, lw=1)
        # axis.plot(x * scalar, (ref - y) * scalar, color="red", alpha=0.5, lw=1)
        axis.set_xlabel("Normalized Cutoff [rad/2π]")
        axis.set_ylabel("ULP Error")
        axis.grid(color="#f4f4f4")
        # axis.set_xscale("log")
        # axis.set_yscale("log")

    prefix = "c"
    target_fn = cutoffToAlphaC
    # prefix = "d"
    # target_fn = cutoffToAlphaD

    fig, ax = plt.subplots(1, 2)
    fig.set_size_inches(8, 4)
    plotAxis(ax[0], np.float32, target_fn, prefix)
    plotAxis(ax[1], np.float64, target_fn, prefix)
    fig.tight_layout()
    plt.show()


def testSingle():
    dtype = np.float32
    # x = 1e-323
    # x = 2.2826844664704e-310
    # x = 2.39138372677854e-310
    # x = 2.1163745981989623e-308
    x = 1.1675765e-41  # Known smallest f32 that produces an error.
    # x = 5.48450035e-316  # Known smallest f64 that produces an error.

    r = cutoffToAlphaReference(x, dtype)
    a = cutoffToAlphaD(x, dtype)
    ulp = ulp_error_array(r, a, dtype)

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
        ulp = ulp_error_array(r, a)
        if ulp > 0:
            rel = relative_error_array(r, a)
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


if __name__ == "__main__":
    # renderPlots()
    # plotError()
    # plotUlp()
    # sweep()
    testSingle()
