from mpmath import mp
import numpy as np

mp.dps = 50


def subnormal_near_zero(n=1000, dtype=np.float64):
    ss = np.finfo(dtype).smallest_subnormal
    return ss * np.arange(1, n + 1, dtype=dtype)


def subnormal_below_normal(n=1000, dtype=np.float64):
    sn = np.finfo(dtype).smallest_normal
    ss = np.finfo(dtype).smallest_subnormal
    return sn - ss * np.arange(0, n + 1, dtype=dtype)


def ulp_f32(x_exact):
    abs_x = abs(x_exact)
    if abs_x < mp.mpf(2) ** -126:
        return mp.mpf(2) ** -149
    else:
        e = mp.floor(mp.log(abs_x, 2))
        return mp.mpf(2) ** (e - 23)


def ulp_f64(x_exact):
    abs_x = abs(x_exact)
    if abs_x < mp.mpf(2) ** -1022:
        return mp.mpf(2) ** -1074
    else:
        e = mp.floor(mp.log(abs_x, 2))
        return mp.mpf(2) ** (e - 52)


def approx_f32(x):
    def fma_f32(a, b, c):
        with mp.workprec(50):
            exact = mp.mpf(float(a)) * mp.mpf(float(b)) + mp.mpf(float(c))
        with mp.workprec(24):
            f32_exact = mp.mpf(exact)
        return np.float32(float(f32_exact))

    # fmt: off
    c_f32 = [
        +1.551207740e-01,
        -1.196484253e+00,
        +5.100139453e+00,
        -1.033541937e+01,
        +6.283185274e+00,
    ]
    # fmt: on

    x_f32 = np.float32(x)
    t = np.float32(x_f32 * x_f32)
    res = np.float32(c_f32[0])
    res = fma_f32(res, t, c_f32[1])
    res = fma_f32(res, t, c_f32[2])
    res = fma_f32(res, t, c_f32[3])
    res = fma_f32(res, t, c_f32[4])
    return np.float32(x_f32 * res)


def approx_f64(x):
    def fma_f64(a, b, c):
        exact = mp.mpf(float(a)) * mp.mpf(float(b)) + mp.mpf(float(c))
        return np.float64(float(exact))

    # fmt: off
    c_f64 = [
        -4.22459059842197741e-05,
        +9.31955234064561289e-04,
        -1.47407221433722875e-02,
        +1.64291756513106096e-01,
        -1.19852905755141936e+00,
        +5.10032807971955293e+00,
        -1.03354255600995053e+01,
        +6.28318530717958534e+00,
    ]
    # fmt: on

    x_f64 = np.float64(x)
    t = np.float64(x_f64 * x_f64)
    res = np.float64(c_f64[0])
    res = fma_f64(res, t, c_f64[1])
    res = fma_f64(res, t, c_f64[2])
    res = fma_f64(res, t, c_f64[3])
    res = fma_f64(res, t, c_f64[4])
    res = fma_f64(res, t, c_f64[5])
    res = fma_f64(res, t, c_f64[6])
    res = fma_f64(res, t, c_f64[7])
    return np.float64(x_f64 * res)


def measure_ulp(dtype, approx_fn, ulp_fn):
    def run(case_name, xs):
        max_ulp = 0
        worst_x = np.nan
        for x in xs:
            y_exact = 2 * mp.sin(mp.pi * mp.mpf(float(x)))
            y_app = mp.mpf(float(approx_fn(x)))

            if y_exact == 0:
                ulp = mp.mpf(0)
            else:
                ulp_val = ulp_fn(y_exact)
                ulp = abs(y_app - y_exact) / ulp_val

            if ulp > max_ulp:
                max_ulp = ulp
                worst_x = x

        worst_x = worst_x if np.isfinite(worst_x) else "N/A"
        print(f"{case_name}, Max ULP: {float(max_ulp):.2e}, x: {worst_x}")

    print(f"--- {np.dtype(dtype).name}")
    rng = np.random.default_rng(628976)
    run(
        "linspace - [0, 0.5)",
        np.linspace(0, 0.5, 100000, dtype=dtype),
    )
    run(
        "Fuzzing  - [0, 0.5)",
        np.array(rng.uniform(0, 0.5, 1000000), dtype=dtype),
    )

    # # mpmath 1.4.1 has a bug on rounding mpf to subnormal floats.
    # run(
    #     "linspace - subnormals (count up)",
    #     subnormal_near_zero(100000, dtype=dtype),
    # )
    # run(
    #     "linspace - subnormals (count down)",
    #     subnormal_below_normal(100000, dtype=dtype),
    # )
    # run(
    #     "Fuzzing  - subnormals",
    #     np.array(rng.uniform(0, np.finfo(dtype).smallest_normal, 100000), dtype=dtype),
    # )


if __name__ == "__main__":
    measure_ulp(np.float32, approx_f32, ulp_f32)
    measure_ulp(np.float64, approx_f64, ulp_f64)
