from mpmath import mp
import numpy as np

mp.dps = 50


def sum_sqrt_mp_array(N: int, dtype=np.float64):
    y = [mp.mpf(0)]
    for i in range(1, N):
        y.append(y[-1] + mp.sqrt(i))
    return np.array([dtype(v) for v in y])


table_threshold = 64
table = sum_sqrt_mp_array(table_threshold)


def sum_of_sqrt_k2(N: int, dtype=np.float32):
    """f32 approximation (table size can be reduced)."""
    if N < table_threshold:
        return dtype(table[N])

    N = dtype(N)
    sqrt_n = np.sqrt(N, dtype=dtype)
    zeta_term = dtype(-0.2078862249773545660)  # zeta(-0.5)

    base = (dtype(2.0) / dtype(3.0) * N + dtype(0.5)) * sqrt_n + zeta_term

    m = dtype(1) / (N * N)

    v = dtype(1) / dtype(24)
    return base + (v / sqrt_n)


def sum_of_sqrt_k6(N: int, dtype=np.float64):
    """f64 approximation."""
    if N < table_threshold:
        return dtype(table[N])

    N = dtype(N)
    sqrt_n = np.sqrt(N, dtype=dtype)
    zeta_term = dtype(-0.2078862249773545660)  # zeta(-0.5)

    base = (dtype(2.0) / dtype(3.0) * N + dtype(0.5)) * sqrt_n + zeta_term

    m = dtype(1) / (N * N)

    v = dtype(1) / dtype(9216)
    v = v * m + dtype(-1) / dtype(1920)
    v = v * m + dtype(1) / dtype(24)

    return base + (v / sqrt_n)


_sum_sqrt_cache = [mp.mpf(0)]


def sum_of_sqrt_mp(N):
    if N < 823:
        curr_len = len(_sum_sqrt_cache)
        if N >= curr_len:
            last_val = _sum_sqrt_cache[-1]
            for i in range(curr_len, N + 1):
                last_val = last_val + mp.sqrt(mp.mpf(i))
                _sum_sqrt_cache.append(last_val)
        return _sum_sqrt_cache[N]

    N_mp = mp.mpf(N)
    sqrt_n = mp.sqrt(N_mp)
    zeta_const = mp.mpf("-0.2078862249773545660173067253970493022262685312876725376")

    base = (mp.mpf("2") / mp.mpf("3") * N_mp + mp.mpf("0.5")) * sqrt_n + zeta_const

    m = mp.mpf("1") / (N_mp * N_mp)

    v = mp.mpf("0.01257060903645608396757216680617559523809523809523809524")
    v = v * m + mp.mpf("-0.0022081132046878337860107421875")
    v = v * m + mp.mpf("0.0005166033903757731119791666666666666666666666666666666667")
    v = v * m + mp.mpf("-0.0001689312950013175843253968253968253968253968253968253968")
    v = v * m + mp.mpf("0.00008265177408854166666666666666666666666666666666666666667")
    v = v * m + mp.mpf("-0.000067138671875")
    v = v * m + mp.mpf("0.0001085069444444444444444444444444444444444444444444444444")
    v = v * m + mp.mpf("-0.0005208333333333333333333333333333333333333333333333333333")
    v = v * m + mp.mpf("0.04166666666666666666666666666666666666666666666666666667")

    return base + (v / sqrt_n)


def relative_error(ref, approx, dtype=np.float64):
    ref_arr = np.asarray(ref, dtype=dtype)
    approx_arr = np.asarray(approx, dtype=dtype)
    with np.errstate(divide="ignore", invalid="ignore"):
        error = np.abs((ref_arr - approx_arr) / ref_arr, dtype=dtype)
    return np.nan_to_num(error, nan=0.0)


def ulp_error(ref, approx, dtype=np.float64):
    approx = np.asanyarray(approx, dtype=dtype)
    ref = np.asanyarray(ref, dtype=dtype)

    eps = np.spacing(ref, dtype=dtype)
    with np.errstate(divide="ignore", invalid="ignore"):
        error = np.abs(approx - ref, dtype=dtype) / eps

    error = np.where(approx == ref, 0.0, error)
    error = np.where(np.isinf(ref) & (approx != ref), np.inf, error)

    return error.item() if error.ndim == 0 else error


def test_from_0(dtype=np.float64):
    N = 2**20
    ref = sum_sqrt_mp_array(N)
    approx = [sum_of_sqrt_k6(n, dtype) for n in range(N)]

    rel = relative_error(ref, approx, dtype=dtype)
    ulp = ulp_error(ref, approx, dtype=dtype)

    def print_max_error(title, error):
        i = np.argmax(ulp)
        print(f"{title:4} {i:8} {error[i]:.2e} {ref[i]:.17e} {approx[i]:.17e}")

    print(f"{"":4} {"N":8} {"Error":8} {"Ref.":22} {"Approx.":22}")
    print_max_error("ULP", ulp)
    print_max_error("Rel", rel)


def test_large(dtype=np.float64):
    upper_bound = int(1 / np.finfo(dtype).eps) - 1
    width = 2
    N = np.arange(upper_bound - width, upper_bound + width)

    ref = [sum_of_sqrt_mp(n) for n in N]
    approx = [sum_of_sqrt_k6(n, dtype) for n in N]

    for n, a in zip(N, approx):
        print(n, a)

    rel = relative_error(ref, approx, dtype=dtype)
    ulp = ulp_error(ref, approx, dtype=dtype)

    def print_max_error(title, error):
        i = np.argmax(ulp)
        print(
            f"{title:4} {N[i]:16} {error[i]:10.2e} {ref[i]:24.17e} {approx[i]:24.17e}"
        )

    print(f"{"":4} {"N":16} {"Error":11} {"Ref.":24} {"Approx.":24}")
    print_max_error("ULP", ulp)
    print_max_error("Rel", rel)


def test_single():
    N = np.int64(4503599627370496)
    print(sum_of_sqrt_k2(N))
    print(sum_of_sqrt_k6(N))
    print(float(sum_of_sqrt_mp(N)))


test_single()
# test_from_0()
# test_large()
