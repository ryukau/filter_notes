"""
- `butter_zpk_mp` is a port of `scipy.signal.butter(..., output="zpk")`.
- `ellip_zpk_mp` is a port of `scipy.signal.ellip(..., output="zpk")`.
- `zpk2sos_mp`
- `tf2sos_mp`
"""

import collections.abc
import mpmath


def product_mp(lst):
    res = mpmath.mpf(1.0)
    for x in lst:
        res *= x
    return res


def buttap_mp(N):
    z = []
    p = []
    for m in range(-N + 1, N, 2):
        theta = mpmath.pi * m / (2 * N)
        pole = -mpmath.exp(1j * theta)
        p.append(pole)
    k = mpmath.mpf(1.0)
    return z, p, k


def _ellipdeg_mp(n, m1):
    """Solve the degree equation using the nome series representation."""
    n = mpmath.mpf(n)
    m1 = mpmath.mpf(m1)
    K1 = mpmath.ellipk(m1)
    K1p = mpmath.ellipk(1 - m1)

    q1 = mpmath.exp(-mpmath.pi * K1p / K1)
    q = q1 ** (1 / n)

    eps = mpmath.mp.eps
    num = mpmath.mpf(0.0)
    r = 0
    while True:
        term = q ** (r * (r + 1))
        num += term
        if term < eps * num and r > 5:
            break
        r += 1

    den = mpmath.mpf(1.0)
    r = 1
    while True:
        term = 2 * (q ** (r**2))
        den += term
        if term < eps * den and r > 5:
            break
        r += 1

    return 16 * q * (num / den) ** 4


def _arc_jac_sn_mp(w, m):
    """Inverse Jacobian elliptic sn function."""
    return mpmath.ellipf(mpmath.asin(w), m)


def _arc_jac_sc1_mp(w, m):
    """Real inverse Jacobian sc function with complementary modulus."""
    zcomplex = _arc_jac_sn_mp(1j * w, m)
    return mpmath.im(zcomplex)


def ellipap_mp(N, rp, rs):
    """Return zeros, poles, and gain of an Nth-order elliptic analog lowpass prototype."""
    if int(N) != N or N < 0:
        raise ValueError("Filter order must be a nonnegative integer")

    rp = mpmath.mpf(rp)
    rs = mpmath.mpf(rs)

    if N == 0:
        return [], [], 10 ** (-rp / 20)
    elif N == 1:
        pow10m1_rp = mpmath.powm1(10, 0.1 * rp)
        p = -mpmath.sqrt(1.0 / pow10m1_rp)
        k = -p
        z = []
        return z, [p], k

    eps_sq = mpmath.powm1(10, 0.1 * rp)
    eps = mpmath.sqrt(eps_sq)
    ck1_sq = eps_sq / mpmath.powm1(10, 0.1 * rs)
    if ck1_sq == 0:
        raise ValueError(
            "Cannot design a filter with the given rp and rs specifications."
        )

    val = (mpmath.ellipk(ck1_sq), mpmath.ellipk(1 - ck1_sq))
    m = _ellipdeg_mp(N, ck1_sq)
    capk = mpmath.ellipk(m)

    j_list = list(range(1 - N % 2, N, 2))

    s = []
    c = []
    d = []
    for j_val in j_list:
        u_val = mpmath.mpf(j_val) * capk / N
        s.append(mpmath.ellipfun("sn", u_val, m))
        c.append(mpmath.ellipfun("cn", u_val, m))
        d.append(mpmath.ellipfun("dn", u_val, m))

    # Threshold small values of sn near zero to avoid division issues
    snew = [val_ for val_ in s if abs(val_) > mpmath.mp.eps]
    sqrt_m = mpmath.sqrt(m)
    z = []
    for val_ in snew:
        z_val = 1j / (sqrt_m * val_)
        z.append(z_val)
    z = z + [mpmath.conj(val_) for val_ in z]

    r_val = _arc_jac_sc1_mp(1.0 / eps, ck1_sq)
    v0 = capk * r_val / (N * val[0])

    sv = mpmath.ellipfun("sn", v0, 1 - m)
    cv = mpmath.ellipfun("cn", v0, 1 - m)
    dv = mpmath.ellipfun("dn", v0, 1 - m)

    p = []
    for s_val, c_val, d_val in zip(s, c, d):
        num_val = -(c_val * d_val * sv * cv + 1j * s_val * dv)
        den_val = 1.0 - (d_val * sv) ** 2.0
        p_val = num_val / den_val
        p.append(p_val)

    if N % 2 != 0:
        # The first pole on the real axis has imaginary part of zero
        p[0] = mpmath.re(p[0])
        p = p + [mpmath.conj(x) for x in p[1:]]
    else:
        p = p + [mpmath.conj(x) for x in p]

    # Calculate system gain
    prod_p = product_mp([-x for x in p])
    prod_z = product_mp([-x for x in z])
    k = mpmath.re(prod_p / prod_z)
    if N % 2 == 0:
        k = k / mpmath.sqrt(1 + eps_sq)

    return z, p, k


def lp2lp_zpk_mp(z, p, k, wo=1.0):
    wo = mpmath.mpf(wo)
    degree = len(p) - len(z)
    z_lp = [wo * x for x in z]
    p_lp = [wo * x for x in p]
    k_lp = k * (wo**degree)
    return z_lp, p_lp, k_lp


def lp2hp_zpk_mp(z, p, k, wo=1.0):
    wo = mpmath.mpf(wo)
    degree = len(p) - len(z)
    z_hp = [wo / x for x in z]
    p_hp = [wo / x for x in p]
    z_hp = z_hp + [mpmath.mpf(0.0)] * degree

    prod_z = product_mp([-x for x in z])
    prod_p = product_mp([-x for x in p])
    k_hp = k * mpmath.re(prod_z / prod_p)
    return z_hp, p_hp, k_hp


def lp2bp_zpk_mp(z, p, k, wo=1.0, bw=1.0):
    wo = mpmath.mpf(wo)
    bw = mpmath.mpf(bw)
    degree = len(p) - len(z)

    z_lp = [x * bw / 2 for x in z]
    p_lp = [x * bw / 2 for x in p]

    z_bp = []
    for x in z_lp:
        sqrt_val = mpmath.sqrt(x**2 - wo**2)
        z_bp.append(x + sqrt_val)
        z_bp.append(x - sqrt_val)

    p_bp = []
    for x in p_lp:
        sqrt_val = mpmath.sqrt(x**2 - wo**2)
        p_bp.append(x + sqrt_val)
        p_bp.append(x - sqrt_val)

    z_bp = z_bp + [mpmath.mpf(0.0)] * degree
    k_bp = k * (bw**degree)
    return z_bp, p_bp, k_bp


def lp2bs_zpk_mp(z, p, k, wo=1.0, bw=1.0):
    wo = mpmath.mpf(wo)
    bw = mpmath.mpf(bw)
    degree = len(p) - len(z)

    z_hp = [(bw / 2) / x for x in z]
    p_hp = [(bw / 2) / x for x in p]

    z_bs = []
    for x in z_hp:
        sqrt_val = mpmath.sqrt(x**2 - wo**2)
        z_bs.append(x + sqrt_val)
        z_bs.append(x - sqrt_val)

    p_bs = []
    for x in p_hp:
        sqrt_val = mpmath.sqrt(x**2 - wo**2)
        p_bs.append(x + sqrt_val)
        p_bs.append(x - sqrt_val)

    z_bs = z_bs + [1j * wo] * degree + [-1j * wo] * degree

    prod_z = product_mp([-x for x in z])
    prod_p = product_mp([-x for x in p])
    k_bs = k * mpmath.re(prod_z / prod_p)
    return z_bs, p_bs, k_bs


def bilinear_zpk_mp(z, p, k, fs):
    fs = mpmath.mpf(fs)
    degree = len(p) - len(z)
    fs2 = 2.0 * fs

    z_z = [(fs2 + x) / (fs2 - x) for x in z]
    p_z = [(fs2 + x) / (fs2 - x) for x in p]
    z_z = z_z + [mpmath.mpf(-1.0)] * degree

    prod_z = product_mp([fs2 - x for x in z])
    prod_p = product_mp([fs2 - x for x in p])
    k_z = k * mpmath.re(prod_z / prod_p)
    return z_z, p_z, k_z


def butter_zpk_mp(N, Wn, btype="low", analog=False, fs=None):
    btype_map = {
        "low": "lowpass",
        "lowpass": "lowpass",
        "high": "highpass",
        "highpass": "highpass",
        "band": "bandpass",
        "bandpass": "bandpass",
        "stop": "bandstop",
        "bandstop": "bandstop",
    }
    btype = btype_map.get(btype.lower())
    if btype is None:
        raise ValueError(f"'{btype}' is an invalid bandtype for filter.")

    if fs is not None:
        if analog:
            raise ValueError("fs cannot be specified for an analog filter")
        fs_mp = mpmath.mpf(fs)
    else:
        fs_mp = None

    if isinstance(Wn, str):
        Wn_mp = [mpmath.mpf(Wn)]
    elif isinstance(Wn, collections.abc.Iterable):
        Wn_mp = [mpmath.mpf(x) for x in Wn]
    else:
        Wn_mp = [mpmath.mpf(Wn)]

    if fs_mp is not None:
        Wn_mp = [w / (fs_mp / 2) for w in Wn_mp]

    if not analog:
        for w in Wn_mp:
            if w <= 0 or w >= 1:
                raise ValueError(
                    "Digital filter critical frequencies must be 0 < Wn < 1"
                )
        fs_internal = mpmath.mpf(2.0)
        warped = [
            2 * fs_internal * mpmath.tan(mpmath.pi * w / fs_internal) for w in Wn_mp
        ]
    else:
        for w in Wn_mp:
            if w <= 0:
                raise ValueError(
                    "Analog filter critical frequencies must be greater than 0"
                )
        warped = Wn_mp

    z, p, k = buttap_mp(N)

    if btype in ("lowpass", "highpass"):
        if len(warped) != 1:
            raise ValueError(
                "Must specify a single critical frequency Wn for lowpass or highpass filter"
            )
        wo = warped[0]
        if btype == "lowpass":
            z, p, k = lp2lp_zpk_mp(z, p, k, wo)
        else:
            z, p, k = lp2hp_zpk_mp(z, p, k, wo)
    elif btype in ("bandpass", "bandstop"):
        if len(warped) != 2:
            raise ValueError(
                "Wn must specify start and stop frequencies for bandpass or bandstop filter"
            )
        if warped[0] >= warped[1]:
            raise ValueError("Wn[0] must be less than Wn[1]")
        bw = warped[1] - warped[0]
        wo = mpmath.sqrt(warped[0] * warped[1])
        if btype == "bandpass":
            z, p, k = lp2bp_zpk_mp(z, p, k, wo, bw)
        else:
            z, p, k = lp2bs_zpk_mp(z, p, k, wo, bw)

    if not analog:
        z, p, k = bilinear_zpk_mp(z, p, k, fs=2.0)

    return z, p, k


def ellip_zpk_mp(N, rp, rs, Wn, btype="low", analog=False, fs=None):
    """Arbitrary precision Elliptic filter design returning zeros, poles, and gain."""
    btype_map = {
        "low": "lowpass",
        "lowpass": "lowpass",
        "high": "highpass",
        "highpass": "highpass",
        "band": "bandpass",
        "bandpass": "bandpass",
        "stop": "bandstop",
        "bandstop": "bandstop",
    }
    btype = btype_map.get(btype.lower())
    if btype is None:
        raise ValueError(f"'{btype}' is an invalid bandtype for filter.")

    if fs is not None:
        if analog:
            raise ValueError("fs cannot be specified for an analog filter")
        fs_mp = mpmath.mpf(fs)
    else:
        fs_mp = None

    if isinstance(Wn, str):
        Wn_mp = [mpmath.mpf(Wn)]
    elif isinstance(Wn, collections.abc.Iterable):
        Wn_mp = [mpmath.mpf(x) for x in Wn]
    else:
        Wn_mp = [mpmath.mpf(Wn)]

    if fs_mp is not None:
        Wn_mp = [w / (fs_mp / 2) for w in Wn_mp]

    if not analog:
        for w in Wn_mp:
            if w <= 0 or w >= 1:
                raise ValueError(
                    "Digital filter critical frequencies must be 0 < Wn < 1"
                )
        fs_internal = mpmath.mpf(2.0)
        warped = [
            2 * fs_internal * mpmath.tan(mpmath.pi * w / fs_internal) for w in Wn_mp
        ]
    else:
        for w in Wn_mp:
            if w <= 0:
                raise ValueError(
                    "Analog filter critical frequencies must be greater than 0"
                )
        warped = Wn_mp

    # Retrieve analog prototype
    z, p, k = ellipap_mp(N, rp, rs)

    if btype in ("lowpass", "highpass"):
        if len(warped) != 1:
            raise ValueError(
                "Must specify a single critical frequency Wn for lowpass or highpass filter"
            )
        wo = warped[0]
        if btype == "lowpass":
            z, p, k = lp2lp_zpk_mp(z, p, k, wo)
        else:
            z, p, k = lp2hp_zpk_mp(z, p, k, wo)
    elif btype in ("bandpass", "bandstop"):
        if len(warped) != 2:
            raise ValueError(
                "Wn must specify start and stop frequencies for bandpass or bandstop filter"
            )
        if warped[0] >= warped[1]:
            raise ValueError("Wn[0] must be less than Wn[1]")
        bw = warped[1] - warped[0]
        wo = mpmath.sqrt(warped[0] * warped[1])
        if btype == "bandpass":
            z, p, k = lp2bp_zpk_mp(z, p, k, wo, bw)
        else:
            z, p, k = lp2bs_zpk_mp(z, p, k, wo, bw)

    if not analog:
        z, p, k = bilinear_zpk_mp(z, p, k, fs=2.0)

    return z, p, k


def _cplxreal_mp(z, tol=None):
    """Split into complex and real parts, combining conjugate pairs in mpmath.

    The 1-D input vector `z` is split up into its complex (`zc`) and real (`zr`)
    elements. Every complex element must be part of a complex-conjugate pair.
    Two complex numbers are considered a conjugate pair if their real and imaginary
    parts differ in magnitude by less than ``tol * abs(z)``.
    """
    if len(z) == 0:
        return [], []

    if tol is None:
        # Scale tolerance with the configured mpmath precision
        tol = mpmath.mpf(100) * (mpmath.mpf(10) ** (-mpmath.mp.dps * 4 // 5))

    # Sort by real part, and then magnitude of imaginary part
    z = sorted(z, key=lambda x: (mpmath.re(x), abs(mpmath.im(x))))

    zr = []
    zc_candidates = []
    for x in z:
        if abs(mpmath.im(x)) <= tol * abs(x):
            zr.append(mpmath.re(x))
        else:
            zc_candidates.append(x)

    if len(zc_candidates) == 0:
        return [], zr

    zp = []
    zn = []
    for x in zc_candidates:
        if mpmath.im(x) > 0:
            zp.append(x)
        else:
            zn.append(x)

    if len(zp) != len(zn):
        raise ValueError("Array contains complex value with no matching conjugate.")

    # Find runs of (approximately) the same real part
    runs = []
    start = 0
    for i in range(len(zp) - 1):
        if not (abs(mpmath.re(zp[i + 1]) - mpmath.re(zp[i])) <= tol * abs(zp[i])):
            runs.append((start, i + 1))
            start = i + 1
    runs.append((start, len(zp)))

    # Sort each run by the magnitude of its imaginary part
    for start, stop in runs:
        zp_slice = zp[start:stop]
        zp[start:stop] = sorted(zp_slice, key=lambda x: abs(mpmath.im(x)))

        zn_slice = zn[start:stop]
        zn[start:stop] = sorted(zn_slice, key=lambda x: abs(mpmath.im(x)))

    # Check that negatives match positives
    for i in range(len(zp)):
        if abs(zp[i] - mpmath.conj(zn[i])) > tol * abs(zn[i]):
            raise ValueError("Array contains complex value with no matching conjugate.")

    # Average out numerical inaccuracies in the real/imag parts
    zc = []
    for i in range(len(zp)):
        zc.append((zp[i] + mpmath.conj(zn[i])) / 2)

    return zc, zr


def _nearest_real_complex_idx_mp(fro, to, which, tol=None):
    """Get the next closest real or complex element based on distance."""
    if tol is None:
        tol = mpmath.mpf(100) * (mpmath.mpf(10) ** (-mpmath.mp.dps))

    candidates = []
    for idx, val in enumerate(fro):
        dist = abs(val - to)
        candidates.append((idx, val, dist))

    # Sort by distance
    candidates_sorted = sorted(candidates, key=lambda x: x[2])

    for idx, val, dist in candidates_sorted:
        is_real = abs(mpmath.im(val)) <= tol * abs(val)
        if which == "any":
            return idx
        elif which == "real" and is_real:
            return idx
        elif which == "complex" and not is_real:
            return idx

    raise ValueError(f"No matching {which} element found.")


def poly_mp(roots):
    """Find the coefficients of a polynomial with the given roots."""
    c = [mpmath.mpf(1.0)]
    for r in roots:
        next_c = [mpmath.mpf(0.0)] * (len(c) + 1)
        for i in range(len(c)):
            next_c[i] += c[i]
            next_c[i + 1] -= c[i] * r
        c = next_c
    return [mpmath.re(x) for x in c]


def zpk2tf_mp(z, p, k):
    """Return polynomial transfer function representation from zeros, poles, and gain."""
    b = [x * k for x in poly_mp(z)]
    a = poly_mp(p)
    return b, a


def _single_zpksos_mp(z, p, k):
    """Create one second-order section from up to two zeros and poles."""
    b, a = zpk2tf_mp(z, p, k)
    sos = [mpmath.mpf(0.0)] * 6
    for i in range(len(b)):
        sos[3 - len(b) + i] = b[i]
    for i in range(len(a)):
        sos[6 - len(a) + i] = a[i]
    return sos


def zpk2sos_mp(z, p, k, pairing=None, *, analog=False):
    """Return second-order sections from zeros, poles, and gain in mpmath."""
    z = [mpmath.mpc(x) for x in z]
    p = [mpmath.mpc(x) for x in p]
    k_val = mpmath.mpc(k)

    tol = mpmath.mpf(100) * (mpmath.mpf(10) ** (-mpmath.mp.dps * 4 // 5))

    if abs(mpmath.im(k_val)) > tol * abs(k_val):
        raise ValueError("k must be real")
    k = mpmath.re(k_val)

    if pairing is None:
        pairing = "minimal" if analog else "nearest"

    valid_pairings = ["nearest", "keep_odd", "minimal"]
    if pairing not in valid_pairings:
        raise ValueError(f"pairing must be one of {valid_pairings}, not {pairing}")

    if analog and pairing != "minimal":
        raise ValueError('for analog zpk2sos conversion, pairing must be "minimal"')

    if len(z) == len(p) == 0:
        if not analog:
            return [
                [
                    k,
                    mpmath.mpf(0.0),
                    mpmath.mpf(0.0),
                    mpmath.mpf(1.0),
                    mpmath.mpf(0.0),
                    mpmath.mpf(0.0),
                ]
            ]
        else:
            return [
                [
                    mpmath.mpf(0.0),
                    mpmath.mpf(0.0),
                    k,
                    mpmath.mpf(0.0),
                    mpmath.mpf(0.0),
                    mpmath.mpf(1.0),
                ]
            ]

    if pairing != "minimal":
        p_len, z_len = len(p), len(z)
        if z_len > p_len:
            p = p + [mpmath.mpc(0.0)] * (z_len - p_len)
        elif p_len > z_len:
            z = z + [mpmath.mpc(0.0)] * (p_len - z_len)

        n_sections = (max(len(p), len(z)) + 1) // 2

        if len(p) % 2 == 1 and pairing == "nearest":
            p = p + [mpmath.mpc(0.0)]
            z = z + [mpmath.mpc(0.0)]
        assert len(p) == len(z)
    else:
        if len(p) < len(z):
            raise ValueError("for analog zpk2sos conversion, must have len(p)>=len(z)")
        n_sections = (len(p) + 1) // 2

    zc_z, zr_z = _cplxreal_mp(z, tol=tol)
    z = zc_z + zr_z

    zc_p, zr_p = _cplxreal_mp(p, tol=tol)
    p = zc_p + zr_p

    if not analog:

        def idx_worst(p_list):
            min_val = None
            min_idx = -1
            for idx, val in enumerate(p_list):
                score = abs(mpmath.mpf(1.0) - abs(val))
                if min_val is None or score < min_val:
                    min_val = score
                    min_idx = idx
            return min_idx

    else:

        def idx_worst(p_list):
            min_val = None
            min_idx = -1
            for idx, val in enumerate(p_list):
                score = abs(mpmath.re(val))
                if min_val is None or score < min_val:
                    min_val = score
                    min_idx = idx
            return min_idx

    sos = [[mpmath.mpf(0.0)] * 6 for _ in range(n_sections)]

    for si in range(n_sections - 1, -1, -1):
        p1_idx = idx_worst(p)
        p1 = p[p1_idx]
        del p[p1_idx]

        def is_real_mp(val):
            return abs(mpmath.im(val)) <= tol * abs(val)

        def sum_real_mp(lst):
            return sum(1 for x in lst if is_real_mp(x))

        if is_real_mp(p1) and sum_real_mp(p) == 0:
            if pairing != "minimal":
                z1_idx = _nearest_real_complex_idx_mp(z, p1, "real", tol=tol)
                z1 = z[z1_idx]
                del z[z1_idx]
                sos[si] = _single_zpksos_mp(
                    [z1, mpmath.mpf(0.0)],
                    [p1, mpmath.mpf(0.0)],
                    mpmath.mpf(1.0),
                )
            elif len(z) > 0:
                z1_idx = _nearest_real_complex_idx_mp(z, p1, "real", tol=tol)
                z1 = z[z1_idx]
                del z[z1_idx]
                sos[si] = _single_zpksos_mp([z1], [p1], mpmath.mpf(1.0))
            else:
                sos[si] = _single_zpksos_mp([], [p1], mpmath.mpf(1.0))

        elif (
            len(p) + 1 == len(z)
            and not is_real_mp(p1)
            and sum_real_mp(p) == 1
            and sum_real_mp(z) == 1
        ):
            z1_idx = _nearest_real_complex_idx_mp(z, p1, "complex", tol=tol)
            z1 = z[z1_idx]
            del z[z1_idx]
            sos[si] = _single_zpksos_mp(
                [z1, mpmath.conj(z1)], [p1, mpmath.conj(p1)], mpmath.mpf(1.0)
            )

        else:
            if is_real_mp(p1):
                real_p_indices = [idx for idx, val in enumerate(p) if is_real_mp(val)]
                real_p_values = [p[idx] for idx in real_p_indices]
                worst_real_p_sub_idx = idx_worst(real_p_values)
                p2_idx = real_p_indices[worst_real_p_sub_idx]
                p2 = p[p2_idx]
                del p[p2_idx]
            else:
                p2 = mpmath.conj(p1)

            if len(z) > 0:
                z1_idx = _nearest_real_complex_idx_mp(z, p1, "any", tol=tol)
                z1 = z[z1_idx]
                del z[z1_idx]

                if not is_real_mp(z1):
                    sos[si] = _single_zpksos_mp(
                        [z1, mpmath.conj(z1)], [p1, p2], mpmath.mpf(1.0)
                    )
                else:
                    if len(z) > 0:
                        z2_idx = _nearest_real_complex_idx_mp(z, p1, "real", tol=tol)
                        z2 = z[z2_idx]
                        assert is_real_mp(z2)
                        del z[z2_idx]
                        sos[si] = _single_zpksos_mp([z1, z2], [p1, p2], mpmath.mpf(1.0))
                    else:
                        sos[si] = _single_zpksos_mp([z1], [p1, p2], mpmath.mpf(1.0))
            else:
                sos[si] = _single_zpksos_mp([], [p1, p2], mpmath.mpf(1.0) or z)

    assert len(p) == len(z) == 0

    sos[0][0] *= k
    sos[0][1] *= k
    sos[0][2] *= k

    sos_final = []
    for row in sos:
        sos_final.append([mpmath.re(x) for x in row])
    return sos_final


def normalize_mp(b, a, tol=None):
    """Normalize numerator/denominator of a continuous-time transfer function in mpmath."""
    if tol is None:
        tol = mpmath.mpf(100) * (mpmath.mpf(10) ** (-mpmath.mp.dps))
    b = [mpmath.mpc(x) for x in b]
    a = [mpmath.mpc(x) for x in a]

    # Trim leading zeros from denominator
    first_a = 0
    while first_a < len(a) and abs(a[first_a]) <= tol:
        first_a += 1
    if first_a == len(a):
        raise ValueError("Denominator must have at least one non-zero element.")
    a = a[first_a:]

    # Trim leading zeros from numerator
    first_b = 0
    while first_b < len(b) and abs(b[first_b]) <= tol:
        first_b += 1
    if first_b == len(b):
        b = [mpmath.mpc(0.0)]
    else:
        b = b[first_b:]

    a0 = a[0]
    b = [x / a0 for x in b]
    a = [x / a0 for x in a]
    return b, a


def polyroots_mp_eig(p, tol=None):
    """
    Polynomial root finder using companion matrix eigenvalues followed by
    Newton-Raphson refinement using Horner's method. Replaces mpmath.polyroots
    for ill-conditioned or widely-scaled polynomials.
    """
    if tol is None:
        # Scale tolerance to account for polynomial conditioning limits
        tol = mpmath.mpf(100) * (mpmath.mpf(10) ** (-mpmath.mp.dps * 4 // 5))

    p = [mpmath.mpc(x) for x in p]

    # Strip leading zeros
    non_zero = []
    for idx, val in enumerate(p):
        if abs(val) > tol:
            non_zero.append(idx)

    if not non_zero:
        return []

    trailing_zeros = len(p) - non_zero[-1] - 1
    p_stripped = p[non_zero[0] : non_zero[-1] + 1]

    N = len(p_stripped)
    if N > 1:
        n_matrix = N - 1
        A = mpmath.matrix(n_matrix, n_matrix)
        p0 = p_stripped[0]
        for j in range(n_matrix):
            A[0, j] = -p_stripped[j + 1] / p0
        for i in range(1, n_matrix):
            A[i, i - 1] = mpmath.mpf(1.0)

        # Robust eigenvalue extraction
        E, L = mpmath.eig(A)
        roots = [x for x in E]

        # Polish roots to full precision using Newton-Raphson
        refined = []
        for r in roots:
            curr = mpmath.mpc(r)
            for _ in range(5):
                val = mpmath.mpc(0.0)
                deriv = mpmath.mpc(0.0)
                for c in p_stripped:
                    deriv = deriv * curr + val
                    val = val * curr + c
                if abs(deriv) == 0:
                    break
                diff = val / deriv
                curr -= diff
                if abs(diff) < 10 ** (-mpmath.mp.dps - 5):
                    break
            refined.append(curr)
        roots = refined
    else:
        roots = []

    roots = roots + [mpmath.mpc(0.0)] * trailing_zeros
    return roots


def tf2zpk_mp(b, a, tol=None):
    """Return zero, pole, gain representation from a transfer function in mpmath."""
    b, a = normalize_mp(b, a, tol=tol)
    k = b[0]
    b_norm = [x / k for x in b]

    # z = mpmath.polyroots(b_norm, maxsteps=256) if len(b_norm) > 1 else []
    # p = mpmath.polyroots(a, maxsteps=256) if len(a) > 1 else []
    z = polyroots_mp_eig(b_norm) if len(b_norm) > 1 else []
    p = polyroots_mp_eig(a) if len(a) > 1 else []

    z = [mpmath.mpc(x) for x in z]
    p = [mpmath.mpc(x) for x in p]
    return z, p, k


def tf2sos_mp(b, a, pairing=None, *, analog=False):
    """Return second-order sections from transfer function coefficients in mpmath."""
    z, p, k = tf2zpk_mp(b, a)
    return zpk2sos_mp(z, p, k, pairing=pairing, analog=analog)


# -- Arbitrary-precision mpmath implementations


def lfilter_mp(b, a, x):
    """
    Direct Form I 1-D filter implementation in mpmath.
    Assumes standard filter difference equation:
    a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ... - a[1]*y[n-1] - ...
    """
    N_x = len(x)
    y = [mpmath.mpf(0.0)] * N_x
    Nb = len(b)
    Na = len(a)

    for n in range(N_x):
        val_b = mpmath.mpf(0.0)
        for i in range(Nb):
            if n - i >= 0:
                val_b += b[i] * x[n - i]

        val_a = mpmath.mpf(0.0)
        for j in range(1, Na):
            if n - j >= 0:
                val_a += a[j] * y[n - j]

        y[n] = (val_b - val_a) / a[0]
    return y


def sosfilt_mp(sos, x):
    """
    Filters a signal x through a cascade of SOS rows in mpmath.
    Each row is represented as [b0, b1, b2, a0, a1, a2].
    """
    curr_x = list(x)
    for section in sos:
        b = section[0:3]
        a = section[3:6]
        curr_x = lfilter_mp(b, a, curr_x)
    return curr_x


def filter_polyphase_ba_mp(q_polyphase, a, inputs):
    M = len(q_polyphase)
    N_samples = len(inputs[0])
    sum_fir = [mpmath.mpf(0.0)] * N_samples
    for k in range(M):
        filtered_k = lfilter_mp(q_polyphase[k], [mpmath.mpf(1.0)], inputs[k])
        for n in range(N_samples):
            sum_fir[n] += filtered_k[n]
    return lfilter_mp([mpmath.mpf(1.0)], a, sum_fir)


def filter_polyphase_sos_mp(sos_polyphase, inputs):
    M = len(sos_polyphase)
    N_samples = len(inputs[0])
    y = [mpmath.mpf(0.0)] * N_samples
    for k in range(M):
        filtered_k = sosfilt_mp(sos_polyphase[k], inputs[k])
        for n in range(N_samples):
            y[n] += filtered_k[n]
    return y


def filter_polyphase_hybrid_mp(q_polyphase, sos_sections, inputs):
    M = len(q_polyphase)
    N_samples = len(inputs[0])
    sum_fir = [mpmath.mpf(0.0)] * N_samples
    for k in range(M):
        filtered_k = lfilter_mp(q_polyphase[k], [mpmath.mpf(1.0)], inputs[k])
        for n in range(N_samples):
            sum_fir[n] += filtered_k[n]

    sos = []
    for sec in sos_sections:
        sos.append(
            [
                mpmath.mpf(1.0),
                mpmath.mpf(0.0),
                mpmath.mpf(0.0),
                mpmath.mpf(1.0),
                sec[0],
                sec[1],
            ]
        )

    return sosfilt_mp(sos, sum_fir)


# --- Tests ---


def test_butter_zpk_mp():
    import numpy as np
    from scipy import signal

    def sort_complex_list(lst):
        """Sort complex values consistently for reliable comparisons."""
        return sorted(
            lst, key=lambda val: (float(mpmath.re(val)), float(mpmath.im(val)))
        )

    # Test cases representing various configurations
    test_cases = [
        # (N, Wn, btype, analog, fs)
        (3, 0.4, "low", False, None),
        (4, 0.3, "high", False, None),
        (2, [0.2, 0.5], "band", False, None),
        (3, [0.1, 0.6], "stop", False, None),
        (5, 100, "low", True, None),
        (4, [10, 50], "band", True, None),
        (3, 1000, "high", False, 44100),
        (2, [50, 200], "stop", False, 1000),
    ]

    print(
        f"{'Case':<5} | {'N':<2} | {'btype':<8} | {'analog':<6} | {'fs':<7} | {'Status'}"
    )
    print("-" * 55)

    for i, (N, Wn, btype, analog, fs) in enumerate(test_cases):
        z_sc, p_sc, k_sc = signal.butter(
            N, Wn, btype=btype, analog=analog, output="zpk", fs=fs
        )

        z_mp, p_mp, k_mp = butter_zpk_mp(N, Wn, btype=btype, analog=analog, fs=fs)

        z_sc_sorted = np.sort(np.round(z_sc, 12))
        z_mp_sorted = np.sort(
            np.round([complex(x) for x in sort_complex_list(z_mp)], 12)
        )

        p_sc_sorted = np.sort(np.round(p_sc, 12))
        p_mp_sorted = np.sort(
            np.round([complex(x) for x in sort_complex_list(p_mp)], 12)
        )

        assert np.allclose(
            z_sc_sorted, z_mp_sorted, atol=1e-12
        ), f"Zeros mismatch in case {i+1}"
        assert np.allclose(
            p_sc_sorted, p_mp_sorted, atol=1e-12
        ), f"Poles mismatch in case {i+1}"
        assert np.allclose(
            k_sc, float(k_mp), atol=1e-12
        ), f"Gain mismatch in case {i+1}"

        print(
            f"#{i+1:<4} | {N:<2} | {btype:<8} | {str(analog):<6} | {str(fs):<7} | Passed"
        )


def test_ellip_zpk_mp():
    import numpy as np
    import scipy.signal as signal

    def sort_complex_list(lst):
        """Sort complex values consistently for reliable comparisons."""
        return sorted(
            lst, key=lambda val: (float(mpmath.re(val)), float(mpmath.im(val)))
        )

    def comp_arrays(a, b, atol=1e-12):
        """Sort and compare two arrays of complex numbers."""
        if len(a) != len(b):
            return False
        if len(a) == 0:
            return True
        return np.allclose(np.sort(a), np.sort(b), atol=atol)

    # Set high precision for mpmath operations
    mpmath.mp.dps = 50

    test_cases = [
        # (N, rp, rs, Wn, btype, analog, fs)
        (3, 1, 40, 0.4, "low", False, None),
        (4, 0.5, 60, 0.3, "high", False, None),
        (2, 1, 50, [0.2, 0.5], "band", False, None),
        (3, 0.1, 80, [0.1, 0.6], "stop", False, None),
        (5, 1, 40, 100, "low", True, None),
        (4, 0.5, 50, [10, 50], "band", True, None),
        (3, 1, 60, 1000, "high", False, 44100),
        (2, 0.5, 40, [50, 200], "stop", False, 1000),
        # Edge cases
        (0, 1, 40, 0.4, "low", False, None),
        (1, 1, 40, 0.4, "low", False, None),
        (1, 1, 40, 0.4, "high", False, None),
        (1, 1, 40, [0.2, 0.5], "band", False, None),
    ]

    print(
        f"{'No.':<4} | {'N':<2} | {'rp':<4} | {'rs':<4} | {'btype':<6} | {'analog':<6} | {'fs':<7} | {'Status'}"
    )
    print("-" * 65)

    for i, (N, rp, rs, Wn, btype, analog, fs) in enumerate(test_cases):
        # Compute using SciPy's float64 implementation
        z_sp, p_sp, k_sp = signal.ellip(
            N, rp, rs, Wn, btype=btype, analog=analog, output="zpk", fs=fs
        )

        # Compute using mpmath arbitrary precision implementation
        z_mp, p_mp, k_mp = ellip_zpk_mp(
            N, rp, rs, Wn, btype=btype, analog=analog, fs=fs
        )

        # Convert high-precision numbers back to python floats/complex for assertion checks
        z_mp_c = [complex(x) for x in sort_complex_list(z_mp)]
        p_mp_c = [complex(x) for x in sort_complex_list(p_mp)]

        # Run comparative assertions
        assert comp_arrays(z_sp, z_mp_c), f"Zeros mismatch in case {i+1}"
        assert comp_arrays(p_sp, p_mp_c), f"Poles mismatch in case {i+1}"
        assert np.allclose(
            k_sp, float(k_mp), atol=1e-12
        ), f"Gain mismatch in case {i+1}"

        print(
            f"#{i+1:<3} | {N:<2} | {rp:<4} | {rs:<4} | {btype:<6} | {str(analog):<6} | {str(fs):<7} | Passed"
        )


def test_tf2sos():
    """
    This also test `zpk2sos` as `tf2sos` is a thin wrapper of `zpk2sos`.
    """

    import numpy as np
    import scipy.signal as signal

    # Set mpmath precision
    mpmath.mp.dps = 50

    # -------------------------------------------------------------
    # Test 1: Digital Filter ZPK conversion (using butter_zpk_mp)
    # -------------------------------------------------------------
    # Design a 6th order Butterworth bandpass digital filter
    z_bp, p_bp, k_bp = butter_zpk_mp(6, [0.2, 0.5], btype="bandpass")

    # Convert using our mpmath implementation
    sos_mp = zpk2sos_mp(z_bp, p_bp, k_bp)

    # Convert using SciPy's float64 implementation
    z_sp, p_sp, k_sp = signal.butter(6, [0.2, 0.5], btype="bandpass", output="zpk")
    sos_sp = signal.zpk2sos(z_sp, p_sp, k_sp)

    # Compare values
    diffs = []
    for r_mp, r_sp in zip(sos_mp, sos_sp):
        r_mp_float = [float(x) for x in r_mp]
        diff = max(abs(a - b) for a, b in zip(r_mp_float, r_sp))
        diffs.append(diff)

    print(
        "Max absolute difference for butter(6, [0.2, 0.5]):",
        max(diffs),
    )
    # Expected Output: Max absolute difference is on the order of 1e-16 (machine epsilon)

    # -------------------------------------------------------------
    # Test 2: Continuous-Time (Analog) TF conversion (tf2sos)
    # -------------------------------------------------------------
    b_ana = [1.0]
    a_ana = [1.0, 2.0, 2.0, 1.0]  # s^3 + 2s^2 + 2s + 1

    sos_ana_sp = signal.tf2sos(b_ana, a_ana, analog=True)
    sos_ana_mp = tf2sos_mp(b_ana, a_ana, analog=True)

    print("\nSciPy Analog SOS:")
    print(sos_ana_sp)

    print("\nmpmath Analog SOS (cast to float):")
    print(np.array([[float(y) for y in x] for x in sos_ana_mp]))


def test_filters():
    import numpy as np
    from scipy.signal import lfilter, sosfilt
    from design import design_polyphase_butterworth, apply

    def filter_polyphase_ba_float(q_polyphase, a_float, inputs):
        """
        Simulates the standard 'ba' representation.
        Each branch is filtered with an FIR filter Q_k(z), the outputs are
        accumulated, and the sum is passed through the low-rate IIR denominator 1/A(z).
        """
        M = len(q_polyphase)
        sum_fir = np.zeros_like(inputs[0], dtype=float)
        for k in range(M):
            sum_fir += lfilter(q_polyphase[k], [1.0], inputs[k])
        return lfilter([1.0], a_float, sum_fir)

    def filter_polyphase_sos_float(sos_polyphase, inputs):
        """
        Simulates the 'sos' representation.
        Each branch filters its input with a distinct complete SOS cascade (Q_k(z)/A(z)),
        and the branch outputs are summed.
        """
        M = len(sos_polyphase)
        y = np.zeros_like(inputs[0], dtype=float)
        for k in range(M):
            y += sosfilt(sos_polyphase[k], inputs[k])
        return y

    def filter_polyphase_hybrid_float(q_polyphase, sos_sections, inputs):
        """
        Simulates the 'hybrid' representation.
        Each branch filters its input with an FIR Q_k(z) and the outputs are summed.
        The sum is passed through the shared all-pole IIR denominator implemented
        as cascaded second-order sections (1 / (1 + a1*z^-1 + a2*z^-2)).
        """
        M = len(q_polyphase)
        sum_fir = np.zeros_like(inputs[0], dtype=float)
        for k in range(M):
            sum_fir += lfilter(q_polyphase[k], [1.0], inputs[k])

        # Formulate Scipy-compatible SOS matrix:
        # Each section is [b0, b1, b2, a0, a1, a2]
        # Here, numerator is [1.0, 0.0, 0.0], denominator is [1.0, a1, a2]
        sos = []
        for sec in sos_sections:
            sos.append([1.0, 0.0, 0.0, 1.0, sec[0], sec[1]])

        return sosfilt(sos, sum_fir)

    # order, cutoff, M = 16, 0.015625, 8  # Target spec, slow to compute.
    order, cutoff, M = 8, 0.2, 4  # Lightweight parameters for testing.
    N_samples = 100

    inputs_mp = [
        [2 * mpmath.rand() - mpmath.mpf(1) for _ in range(N_samples)] for _ in range(M)
    ]

    q_ba_mp, a_ba_mp = design_polyphase_butterworth(
        order, cutoff, M, "ba", as_float=False, workdps=2000
    )
    sos_poly_mp = design_polyphase_butterworth(
        order, cutoff, M, "sos", as_float=False, workdps=2000
    )
    q_hyb_mp, sos_hyb_mp = design_polyphase_butterworth(
        order, cutoff, M, "hybrid", as_float=False, workdps=2000
    )

    mpmath.mp.dps = 100
    y_ba_mp = filter_polyphase_ba_mp(q_ba_mp, a_ba_mp, inputs_mp)
    y_sos_mp = filter_polyphase_sos_mp(sos_poly_mp, inputs_mp)
    y_hyb_mp = filter_polyphase_hybrid_mp(q_hyb_mp, sos_hyb_mp, inputs_mp)

    def to_f64(v):
        return apply(v, float)

    inputs_f64 = to_f64(inputs_mp)
    q_ba_f64, a_ba_f64 = to_f64(q_ba_mp), to_f64(a_ba_mp)
    sos_poly_f64 = to_f64(sos_poly_mp)
    q_hyb_f64, sos_hyb_f64 = to_f64(q_hyb_mp), to_f64(sos_hyb_mp)

    y_ba_flt = filter_polyphase_ba_float(q_ba_f64, a_ba_f64, inputs_f64)
    y_sos_flt = filter_polyphase_sos_float(sos_poly_f64, inputs_f64)
    y_hyb_flt = filter_polyphase_hybrid_float(q_hyb_f64, sos_hyb_f64, inputs_f64)

    def diff_np(x, y):
        return np.max(np.abs(x - y))

    print("--- 1. Testing Float Implementations (Double Precision) ---")
    print("Max absolute difference (BA vs SOS):   ", diff_np(y_ba_flt, y_sos_flt))
    print("Max absolute difference (BA vs Hybrid):", diff_np(y_ba_flt, y_hyb_flt))
    print("Max absolute difference (SOS vs Hybrid):", diff_np(y_sos_flt, y_hyb_flt))

    def diff_mp(x, y):
        return mpmath.nstr(max(abs(a - b) for a, b in zip(x, y)))

    print("\n--- 2. Testing Arbitrary-Precision mpmath Implementations (100 dps) ---")
    print("Max absolute difference (BA vs SOS):   ", diff_mp(y_ba_mp, y_sos_mp))
    print("Max absolute difference (BA vs Hybrid):", diff_mp(y_ba_mp, y_hyb_mp))
    print("Max absolute difference (SOS vs Hybrid):", diff_mp(y_sos_mp, y_hyb_mp))


if __name__ == "__main__":
    test_butter_zpk_mp()
    test_ellip_zpk_mp()
    test_tf2sos()
    test_filters()
