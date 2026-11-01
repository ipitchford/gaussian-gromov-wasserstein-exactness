#!/usr/bin/env python3
"""Certified bracket for squared-distance Gaussian Gromov-Wasserstein.

Standard library only. Given the covariance eigenvalues of two centred
Gaussians, ``bracket(a, b)`` returns the Gaussian-restricted optimum Q_G of
eq. (4) and a certified interval for the unrestricted optimum Q_*:

    Q_G - min(8(A2 B2 - a.b), Q_G/2)  <=  Q_*  <=  Q_G - lower_gap,

where lower_gap is the largest of the proved certificates S/(A2 B2),
32 S^2 / K_C and max_ij 32 Delta_ij^2 / K_ij (manuscript Section 6). The upper
gap 8(A2 B2 - a.b) is the comparison of Delon, Desolneux and Salmona (2022);
the cap Q_G/2 follows from Q_G = 4(alpha-beta)^2 + 8||a-b||^2 >= 16(A2 B2 - a.b).

Floating-point arithmetic is used; the certificates are proved inequalities
evaluated in double precision, not interval arithmetic.
"""
from __future__ import annotations

import math


def _prepare(a, b):
    a = sorted((float(x) for x in a), reverse=True)
    b = sorted((float(x) for x in b), reverse=True)
    if any(x < 0 or not math.isfinite(x) for x in a + b):
        raise ValueError("eigenvalues must be finite and nonnegative")
    d = max(len(a), len(b))
    a += [0.0] * (d - len(a))
    b += [0.0] * (d - len(b))
    return a, b


def bracket(a, b, pair_certificate=True):
    """Return Q_G, a certified interval for Q_*, cos(theta) and relative gaps.

    ``a`` and ``b`` are covariance eigenvalue lists (any order, any lengths;
    they are sorted nonincreasingly and zero-padded). Both must have positive
    rank. The pairwise certificate costs O(d^2); disable it for very large d.
    """
    a, b = _prepare(a, b)
    d = len(a)
    alpha, beta = math.fsum(a), math.fsum(b)
    A2sq = math.fsum(x * x for x in a)
    B2sq = math.fsum(x * x for x in b)
    if A2sq == 0 or B2sq == 0:
        raise ValueError("both spectra need positive rank (Theorem 8)")
    A2, B2 = math.sqrt(A2sq), math.sqrt(B2sq)
    ab = math.fsum(x * y for x, y in zip(a, b))
    q_gauss = 4 * ((alpha - beta) ** 2 + 2 * math.fsum((x - y) ** 2 for x, y in zip(a, b)))
    # S = sum_{i<j} C_ij^2 with C = b a^T - a b^T; summed directly (no cancellation).
    S = 0.0
    curv_cov = 0.0       # 8 sum_{i<j} C_ij^2 (a_i b_i + a_j b_j)
    spread_a = 0.0       # sum_{i<j} C_ij^2 (a_i - a_j)^2
    spread_b = 0.0
    best_pair = 0.0
    for i in range(d):
        ai, bi = a[i], b[i]
        for j in range(i + 1, d):
            aj, bj = a[j], b[j]
            c = bi * aj - ai * bj
            if c == 0.0:
                continue
            c2 = c * c
            S += c2
            curv_cov += c2 * (ai * bi + aj * bj)
            spread_a += c2 * (ai - aj) ** 2
            spread_b += c2 * (bi - bj) ** 2
            if pair_certificate:
                pa = 24 * (ai * ai + aj * aj) - 32 * ai * aj
                pb = 24 * (bi * bi + bj * bj) - 32 * bi * bj
                k_ij = math.sqrt(max(pa, 0.0) * max(pb, 0.0)) + 8 * (ai * bi + aj * bj)
                if k_ij > 0:
                    best_pair = max(best_pair, 32 * c2 / k_ij)
    # ||C a||^2 = A2^2 S and ||C b||^2 = B2^2 S (Lagrange identity), eq. (16).
    p_a = 8 * A2sq * S + 16 * spread_a
    p_b = 8 * B2sq * S + 16 * spread_b
    k_c = math.sqrt(p_a * p_b) + 8 * curv_cov
    certificates = {"S_over_A2B2": S / (A2 * B2),
                    "sharp_32S2_over_KC": 32 * S * S / k_c if k_c > 0 else 0.0}
    if pair_certificate:
        certificates["best_pair_32D2_over_Kij"] = best_pair
    lower_gap = max(certificates.values())
    upper_dds = 8 * (A2 * B2 - ab)
    upper_gap = min(upper_dds, q_gauss / 2)
    cos_theta = ab / (A2 * B2)
    rel = (lambda x: x / q_gauss) if q_gauss > 0 else (lambda x: 0.0)
    return {
        "Q_G": q_gauss,
        "Q_star_interval": (q_gauss - upper_gap, q_gauss - lower_gap),
        "lower_gap": lower_gap,
        "upper_gap": upper_gap,
        "upper_gap_source": "8(A2B2-a.b)" if upper_dds <= q_gauss / 2 else "Q_G/2",
        "relative_lower_gap": rel(lower_gap),
        "relative_upper_gap": rel(upper_gap),
        "cos_theta": cos_theta,
        "S": S,
        "K_C": k_c,
        "certificates": certificates,
        "certified_time_4S_over_KC": 4 * S / k_c if k_c > 0 else 0.0,
    }


if __name__ == "__main__":
    import json
    result = bracket([1, 1], [1])
    # Plane-to-line (manuscript Section 7): Q_G = 12 and the sharp certificate
    # 32/(8+8*sqrt(6)) = 4/(1+sqrt(6)) reproduce the paper's value.
    if not (result["Q_G"] == 12.0
            and abs(result["lower_gap"] - 4 / (1 + math.sqrt(6))) < 1e-12
            and abs(result["upper_gap"] - 8 * (math.sqrt(2) - 1)) < 1e-12):
        raise SystemExit("CHECK FAILED: plane-to-line demo does not reproduce manuscript values")
    print("Plane-to-line, a=(1,1), b=(1):")
    print(json.dumps(result, indent=2))
