#!/usr/bin/env python3
"""Extension checks for the squared-distance Gaussian GW paper.

Sections (keys of extension_checks.json):

  a_factor_two      Q_G >= 16(A2B2 - a.b), i.e. the known upper gap is at most
                    Q_G/2 (exact rational check on an integer grid, and a
                    floating-point scan of random spectra).
  b_comonotone      Comonotone (Hoeffding-Frechet) refinement of the upper gap,
                    D <= 4(M_co - alpha beta - 2 a.b). Plane-to-line value by
                    mpmath (two independent integral representations);
                    general weighted chi-square cases by numerical quantiles.
  c_plane_to_line   Bracket table for N(0,I_2) versus N(0,1).
  d_equal_dim       Worked example a=(2,1), b=(1,1).
  e_realistic       Trace-normalised realistic spectra table.
  f_path_curvature  Whole-path check |G''(t)| <= K_C and the Taylor lower
                    bound along the W flow (methodology review W3).
  g_near_isotropic  Sharpness of the factor-two corollary near isotropy:
                    (32 S^2/K_C)/Q_G versus d/(2(d+2)), and upper ratio <= 1/2.

Labels: "certified" means a proved inequality evaluated in floating point
(or exactly where stated); "numerical" means quadrature/ODE evaluation
without interval certification. Every check raises SystemExit explicitly,
so the script is also correct under ``python -O``.
"""
from __future__ import annotations

import argparse
import gc
import importlib.util
import itertools
import json
import math
import platform
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath
import numpy as np
import scipy
from scipy import optimize, special, stats
from scipy.integrate import quad, solve_ivp
from scipy.special import roots_hermitenorm

HERE = Path(__file__).resolve().parent
CHECK_COUNT = 0


def check(condition, message):
    """Explicit check; raises SystemExit so that -O cannot remove it."""
    global CHECK_COUNT
    if not bool(condition):
        raise SystemExit(f"CHECK FAILED: {message}")
    CHECK_COUNT += 1


def load_bracket():
    # Loaded by path so the script also works under ``python -I``; no bytecode
    # cache is written next to the shipped sources.
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("gw_gaussian_bracket",
                                                  HERE / "gw_gaussian_bracket.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# Spectral quantities (independent NumPy implementation)
# ---------------------------------------------------------------------------

def pad(a, b):
    a = np.sort(np.asarray(a, float))[::-1]
    b = np.sort(np.asarray(b, float))[::-1]
    d = max(len(a), len(b))
    return np.pad(a, (0, d - len(a))), np.pad(b, (0, d - len(b)))


def spectral(a, b, pairs=True):
    a, b = pad(a, b)
    alpha, beta = a.sum(), b.sum()
    A2, B2 = np.linalg.norm(a), np.linalg.norm(b)
    ab = float(a @ b)
    QG = 4 * ((alpha - beta) ** 2 + 2 * np.sum((a - b) ** 2))
    C = np.outer(b, a) - np.outer(a, b)
    iu = np.triu_indices(len(a), 1)
    C2 = C[iu] ** 2
    S = float(C2.sum())

    def P(x):
        Cx = C @ x
        return 8 * Cx @ Cx + 16 * np.sum(C2 * (x[iu[0]] - x[iu[1]]) ** 2)

    KC = math.sqrt(P(a) * P(b)) + 8 * np.sum(C2 * (a[iu[0]] * b[iu[0]] + a[iu[1]] * b[iu[1]]))
    out = dict(a=a, b=b, alpha=alpha, beta=beta, A2=A2, B2=B2, ab=ab, QG=QG, S=S, KC=KC,
               C=C, cos=ab / (A2 * B2), low=S / (A2 * B2),
               sharp=32 * S * S / KC if KC > 0 else 0.0,
               up_dds=8 * (A2 * B2 - ab))
    if pairs:
        i, j = iu
        Dl = b[i] * a[j] - b[j] * a[i]
        Pf = lambda x, y: 24 * (x * x + y * y) - 32 * x * y
        K = np.sqrt(Pf(a[i], a[j]) * Pf(b[i], b[j])) + 8 * (a[i] * b[i] + a[j] * b[j])
        with np.errstate(divide="ignore", invalid="ignore"):
            cert = np.where((Dl != 0) & (K > 0), 32 * Dl * Dl / K, 0.0)
        out["pair"] = float(cert.max()) if len(cert) else 0.0
    return out


# ---------------------------------------------------------------------------
# (a) Factor-two corollary
# ---------------------------------------------------------------------------

def part_a(named_cases):
    # Exact: Q_G >= 16(A2B2 - a.b)  <=>  (Q_G/16 + a.b)^2 >= A2^2 B2^2 since both
    # sides of Q_G/16 + a.b >= A2B2 are nonnegative. Same for ||a-b||^2 >= 2(A2B2-a.b).
    exact_cases = 0
    triples = [t for t in itertools.product(range(5), repeat=3) if t[0] >= t[1] >= t[2]]
    for a in triples:
        for b in triples:
            if not any(a) or not any(b):
                continue
            A2sq = sum(x * x for x in a); B2sq = sum(x * x for x in b)
            ab = sum(x * y for x, y in zip(a, b))
            QG = 4 * ((sum(a) - sum(b)) ** 2 + 2 * sum((x - y) ** 2 for x, y in zip(a, b)))
            nab = sum((x - y) ** 2 for x, y in zip(a, b))
            check(QG == 4 * (sum(a) - sum(b)) ** 2 + 8 * nab, f"Q_G identity {a},{b}")
            check((Fraction(QG, 16) + ab) ** 2 >= A2sq * B2sq, f"exact factor-two {a},{b}")
            check((Fraction(nab, 2) + ab) ** 2 >= A2sq * B2sq, f"exact ||a-b||^2 bound {a},{b}")
            exact_cases += 1
    rng = np.random.default_rng(20261031)
    worst = 0.0
    worst_case = None
    n_random = 0
    for trial in range(20000):
        d = int(rng.choice([1, 2, 3, 5, 10, 50]))
        da, db = int(rng.integers(1, d + 1)), int(rng.integers(1, d + 1))
        kind = trial % 3
        if kind == 0:
            a, b = rng.random(da), rng.random(db)
        elif kind == 1:
            a, b = rng.exponential(size=da) ** 3, rng.exponential(size=db) ** 3
        else:   # near-proportional / near-isotropic
            base = rng.random(d) + 0.5
            a = base[:da] * (1 + 0.05 * rng.random(da))
            b = (rng.random() + 0.2) * base[:db] * (1 + 0.05 * rng.random(db))
        s = spectral(a, b, pairs=False)
        a_, b_ = s["a"], s["b"]
        ident = 4 * (s["alpha"] - s["beta"]) ** 2 + 8 * np.sum((a_ - b_) ** 2)
        scale = max(1.0, s["QG"], s["A2"] * s["B2"])
        check(abs(ident - s["QG"]) <= 1e-12 * scale, f"Q_G identity random trial {trial}")
        check(np.sum((a_ - b_) ** 2) >= 2 * (s["A2"] * s["B2"] - s["ab"]) - 1e-12 * scale,
              f"||a-b||^2 >= 2(A2B2-a.b) random trial {trial}")
        check(s["QG"] >= 16 * (s["A2"] * s["B2"] - s["ab"]) - 1e-12 * scale,
              f"Q_G >= 16(A2B2-a.b) random trial {trial}")
        if s["QG"] > 1e-9 * scale:
            r = s["up_dds"] / s["QG"]
            if r > worst:
                worst, worst_case = r, (d, da, db, kind)
        n_random += 1
    check(worst <= 0.5 + 1e-12, f"upper gap / Q_G exceeds 1/2: {worst}")
    named = {}
    for name, (a, b) in named_cases.items():
        s = spectral(a, b, pairs=False)
        check(s["QG"] >= 16 * (s["A2"] * s["B2"] - s["ab"]) * (1 - 1e-12), f"named {name}")
        named[name] = {"upper_dds_over_QG": s["up_dds"] / s["QG"]}
    return {"statement": "Q_G = 4(alpha-beta)^2 + 8||a-b||^2 >= 16(A2B2 - a.b), hence "
                         "8(A2B2 - a.b) <= Q_G/2 and Q_G/2 <= Q_* <= Q_G",
            "exact_integer_pairs_d3_entries_0_to_4": exact_cases,
            "random_trials": n_random, "seed": 20261031,
            "max_upper_gap_over_QG_random": worst, "worst_random_case(d,da,db,kind)": worst_case,
            "named_cases": named, "status": "passed"}


# ---------------------------------------------------------------------------
# Quantiles of weighted chi-square laws R = sum_i lam_i Z_i^2
# ---------------------------------------------------------------------------

class EqualWeights:
    """lam * chi^2_k: exact SciPy quantiles."""
    def __init__(self, k, lam):
        self.k, self.lam = int(k), float(lam)
        self.mean, self.var = self.k * self.lam, 2 * self.k * self.lam ** 2
        self.lam_max = self.lam
        self.method = f"exact scipy chi2({self.k}) quantiles"

    def quantile_s(self, s):        # quantile at u = 1 - exp(-s)
        s = np.asarray(s, float)
        upper = stats.chi2.isf(np.exp(-s), self.k)
        lower = stats.chi2.ppf(-np.expm1(-s), self.k)
        return self.lam * np.where(s > math.log(2), upper, lower)


class _Bisect:
    def quantile_s(self, s):
        s = np.asarray(s, float)
        tsf = np.exp(-s); tcdf = -np.expm1(-s)
        use_sf = s > math.log(2)
        lo = np.zeros_like(s)
        hi = self.mean + 2 * self.lam_max * (s + 8) + 10 * math.sqrt(self.var)
        for _ in range(60):
            bad = self.sf(hi) > tsf
            if not bad.any():
                break
            hi = np.where(bad, 2 * hi, hi)
        check(not (self.sf(hi) > tsf).any(), "quantile bracket expansion")
        for _ in range(64):
            mid = 0.5 * (lo + hi)
            below = np.where(use_sf, self.sf(mid) > tsf, self.cdf(mid) < tcdf)
            lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
        return 0.5 * (lo + hi)


class TwoWeights(_Bisect):
    """lam1 Z1^2 + lam2 Z2^2 via z1 = sqrt(x/lam1) sin t (smooth 1-D integrals)."""
    def __init__(self, lam1, lam2, nodes=256):
        self.l1, self.l2 = float(lam1), float(lam2)
        t, w = np.polynomial.legendre.leggauss(nodes)
        self.t, self.w = (t + 1) * math.pi / 4, w * math.pi / 4
        self.mean = self.l1 + self.l2
        self.var = 2 * (self.l1 ** 2 + self.l2 ** 2)
        self.lam_max = max(self.l1, self.l2)
        self.method = f"Gauss-Legendre ({nodes}) after z1=sqrt(x/l1) sin t; bisection"

    def _parts(self, x):
        x = np.asarray(x, float)[..., None]
        c = np.sqrt(x / self.l1)
        jac = np.exp(-0.5 * (c * np.sin(self.t)) ** 2) / math.sqrt(2 * math.pi)
        jac = 2 * jac * c * np.cos(self.t) * self.w
        arg = np.sqrt(x / self.l2) * np.cos(self.t)
        return c[..., 0], jac, arg

    def cdf(self, x):
        c, jac, arg = self._parts(x)
        return np.sum(jac * (2 * special.ndtr(arg) - 1), axis=-1)

    def sf(self, x):
        c, jac, arg = self._parts(x)
        return 2 * special.ndtr(-c) + np.sum(jac * 2 * special.ndtr(-arg), axis=-1)


class Imhof(_Bisect):
    """Imhof (1961) inversion on a fixed composite Gauss-Legendre grid.

    P(R > x) = 1/2 + (1/pi) int_0^inf sin(phi(u) - x u/2) / (u rho(u)) du,
    phi = (1/2) sum arctan(lam u), rho = prod (1 + lam^2 u^2)^(1/4).
    The grid is truncated where the envelope tail int_U^inf du/(u rho) < tail_tol.
    Suitable when many eigenvalues are nonnegligible (fast envelope decay).
    """
    def __init__(self, lam, tail_tol=1e-12, x_max=None, nodes_per_panel=16):
        lam = np.asarray(lam, float)
        lam = lam[lam > 0]
        self.lam = lam
        self.mean, self.var = lam.sum(), 2 * np.sum(lam ** 2)
        self.lam_max = lam.max()
        if x_max is None:
            x_max = self.mean + 2 * self.lam_max * 40 + 12 * math.sqrt(self.var)
        self.x_max = x_max

        def log_env(u):
            return -math.log(u) - 0.25 * np.sum(np.log1p((lam * u) ** 2))

        def tail(U):
            val, _ = quad(lambda u: math.exp(log_env(u)), U, np.inf, limit=500)
            return val
        U = 10 / self.lam_max
        while tail(U) > tail_tol:
            U *= 1.5
            check(U < 1e7 / self.lam_max, "Imhof truncation too large; use another method")
        self.U, self.tail_bound = U, tail(U)
        width = min(math.pi / (x_max / 2), 0.5 / self.lam_max)
        n_panels = int(math.ceil(U / width))
        check(n_panels < 200000, "Imhof grid too fine")
        g, gw = np.polynomial.legendre.leggauss(nodes_per_panel)
        edges = np.linspace(0, U, n_panels + 1)
        mid, half = (edges[1:] + edges[:-1]) / 2, (edges[1:] - edges[:-1]) / 2
        u = (mid[:, None] + half[:, None] * g).ravel()
        w = (half[:, None] * gw).ravel()
        phi = np.zeros_like(u); lrho = np.zeros_like(u)
        for chunk in np.array_split(np.arange(len(lam)), max(1, len(lam) // 64)):
            lu = np.outer(u, lam[chunk])
            phi += 0.5 * np.arctan(lu).sum(axis=1)
            lrho += 0.25 * np.log1p(lu ** 2).sum(axis=1)
        self.u, self.phi, self.amp = u, phi, w * np.exp(-lrho) / u / math.pi
        self.nodes = len(u)
        self.method = (f"Imhof inversion, {self.nodes} Gauss-Legendre nodes on [0,{U:.4g}], "
                       f"envelope tail bound {self.tail_bound:.1e}; bisection")

    def _integral(self, x):
        x = np.atleast_1d(np.asarray(x, float))
        out = np.empty_like(x)
        for k in range(0, len(x), 64):
            xs = x[k:k + 64]
            out[k:k + 64] = np.sin(self.phi[None, :] - 0.5 * xs[:, None] * self.u[None, :]) @ self.amp
        return out

    def sf(self, x):
        return 0.5 + self._integral(x)

    def cdf(self, x):
        return 0.5 - self._integral(x)


def chi_law(lam):
    lam = np.asarray(lam, float)
    lam = np.sort(lam[lam > 0])[::-1]
    if np.allclose(lam, lam[0], rtol=1e-15, atol=0):
        return EqualWeights(len(lam), lam[0])
    if len(lam) == 2:
        return TwoWeights(lam[0], lam[1])
    return Imhof(lam)


def s_grid(s_max, k=16):
    edges = [0, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.2, 0.5, 1, 2, 3.5, 5, 7, 9.5, 12.5, 16, 20, 25, 30, 35, 40, 45]
    edges = [e for e in edges if e < s_max] + [s_max]
    g, gw = np.polynomial.legendre.leggauss(k)
    s, w = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s.append((hi + lo) / 2 + (hi - lo) / 2 * g); w.append((hi - lo) / 2 * gw)
    return np.concatenate(s), np.concatenate(w)


def comonotone(lawR, lawS, s_max):
    """M_co - alpha beta = int_0^1 (F_R^{-1}(u) - alpha)(F_S^{-1}(u) - beta) du.

    Substitution u = 1 - exp(-s); composite Gauss-Legendre in s (16 and 24
    nodes per panel; the difference is reported as an error estimate);
    truncation at s_max with a linear-in-s tail extrapolation using the
    asymptotic quantile slope 2 lam_max.
    """
    q1 = lawR.quantile_s([s_max])[0] - lawR.mean
    q2 = lawS.quantile_s([s_max])[0] - lawS.mean
    k1, k2 = 2 * lawR.lam_max, 2 * lawS.lam_max
    tail = math.exp(-s_max) * (q1 * q2 + q1 * k2 + q2 * k1 + 2 * k1 * k2)
    values = {}
    for k in (16, 24):
        s, w = s_grid(s_max, k)
        qR, qS = lawR.quantile_s(s), lawS.quantile_s(s)
        e = np.exp(-s) * w
        cR, cS = qR - lawR.mean, qS - lawS.mean
        values[k] = float(np.sum(cR * cS * e)) + tail
        if k == 24:
            diag = {"int_QR_minus_mean": float(np.sum(cR * e)), "int_QS_minus_mean": float(np.sum(cS * e)),
                    "int_QR2_minus_second_moment": float(np.sum(qR ** 2 * e)) - (lawR.var + lawR.mean ** 2),
                    "int_QS2_minus_second_moment": float(np.sum(qS ** 2 * e)) - (lawS.var + lawS.mean ** 2),
                    "s_grid_16_vs_24_difference": abs(values[24] - values[16])}
    return values[24], {"s_max": s_max, "tail_estimate": tail, "methods": [lawR.method, lawS.method],
                        "quadrature_diagnostics": diag}


def mpmath_plane_to_line(dps=30):
    mpmath.mp.dps = dps
    # Form 1: int_0^1 (-2 ln(1-u)) * (Phi^{-1}((1+u)/2))^2 du, Phi^{-1}((1+u)/2) = sqrt2 erfinv(u)
    f1 = lambda u: -2 * mpmath.log(1 - u) * 2 * mpmath.erfinv(u) ** 2
    m1 = mpmath.quad(f1, [0, 0.5, 0.9, 0.99, 0.999, 0.9999, 1])
    # Form 2: u = 2 Phi(y) - 1:  int_0^inf y^2 (-2 ln erfc(y/sqrt2)) 2 phi(y) dy
    f2 = lambda y: y ** 2 * (-2 * mpmath.log(mpmath.erfc(y / mpmath.sqrt(2)))) * 2 * mpmath.npdf(y)
    m2 = mpmath.quad(f2, [0, 1, 2, 4, 8, 16, mpmath.inf])
    return m1, m2


def part_b():
    m1, m2 = mpmath_plane_to_line()
    check(abs(m1 - m2) < mpmath.mpf("1e-20"), f"mpmath forms disagree: {m1} vs {m2}")
    M = m1
    # Plane-to-line: alpha=2, beta=1, a.b=1, C = 36, Q_G = 12.
    lower_co = 36 - 4 * (M + 2)
    lower_cs = 12 - 8 * (mpmath.sqrt(2) - 1)
    check(lower_co > lower_cs, "comonotone bound must improve Cauchy-Schwarz")
    check(M <= 2 + mpmath.sqrt(8), "comonotone M exceeds Cauchy-Schwarz value")
    # Validate the numerical machinery against the mpmath value.
    exact_laws, _ = comonotone(EqualWeights(2, 1.0), EqualWeights(1, 1.0), 45)
    two_law, _ = comonotone(TwoWeights(1.0, 1.0), EqualWeights(1, 1.0), 40)
    val_err_exact = abs(exact_laws + 2 - float(M))
    val_err_two = abs(two_law + 2 - float(M))
    check(val_err_exact < 1e-9, f"s-quadrature validation error {val_err_exact}")
    check(val_err_two < 1e-7, f"TwoWeights quantile validation error {val_err_two}")
    # Imhof machinery against exact chi2_8 (lam = 1) and a scaled chi2_20.
    imhof_val = {}
    for k, lam in [(8, 1.0), (20, 0.05)]:
        law = Imhof(np.full(k, lam))
        s = np.array([0.01, 0.3, 1.0, 3.0, 8.0, 15.0])
        err = np.max(np.abs(law.quantile_s(s) - EqualWeights(k, lam).quantile_s(s))) / lam
        check(err < 1e-6, f"Imhof quantile validation chi2_{k}: {err}")
        imhof_val[f"chi2_{k}_x{lam}"] = float(err)
    general = {}
    # 1 versus n isotropic: a=1_n, b=(1). alpha=n, beta=1, a.b=1, A2B2=sqrt(n).
    for n in range(2, 7):
        mc, info = comonotone(EqualWeights(n, 1.0), EqualWeights(1, 1.0), 45)
        Mco = mc + n
        QG = 4 * ((n - 1) ** 2 + 2 * (n - 1))
        lb_co = QG - 4 * (Mco - n - 2)
        lb_cs = QG - 8 * (math.sqrt(n) - 1)
        check(lb_co >= lb_cs - 1e-9, f"comonotone vs CS, n={n}")
        general[f"line_vs_R{n}"] = {"a": [1.0] * n, "b": [1.0], "M_co": Mco, "Q_G": QG,
                                    "lower_Q_comonotone": lb_co, "lower_Q_cauchy_schwarz": lb_cs,
                                    "improvement": lb_co - lb_cs, "method": info["methods"]}
    for name, a, b in [("a=(2,1),b=(1,1)", [2., 1.], [1., 1.]),
                       ("a=(3,1),b=(1)", [3., 1.], [1.]),
                       ("a=(1,0.5),b=(1,1)", [1., .5], [1., 1.]),
                       ("a=(1,1,1,1),b=(2,1)", [1.] * 4, [2., 1.])]:
        sp = spectral(a, b, pairs=False)
        mc, info = comonotone(chi_law(a), chi_law(b), 40)
        up_co = 4 * (mc - 2 * sp["ab"])
        check(mc <= 2 * sp["A2"] * sp["B2"] * (1 + 1e-9), f"comonotone <= CS for {name}")
        general[name] = {"a": a, "b": b, "M_co_minus_alpha_beta": mc, "Q_G": sp["QG"],
                         "upper_gap_comonotone": up_co, "upper_gap_dds": sp["up_dds"],
                         "lower_Q_comonotone": sp["QG"] - up_co,
                         "lower_Q_cauchy_schwarz": sp["QG"] - sp["up_dds"],
                         "method": info["methods"], "diagnostics": info["quadrature_diagnostics"]}
    return {"statement": "Q_* >= C - 4(M_co + 2 a.b), i.e. D <= 4(M_co - alpha beta - 2 a.b), "
                         "M_co = int_0^1 F^{-1}_{|X|^2}(u) F^{-1}_{|Y|^2}(u) du",
            "plane_to_line": {"M_co_form1_mpmath": mpmath.nstr(m1, 20),
                              "M_co_form2_mpmath": mpmath.nstr(m2, 20),
                              "mpmath_dps": 30,
                              "cauchy_schwarz_value_2_plus_sqrt8": mpmath.nstr(2 + mpmath.sqrt(8), 20),
                              "lower_bound_Q_star_comonotone": mpmath.nstr(lower_co, 16),
                              "lower_bound_Q_star_cauchy_schwarz": mpmath.nstr(lower_cs, 16),
                              "improvement": mpmath.nstr(lower_co - lower_cs, 12),
                              "label": "certified inequality; integral evaluated by mpmath "
                                       "tanh-sinh quadrature in two independent forms "
                                       "(agreement < 1e-20), not interval arithmetic"},
            "machinery_validation": {"exact_quantiles_vs_mpmath": val_err_exact,
                                     "two_weight_quantiles_vs_mpmath": val_err_two,
                                     "imhof_vs_exact_chi2_max_quantile_error": imhof_val},
            "general_cases": general,
            "general_label": "numerical evaluation (quantile bisection + Gauss-Legendre in s), "
                             "not interval-certified",
            "M_co_plane_to_line": float(M), "status": "passed"}


# ---------------------------------------------------------------------------
# Gauss-Hermite W-flow utilities
# ---------------------------------------------------------------------------

def gh_tensor(order, d):
    """Tensor probabilists' Gauss-Hermite rule for N(0, I_d) (as in verify_gw.py)."""
    x, w = roots_hermitenorm(order)
    w = w / math.sqrt(2 * math.pi)
    grids = np.meshgrid(*([x] * d), indexing="ij")
    Z = np.stack([g.ravel() for g in grids], axis=1)
    W = np.prod(np.stack(np.meshgrid(*([w] * d), indexing="ij")), axis=0).ravel()
    return Z, W, f"tensor Gauss-Hermite order {order}, full space"


def orthant_rule(d, L, width, k):
    """Composite Gauss-Legendre rule on [0, L]^d for N(0, I_d), weights x 2^d.

    The W flow is sign-equivariant (W_k is odd in z_k and even in z_j, j != k),
    so every integrand used below (H_a H_b(T), Z_i T_i, generator terms) is even
    in each coordinate; off-diagonal correlations E[Z_i T_j] vanish exactly by
    the same symmetry. Mass outside the box is below 2 d Phi(-L).
    """
    g, gw = np.polynomial.legendre.leggauss(k)
    edges = np.linspace(0, L, int(round(L / width)) + 1)
    x = np.concatenate([(lo + hi) / 2 + (hi - lo) / 2 * g for lo, hi in zip(edges[:-1], edges[1:])])
    w = np.concatenate([(hi - lo) / 2 * gw for lo, hi in zip(edges[:-1], edges[1:])])
    w = w * 2 * np.exp(-x * x / 2) / math.sqrt(2 * math.pi)
    grids = np.meshgrid(*([x] * d), indexing="ij")
    Z = np.stack([gg.ravel() for gg in grids], axis=1)
    W = np.prod(np.stack(np.meshgrid(*([w] * d), indexing="ij")), axis=0).ravel()
    return Z, W, f"orthant composite Gauss-Legendre [0,{L}]^{d}, panel {width}, {k} nodes/panel"


def flow_path(a, b, rule, times, callback, rtol=1e-11, atol=1e-12):
    """Integrate dz/dt = z * (C (z^2 - 1)) from the nodes of ``rule``.

    ``times`` must contain 0; the positive and negative branches are integrated
    outward segment by segment (to bound memory) and ``callback(t, T_t(Z))`` is
    called at every requested time, with T of shape (d, n).
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = len(a)
    C = np.outer(b, a) - np.outer(a, b)
    Z = rule[0]
    n = len(Z)

    def rhs(_, y):
        z = y.reshape(d, n)
        return (z * (C @ (z * z - 1))).ravel()
    times = np.unique(np.asarray(times, float))
    check(np.any(times == 0), "flow_path needs t=0")
    callback(0.0, Z.T.copy())
    for branch in (times[times > 0], times[times < 0][::-1]):
        y, t_prev = Z.T.ravel().copy(), 0.0
        for t in branch:
            sol = solve_ivp(rhs, (t_prev, float(t)), y, method="DOP853", t_eval=[float(t)],
                            rtol=rtol, atol=atol)   # t_eval: keep only the end state (memory)
            check(sol.success, f"ODE failure: {sol.message}")
            y, t_prev = sol.y[:, -1].copy(), float(t)
            del sol
            gc.collect()   # solver objects sit in reference cycles; free their arrays now
            callback(float(t), y.reshape(d, n))
    return C


def G_value(a, b, Z, W, T):
    """G = E|X|^2|Y|^2 + 2||K||_F^2 with X = D_a^{1/2} Z, Y = D_b^{1/2} T(Z).

    Off-diagonal E[Z_i T_j] vanish by sign symmetry, so ||K||^2 = sum a_i b_i r_i^2.
    """
    Ha = (Z * Z) @ a
    Hb = b @ (T * T)
    r = (Z.T * T) @ W
    return W @ (Ha * Hb) + 2 * np.sum(a * b * r * r)


def G_second_derivative_generator(a, b, C, Z, W, T):
    """Analytic G''(t) using the generator: d/dt f(T_t) = (Lf)(T_t)."""
    eta = T * T - 1
    Ce = C @ eta                                  # (C eta)_m at T
    Ha = (Z * Z) @ a
    bz2 = b[:, None] * T * T
    L2Hb = 4 * np.sum(T * T * Ce * (b[:, None] * Ce + C.T @ bz2), axis=0)
    m2 = W @ (Ha * L2Hb)
    Wv = T * Ce                                   # W(T)
    LW = T * (Ce * Ce + 2 * (C @ (T * T * Ce)))   # (W . grad W_i)(T)
    Zt = Z.T
    r = (Zt * T) @ W                              # r_i = E[Z_i T_i]
    r1 = (Zt * Wv) @ W
    r2 = (Zt * LW) @ W
    return m2 + 2 * np.sum(a * b * (2 * r1 ** 2 + 2 * r * r2))


def Q_from_G(a, b, G):
    alpha, beta = a.sum(), b.sum()
    C0 = 4 * (alpha ** 2 + 2 * a @ a) + 4 * (beta ** 2 + 2 * b @ b) - 4 * alpha * beta
    return C0 - 4 * G


# ---------------------------------------------------------------------------
# (c) Plane-to-line bracket
# ---------------------------------------------------------------------------

class Rearrangement:
    """Y = Phi^{-1}(F_h(h(X))), h(x) = x1 + sign(x1)(lam x2^2 + mu x2^4 + kap x1^2).

    Q evaluated by deterministic polar Gauss-Legendre quadrature on the half
    plane x1 > 0 (the integrands are even under x1 -> -x1 since Y is odd), and
    F_h by Gauss-Legendre in x2 with the x1 integral in closed form.
    """
    def __init__(self, r_nodes=24, th_nodes=96, x2_nodes=64):
        edges = [0, 1e-3, 1e-2, 0.05, 0.2, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8.5, 10]
        g, gw = np.polynomial.legendre.leggauss(r_nodes)
        r = np.concatenate([(h + l) / 2 + (h - l) / 2 * g for l, h in zip(edges[:-1], edges[1:])])
        wr = np.concatenate([(h - l) / 2 * gw for l, h in zip(edges[:-1], edges[1:])])
        t, wt = np.polynomial.legendre.leggauss(th_nodes)
        th, wth = t * math.pi / 2, wt * math.pi / 2
        R, TH = np.meshgrid(r, th, indexing="ij")
        self.x1, self.x2 = (R * np.cos(TH)).ravel(), (R * np.sin(TH)).ravel()
        self.w = (2 * np.outer(wr * r * np.exp(-r * r / 2), wth) / (2 * math.pi)).ravel()
        self.g, self.gw = np.polynomial.legendre.leggauss(x2_nodes)
        self.label = f"polar GL r:{len(r)} x theta:{th_nodes}, F_h GL {x2_nodes}"

    def sf_h(self, s, lam, mu, kap):
        """P(h(X) > s) for s > 0."""
        if lam == 0 and mu == 0:
            ystar = np.full_like(s, np.inf)
        else:
            ystar = 2 * s / (lam + np.sqrt(lam * lam + 4 * mu * s))
        L = np.minimum(np.sqrt(ystar), 10.0)
        x2 = (L[:, None] / 2) * (self.g[None, :] + 1)
        c = s[:, None] - lam * x2 ** 2 - mu * x2 ** 4
        c = np.maximum(c, 0.0)
        gx = 2 * c / (1 + np.sqrt(1 + 4 * kap * c))
        inner = np.sum((L[:, None] / 2) * self.gw * np.exp(-x2 ** 2 / 2) / math.sqrt(2 * math.pi)
                       * special.ndtr(-gx), axis=1)
        outside = special.ndtr(-np.sqrt(ystar))      # P(|x2| > sqrt(y*), x1 > 0); 0 if y* = inf
        return outside + 2 * inner

    def evaluate(self, lam, mu, kap):
        x1, x2, w = self.x1, self.x2, self.w
        h = x1 + lam * x2 ** 2 + mu * x2 ** 4 + kap * x1 ** 2
        Y = -special.ndtri(self.sf_h(h, lam, mu, kap))
        R2 = x1 ** 2 + x2 ** 2
        K1, K2 = w @ (x1 * Y), w @ (x2 * Y)   # K2 is exactly zero by symmetry
        G = w @ (R2 * Y * Y) + 2 * (K1 ** 2 + K2 ** 2)
        return {"Q": 36 - 4 * G, "E_Y2_minus_1": w @ (Y * Y) - 1, "E_Y4_minus_3": w @ Y ** 4 - 3,
                "E_X1Y": K1, "E_X2Y": K2}


def part_c(M_co):
    lower_dds = 12 - 8 * (math.sqrt(2) - 1)
    lower_co = 28 - 4 * M_co
    a, b = np.array([1., 1.]), np.array([1., 0.])
    # Numerical W flow (= pair flow here, C_12 = 1). Composite Gauss-Legendre on
    # the orthant is the headline rule; tensor Gauss-Hermite (as in verify_gw.py)
    # is reported for reference: it under-resolves the fast winding of large
    # Lotka-Volterra orbits once t > 0.15.
    rules = {"GL_fine": orthant_rule(2, 8.0, 0.125, 12), "GL_coarse": orthant_rule(2, 8.0, 0.125, 8),
             "GH_256": gh_tensor(256, 2), "GH_384": gh_tensor(384, 2)}
    grid = np.round(np.arange(0, 0.6001, 0.01), 10)
    flow = {}
    for key, rule in rules.items():
        Z, W, label = rule
        qrec = {}      # scalars only, to bound memory

        def keep(t, T, qrec=qrec, Z=Z, W=W):
            qrec[t] = (Q_from_G(a, b, G_value(a, b, Z, W, T)), W @ T[0] ** 2 - 1, W @ T[0] ** 4 - 3)
        flow_path(a, b, rule, np.append(grid, 1 / 7), keep)
        Q = np.array([qrec[float(t)][0] for t in grid])
        k = int(np.argmin(Q))
        q17 = qrec[float(1 / 7)]
        entry = {"rule": label, "nodes": len(W), "Q_grid": Q.tolist(), "t_grid_min": float(grid[k]),
                 "Q_grid_min": float(Q[k]),
                 "Q_at_1_over_7": float(q17[0]),
                 "marginal_Y2_error_at_1_over_7": float(q17[1]),
                 "marginal_Y4_error_at_1_over_7": float(q17[2])}
        if key.startswith("GL"):
            tl, th = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]

            def q_at(t, rule=rule):
                box = {}
                flow_path(a, b, rule, [0.0, t], lambda tt, T: box.__setitem__(tt, T))
                return Q_from_G(a, b, G_value(a, b, rule[0], rule[1], box[float(t)]))
            res = optimize.minimize_scalar(q_at, bounds=(tl, th), method="bounded",
                                           options={"xatol": 1e-7})
            box = {}
            flow_path(a, b, rule, [0.0, float(res.x)],
                      lambda tt, T: box.__setitem__(tt, (W @ T[0] ** 2 - 1, W @ T[0] ** 4 - 3)))
            m2, m4 = box[float(res.x)]
            entry.update({"t_min": float(res.x), "Q_min": float(res.fun),
                          "marginal_Y2_error_at_t_min": float(m2),
                          "marginal_Y4_error_at_t_min": float(m4)})
        flow[key] = entry
    fine_c, coarse_c = flow["GL_fine"], flow["GL_coarse"]
    check(abs(fine_c["Q_min"] - coarse_c["Q_min"]) < 1e-7, "flow minimum quadrature convergence (GL)")
    check(abs(fine_c["marginal_Y2_error_at_t_min"]) < 1e-9 and abs(fine_c["marginal_Y4_error_at_t_min"]) < 1e-8,
          "flow marginal residuals at t_min")
    check(abs(fine_c["Q_at_1_over_7"] - flow["GH_384"]["Q_at_1_over_7"]) < 5e-4,
          "GL and Gauss-Hermite agree at t=1/7")
    check(fine_c["Q_at_1_over_7"] < 76 / 7, "numerical flow at 1/7 below certified 76/7")
    check(fine_c["Q_min"] < 12 - 4 / (1 + math.sqrt(6)), "flow minimum below optimised certificate")
    # Monotone rearrangement couplings.
    rea = Rearrangement()
    lin = rea.evaluate(0.0, 0.0, 0.0)
    check(abs(lin["Q"] - 12) < 1e-8, f"linear control Q={lin['Q']}")
    fun2 = lambda p: rea.evaluate(abs(p[0]), 0.0, abs(p[1]))["Q"]
    r2 = optimize.minimize(fun2, [0.37, 0.17], method="Nelder-Mead",
                           options={"xatol": 1e-6, "fatol": 1e-10, "maxiter": 400})
    fun3 = lambda p: rea.evaluate(abs(p[0]), abs(p[2]), abs(p[1]))["Q"]
    r3 = optimize.minimize(fun3, [abs(r2.x[0]), abs(r2.x[1]), 0.0], method="Nelder-Mead",
                           options={"xatol": 1e-6, "fatol": 1e-10, "maxiter": 800})
    fine = Rearrangement(r_nodes=48, th_nodes=192, x2_nodes=128)
    rearr = {}
    for name, (lam, mu, kap) in {"h=x1+sgn(x1)(lam x2^2+kap x1^2)": (abs(r2.x[0]), 0.0, abs(r2.x[1])),
                                 "h=x1+sgn(x1)(lam x2^2+mu x2^4+kap x1^2)":
                                     (abs(r3.x[0]), abs(r3.x[2]), abs(r3.x[1]))}.items():
        base, ref = rea.evaluate(lam, mu, kap), fine.evaluate(lam, mu, kap)
        check(abs(base["Q"] - ref["Q"]) < 1e-6, f"rearrangement quadrature convergence {name}")
        check(abs(ref["E_Y2_minus_1"]) < 1e-7 and abs(ref["E_Y4_minus_3"]) < 1e-6,
              f"rearrangement marginal residuals {name}")
        check(ref["Q"] >= lower_co, f"rearrangement respects comonotone lower bound {name}")
        rearr[name] = {"lam": lam, "mu": mu, "kap": kap, "Q": ref["Q"],
                       "Q_coarse_grid": base["Q"], "E_Y2_minus_1": ref["E_Y2_minus_1"],
                       "E_Y4_minus_3": ref["E_Y4_minus_3"], "E_X1Y": ref["E_X1Y"],
                       "E_X2Y": ref["E_X2Y"], "quadrature": [rea.label, fine.label]}
    best = min(r["Q"] for r in rearr.values())
    check(lower_dds < lower_co <= best < 12, "bracket ordering DDS < comonotone <= best numerical < 12")
    check(fine_c["Q_min"] >= lower_co, "flow minimum respects comonotone lower bound")
    table = [
        {"bound": "lower", "value": lower_dds, "label": "certified (DDS / Cauchy-Schwarz): 12-8(sqrt2-1)"},
        {"bound": "lower", "value": lower_co, "label": "certified inequality, mpmath-evaluated: comonotone"},
        {"bound": "upper", "value": 76 / 7, "label": "certified: flow at t=1/7, K<28"},
        {"bound": "upper", "value": 12 - 4 / (1 + math.sqrt(6)), "label": "certified: 12-4/(1+sqrt6)"},
        {"bound": "upper", "value": fine_c["Q_at_1_over_7"], "label": "numerical: W flow at t=1/7"},
        {"bound": "upper", "value": fine_c["Q_min"],
         "label": f"numerical: W/pair flow minimum at t={fine_c['t_min']:.6f}"},
    ] + [{"bound": "upper", "value": r["Q"], "label": f"numerical: monotone rearrangement {k}"}
         for k, r in rearr.items()]
    return {"Q_G": 12.0, "table": table, "flow": flow,
            "rearrangement": rearr, "linear_control": lin,
            "enclosure_numerical": [lower_co, best],
            "label": "Upper values marked numerical are deterministic quadrature of feasible "
                     "couplings, not interval-certified.", "status": "passed"}


# ---------------------------------------------------------------------------
# (d) equal-dimension worked example and (e) realistic spectra
# ---------------------------------------------------------------------------

def row(name, a, b, bracket, comonotone_ok=True):
    sp = spectral(a, b)
    br = bracket(list(map(float, a)), list(map(float, b)))
    scale = max(1.0, sp["QG"])
    for key, mine in [("Q_G", sp["QG"]), ("K_C", sp["KC"]), ("S", sp["S"]), ("cos_theta", sp["cos"])]:
        check(abs(br[key] - mine) <= 1e-9 * max(1.0, abs(mine)), f"bracket vs numpy {key} {name}")
    check(abs(br["certificates"]["best_pair_32D2_over_Kij"] - sp["pair"]) <= 1e-9 * scale,
          f"pair certificate {name}")
    check(sp["low"] <= sp["sharp"] * (1 + 1e-12), f"K_C certificate sharper than S/(A2B2) {name}")
    check(sp["KC"] <= 32 * sp["A2"] * sp["B2"] * sp["S"] * (1 + 1e-12), f"K_C <= 32A2B2S {name}")
    out = {"d": len(sp["a"]), "cos_theta": sp["cos"], "Q_G": sp["QG"], "S": sp["S"], "K_C": sp["KC"],
           "lower_S_over_A2B2": sp["low"], "lower_32S2_over_KC": sp["sharp"],
           "lower_best_pair": sp["pair"], "upper_dds_8(A2B2-ab)": sp["up_dds"],
           "upper_half_QG": sp["QG"] / 2,
           "rel_lower_S_over_A2B2": sp["low"] / sp["QG"], "rel_lower_32S2_over_KC": sp["sharp"] / sp["QG"],
           "rel_lower_best_pair": sp["pair"] / sp["QG"],
           "rel_upper_dds": sp["up_dds"] / sp["QG"],
           "bracket_interval": list(br["Q_star_interval"])}
    if comonotone_ok:
        law_a, law_b = chi_law(sp["a"]), chi_law(sp["b"])
        mc, info = comonotone(law_a, law_b, 20 if isinstance(law_a, Imhof) or isinstance(law_b, Imhof) else 40)
        up_co = 4 * (mc - 2 * sp["ab"])
        err = 4 * info["quadrature_diagnostics"]["s_grid_16_vs_24_difference"]
        check(err < 1e-3 * max(up_co, 1e-300), f"comonotone s-grid convergence {name}: {err}")
        check(up_co <= sp["up_dds"] + 1e-9 * scale, f"comonotone upper gap <= DDS {name}")
        check(up_co >= sp["sharp"] - 1e-9 * scale, f"comonotone upper gap >= certified lower gap {name}")
        out.update({"upper_comonotone": up_co, "rel_upper_comonotone": up_co / sp["QG"],
                    "comonotone_method": info["methods"],
                    "comonotone_diagnostics": info["quadrature_diagnostics"],
                    "comonotone_tail_estimate": info["tail_estimate"]})
    return out


def part_d(bracket):
    a, b = [2.0, 1.0], [1.0, 1.0]
    r = row("a=(2,1),b=(1,1)", a, b, bracket)
    check(abs(r["Q_G"] - 12) < 1e-12 and abs(r["S"] - 1) < 1e-12, "worked example Q_G=12, S=1")
    check(abs(r["K_C"] - (math.sqrt(896) + 24)) < 1e-10, "worked example K_C = sqrt(896)+24")
    upper = min(r["upper_dds_8(A2B2-ab)"], r["upper_comonotone"], r["upper_half_QG"])
    lower = max(r["lower_32S2_over_KC"], r["lower_best_pair"], r["lower_S_over_A2B2"])
    r.update({"best_lower_gap": lower, "best_upper_gap": upper,
              "rel_best_lower": lower / r["Q_G"], "rel_best_upper": upper / r["Q_G"],
              "Q_star_interval": [r["Q_G"] - upper, r["Q_G"] - lower],
              "exact_forms": {"Q_G": "12", "S": "1", "A2B2": "sqrt(10)", "a.b": "3",
                              "S_over_A2B2": "1/sqrt(10)", "K_C": "sqrt(896)+24 = 8 sqrt(14)+24",
                              "32S2_over_KC": "32/(24+8 sqrt(14)) = 4/(3+sqrt(14))",
                              "upper_dds": "8(sqrt(10)-3)"},
              "label": "lower gaps certified; comonotone upper gap is a certified inequality "
                       "evaluated numerically", "status": "passed"})
    return r


def part_e(bracket):
    i = np.arange(1, 769, dtype=float)
    nrm = lambda x: np.asarray(x, float) / np.sum(x)
    rng = np.random.default_rng(1)
    cases = {
        "plane->line a=(1,1), b=(1)": ([1., 1.], [1.]),
        "power law d=50: i^-1 vs i^-1.5 (trace 1)": (nrm(i[:50] ** -1), nrm(i[:50] ** -1.5)),
        "exp decay d=50: e^-0.1i vs e^-0.2i (trace 1)": (nrm(np.exp(-.1 * i[:50])), nrm(np.exp(-.2 * i[:50]))),
        "PCA truncation: a=i^-1 (d=200), b=a[:10] (raw)": (i[:200] ** -1, (i[:200] ** -1)[:10]),
        "near-isotropic d=50: a=1, b=1+U(0,0.2), seed 1 (trace 1)":
            (nrm(np.ones(50)), nrm(np.sort(1 + 0.2 * rng.random(50))[::-1])),
        "embeddings 768 vs 512, both i^-1 (trace 1)": (nrm(i ** -1), nrm(i[:512] ** -1)),
    }
    rows = {}
    for name, (a, b) in cases.items():
        rows[name] = row(name, a, b, bracket)
    return {"rows": rows, "seed_near_isotropic": 1,
            "label": "certificates are proved inequalities in floating point; comonotone column is "
                     "a numerical evaluation (Imhof inversion for d>2)", "status": "passed"}


# ---------------------------------------------------------------------------
# (f) whole-path curvature check
# ---------------------------------------------------------------------------

def part_f():
    d2_fine, d2_coarse = (8.0, 0.125, 12), (8.0, 0.125, 8)
    d3_fine, d3_coarse = (6.5, 0.5, 6), (6.5, 0.65, 6)
    cases = [("plane->line a=(1,1), b=(1,0)", [1., 1.], [1., 0.]),
             ("a=(2,1), b=(1,1)", [2., 1.], [1., 1.]),
             ("a=(3,1), b=(1,1)", [3., 1.], [1., 1.]),
             ("a=(3,2,1), b=(1,1,1)", [3., 2., 1.], [1., 1., 1.]),
             ("zero eigenvalue a=(2,1,0), b=(1,1,1)", [2., 1., 0.], [1., 1., 1.]),
             ("a=(5,1,0.5), b=(1,0.8,0.1)", [5., 1., .5], [1., .8, .1]),
             ("line vs R^3 a=(1,1,1), b=(1,0,0)", [1., 1., 1.], [1., 0., 0.])]
    results = {}
    for name, a, b in cases:
        a, b = np.array(a), np.array(b)
        d = len(a)
        sp = spectral(a, b, pairs=False)
        S, KC = sp["S"], sp["KC"]
        tstar = 4 * S / KC
        n = 80 if d == 2 else 40
        h = 2 * tstar / n
        grid = np.arange(-n, n + 1) * h
        values = {}
        for tag, spec in (("fine", d2_fine if d == 2 else d3_fine), ("coarse", d2_coarse if d == 2 else d3_coarse)):
            rule = orthant_rule(d, *spec)
            Zg, Wg, label = rule
            Cg = np.outer(b, a) - np.outer(a, b)
            rec = {}

            def cb(t, T, rec=rec, Zg=Zg, Wg=Wg, Cg=Cg):
                rec[t] = (G_value(a, b, Zg, Wg, T), G_second_derivative_generator(a, b, Cg, Zg, Wg, T),
                          float(Wg @ (b @ (T * T))) - b.sum())
            flow_path(a, b, rule, grid, cb)
            keys = sorted(rec)
            check(len(keys) == len(grid), f"path grid bookkeeping {name}")
            values[tag] = (np.array([rec[k][0] for k in keys]), np.array([rec[k][1] for k in keys]),
                           max(abs(rec[k][2]) for k in keys), label, len(Wg))
        G, g2, marg, label, nodes = values["fine"]
        Gc, g2c, _, label_c, _ = values["coarse"]
        fd = (G[2:] - 2 * G[1:-1] + G[:-2]) / h ** 2
        max_fd, max_gen = float(np.max(np.abs(fd))), float(np.max(np.abs(g2)))
        G0 = G[n]
        taylor_margin = (G - G0) - (4 * S * grid - KC * grid ** 2 / 2)
        quad_change = float(np.max(np.abs(G - Gc)))
        quad_change_g2 = float(np.max(np.abs(g2 - g2c)))
        fd_vs_gen = float(np.max(np.abs(fd - g2[1:-1])))
        slope0 = (8 * (G[n + 1] - G[n - 1]) - (G[n + 2] - G[n - 2])) / (12 * h)   # fourth order
        nz = grid != 0
        margin_over_t2 = float(np.min(taylor_margin[nz] / grid[nz] ** 2))
        # Tolerances: rounding allowance plus twice the observed change between
        # the fine and coarse quadrature rules (an empirical error estimate).
        tol = 1e-6 * KC + 2 * quad_change_g2
        taylor_tol = 1e-9 * max(1.0, abs(G0)) + 2 * quad_change
        check(quad_change < 1e-4 * max(1.0, abs(G0)), f"path quadrature convergence {name}: {quad_change}")
        check(quad_change_g2 < 1e-3 * KC, f"path G'' quadrature convergence {name}: {quad_change_g2}")
        check(fd_vs_gen < 1e-3 * KC, f"finite-difference vs generator G'' {name}: {fd_vs_gen}")
        check(abs(slope0 - 4 * S) < 1e-4 * max(1.0, 4 * S), f"G'(0) = 4S {name}: {slope0}")
        check(max_fd <= KC + tol, f"|G''| <= K_C along path {name}: {max_fd} vs {KC}")
        check(max_gen <= KC + tol, f"generator |G''| <= K_C along path {name}")
        check(np.min(taylor_margin) >= -taylor_tol, f"Taylor lower bound {name}")
        check(marg < 1e-5 * max(1.0, b.sum()), f"target marginal E|Y|^2 preserved {name}: {marg}")
        results[name] = {"a": a.tolist(), "b": b.tolist(), "S": S, "K_C": KC, "t_star_4S_over_KC": tstar,
                         "t_range": [-2 * tstar, 2 * tstar], "grid_points": len(grid), "fd_step": h,
                         "quadrature": label, "nodes": nodes, "convergence_rule": label_c,
                         "max_abs_G2_finite_difference": max_fd, "max_abs_G2_generator": max_gen,
                         "ratio_max_G2_over_KC": max_gen / KC,
                         "max_fd_minus_generator": fd_vs_gen, "G_change_between_rules": quad_change,
                         "G2_change_between_rules": quad_change_g2,
                         "G_prime_0_fourth_order_difference": slope0, "four_S": 4 * S,
                         "min_taylor_margin_over_t2_excluding_0": margin_over_t2,
                         "min_taylor_margin": float(np.min(taylor_margin)),
                         "max_target_second_moment_error": marg,
                         "tolerance_G2": tol, "tolerance_taylor": taylor_tol}
    return {"cases": results,
            "max_ratio_G2_over_KC": max(r["ratio_max_G2_over_KC"] for r in results.values()),
            "tolerance": "max|G''| <= K_C + 1e-6 K_C + 2 x (fine-coarse change in G''); Taylor margin "
                         ">= -(1e-9 max(1,|G(0)|) + 2 x (fine-coarse change in G))",
            "ode": "DOP853 rtol 1e-11 atol 1e-12, integrated outward from t=0 on each branch",
            "symmetry": "orthant quadrature uses exact sign-flip symmetry; off-diagonal E[Z_i T_j] = 0",
            "label": "numerical (deterministic quadrature + ODE), not interval-certified",
            "status": "passed"}


def part_g(c=0.1, seed=20261101):
    """Near-isotropic spectra a = sort(1 + eps u), b = sort(1 + eps v), u, v centred."""
    rng = np.random.default_rng(seed)
    rows = []
    for d in (2, 5, 50, 500):
        for eps in (1e-2, 1e-3):
            u, v = rng.standard_normal(d), rng.standard_normal(d)
            u -= u.mean(); v -= v.mean()
            a, b = np.sort(1 + eps * u)[::-1], np.sort(1 + eps * v)[::-1]
            sp = spectral(a, b, pairs=False)
            r = sp["sharp"] / sp["QG"]
            upper = sp["up_dds"] / sp["QG"]
            target = d / (2 * (d + 2))
            check(r >= target - c * eps, f"near-isotropic lower ratio d={d} eps={eps}: {r} vs {target}")
            check(upper <= 0.5 + 1e-12, f"near-isotropic upper ratio d={d} eps={eps}: {upper}")
            check(r <= upper * (1 + 1e-12), f"lower ratio <= upper ratio d={d} eps={eps}")
            rows.append({"d": d, "eps": eps, "ratio_32S2_over_KC_over_QG": r, "d_over_2(d+2)": target,
                         "deviation": r - target, "deviation_over_eps": (r - target) / eps,
                         "upper_ratio_8(A2B2-ab)_over_QG": upper, "Q_G": sp["QG"], "cos_theta": sp["cos"]})
    return {"rows": rows, "seed": seed, "constant_c": c,
            "check": "r >= d/(2(d+2)) - c*eps and 8(A2B2-a.b)/Q_G <= 1/2",
            "max_abs_deviation_over_eps": max(abs(x["deviation_over_eps"]) for x in rows),
            "label": "certified inequalities evaluated in floating point; random near-isotropic draws",
            "status": "passed"}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=".", help="directory for extension_checks.json")
    parser.add_argument("--skip", default="", help="comma-separated parts to skip (a..g), for debugging")
    args = parser.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    skip = set(filter(None, args.skip.split(",")))
    started = time.perf_counter()
    bracket_mod = load_bracket()
    named = {"plane->line": ([1., 1.], [1.]), "(2,1) vs (1,1)": ([2., 1.], [1., 1.])}
    result = {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
              "mpmath": mpmath.__version__}
    timings = {}
    for key, fn in [("a_factor_two", lambda: part_a(named)), ("b_comonotone", part_b)]:
        if key[0] in skip:
            continue
        t0 = time.perf_counter(); result[key] = fn(); timings[key] = time.perf_counter() - t0
    M_co = result.get("b_comonotone", {}).get("M_co_plane_to_line")
    if M_co is None:
        M_co = float(mpmath_plane_to_line()[0])
    for key, fn in [("c_plane_to_line", lambda: part_c(M_co)),
                    ("d_equal_dim", lambda: part_d(bracket_mod.bracket)),
                    ("e_realistic", lambda: part_e(bracket_mod.bracket)),
                    ("f_path_curvature", part_f),
                    ("g_near_isotropic", part_g)]:
        if key[0] in skip:
            continue
        t0 = time.perf_counter(); result[key] = fn(); timings[key] = time.perf_counter() - t0
    result["timings_seconds"] = timings
    result["elapsed_seconds"] = time.perf_counter() - started
    result["explicit_checks_passed"] = CHECK_COUNT
    result["optimisation_flag"] = "-O (explicit checks active)" if not __debug__ else "none"
    result["skipped_parts"] = sorted(skip)
    # Status is assigned only after every explicit check above has passed.
    result["status"] = "PASS" if not skip else "PARTIAL (parts skipped)"
    (out / "extension_checks.json").write_text(json.dumps(result, indent=2, default=float) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k in
                      ("status", "elapsed_seconds", "explicit_checks_passed", "timings_seconds")}, indent=2))


if __name__ == "__main__":
    main()
