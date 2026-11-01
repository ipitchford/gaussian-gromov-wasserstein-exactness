#!/usr/bin/env python3
"""Reproducible exact algebra and numerical illustrations for Gaussian GW rigidity.

The polynomial checker uses only the Python standard library. Numerical flow
integration and quadrature require NumPy and SciPy; plots require Matplotlib.
None of these finite computations substitutes for the proofs in manuscript.md.

Every check uses an explicit test that raises SystemExit, so the script also
fails correctly under ``python -O``; status fields are written only after the
corresponding checks have run. negative_controls.py demonstrates detection.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import time
from fractions import Fraction
from pathlib import Path


def check(condition, message):
    """Raise SystemExit (not AssertionError) so -O cannot strip the test."""
    if not condition:
        raise SystemExit(f"CHECK FAILED: {message}")


def add(*polys):
    out = {}
    for p in polys:
        for e, c in p.items():
            out[e] = out.get(e, 0) + c
    return {e: c for e, c in out.items() if c}


def scale(p, c):
    return {e: c * v for e, v in p.items() if c * v}


def mul(p, q):
    out = {}
    for e, c in p.items():
        for f, d in q.items():
            g = tuple(x + y for x, y in zip(e, f))
            out[g] = out.get(g, 0) + c * d
    return {e: c for e, c in out.items() if c}


def diff(p, i):
    out = {}
    for e, c in p.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            out[tuple(f)] = c * e[i]
    return out


def normal_moment(n):
    if n % 2:
        return 0
    return math.prod(range(1, n, 2))


def expectation(p):
    return sum(c * math.prod(normal_moment(k) for k in e) for e, c in p.items())


def exact_checks():
    u, v, one = {(1, 0): 1}, {(0, 1): 1}, {(0, 0): 1}
    u2, v2 = mul(u, u), mul(v, v)
    V1 = mul(u, add(v2, scale(one, -1)))
    V2 = mul(v, add(one, scale(u2, -1)))
    def L(p):
        return add(mul(V1, diff(p, 0)), mul(V2, diff(p, 1)))
    checks = []
    def check_equal(name, actual, wanted):
        check(actual == wanted, f"{name}: got {actual!r}, wanted {wanted!r}")
        checks.append({"name": name, "status": "passed", "value": str(actual)})
    check_equal("Gaussian divergence identity", add(diff(V1, 0), diff(V2, 1),
          scale(add(mul(u, V1), mul(v, V2)), -1)), {})
    check_equal("radial derivative", L(add(u2, v2)), scale(add(v2, scale(u2, -1)), 2))
    for name, p in [("uV1", mul(u, V1)), ("vV1", mul(v, V1)),
                    ("uV2", mul(u, V2)), ("vV2", mul(v, V2))]:
        check_equal("covariance derivative " + name, expectation(p), 0)
    check_equal("V1 squared norm", expectation(mul(V1, V1)), 2)
    check_equal("V2 squared norm", expectation(mul(V2, V2)), 2)
    check_equal("L u2 squared norm", expectation(mul(L(u2), L(u2))), 24)
    check_equal("L v2 squared norm", expectation(mul(L(v2), L(v2))), 24)
    check_equal("L u2 L v2 cross moment", expectation(mul(L(u2), L(v2))), -16)
    check_equal("u L2 u", expectation(mul(u, L(L(u)))), -2)
    check_equal("v L2 v", expectation(mul(v, L(L(v)))), -2)
    # Four coefficients verify the general bilinear derivative, not just examples.
    check_equal("M' coefficient ai bi", expectation(mul(u2, L(u2))), 0)
    check_equal("M' coefficient ai bj", expectation(mul(u2, L(v2))), -4)
    check_equal("M' coefficient aj bi", expectation(mul(v2, L(u2))), 4)
    check_equal("M' coefficient aj bj", expectation(mul(v2, L(v2))), 0)
    check_equal("M'' coefficient ai bi", expectation(mul(u2, L(L(u2)))), -24)
    check_equal("M'' coefficient ai bj", expectation(mul(u2, L(L(v2)))), 16)
    check_equal("M'' coefficient aj bi", expectation(mul(v2, L(L(u2)))), 16)
    check_equal("M'' coefficient aj bj", expectation(mul(v2, L(L(v2)))), -24)
    # Exhaust all nonnegative integer two-coordinate spectra up to 5, including
    # zeros and repeated values. Check polynomial identities and gap-bound logic.
    cases = 0
    for ai in range(6):
        for aj in range(ai + 1):
            for bi in range(6):
                for bj in range(bi + 1):
                    P = add(scale(u2, ai), scale(v2, aj))
                    F = add(scale(u2, bi), scale(v2, bj))
                    delta = bi * aj - bj * ai
                    case = (ai, aj, bi, bj)
                    check(expectation(mul(P, L(F))) == 4 * delta,
                          f"grid mixed-moment derivative {case}")
                    check(expectation(mul(L(P), L(P))) == 24*(ai*ai+aj*aj)-32*ai*aj,
                          f"grid generator norm {case}")
                    check(expectation(mul(P, L(L(F)))) == 16*(ai*bj+aj*bi)-24*(ai*bi+aj*bj),
                          f"grid second derivative {case}")
                    cases += 1
    check(cases == 441, f"grid case count {cases} != 441")
    # Recorded only after every grid case has passed its three explicit checks.
    checks.append({"name": "sorted integer spectrum pairs, coordinates 0 through 5",
                   "status": "passed", "cases": cases})
    # Rational conservative curvature K<=28 yields the explicit counterexample.
    t, K = Fraction(1, 7), Fraction(28)
    score_gain = 4*t - K*t*t/2
    check_equal("certified score gain at t=1/7", score_gain, Fraction(2, 7))
    check_equal("certified objective upper bound", Fraction(12)-4*score_gain, Fraction(76, 7))
    return checks


def flow_quadrature(order, times, a=(1., 1.), b=(1., 0.)):
    import numpy as np
    from scipy.integrate import solve_ivp
    from scipy.special import roots_hermitenorm
    nodes, weights = roots_hermitenorm(order)
    weights = weights / np.sqrt(2*np.pi)
    u, v = np.meshgrid(nodes, nodes, indexing="ij")
    w = np.outer(weights, weights).ravel()
    u, v = u.ravel(), v.ravel()
    initial = np.concatenate((u, v))
    count = len(u)
    def rhs(_, z):
        x, y = z[:count], z[count:]
        return np.concatenate((x*(y*y-1), y*(1-x*x)))
    sol = solve_ivp(rhs, (0., max(times)), initial, method="DOP853",
                    t_eval=times, rtol=2e-11, atol=2e-12)
    check(sol.success, f"ODE solver failed: {sol.message}")
    P = a[0]*u*u+a[1]*v*v
    C0 = 4*((sum(a))**2+(sum(b))**2-sum(a)*sum(b)) + 8*(sum(x*x for x in a)+sum(x*x for x in b))
    rows = []
    for k, t in enumerate(sol.t):
        U, W = sol.y[:count, k], sol.y[count:, k]
        ru, rv = w@(u*U), w@(v*W)
        mixed = w@(P*(b[0]*U*U+b[1]*W*W))
        g = mixed + 2*(a[0]*b[0]*ru*ru+a[1]*b[1]*rv*rv)
        rows.append({"order": order, "time": float(t), "objective": float(C0-4*g),
                     "mixed_moment": float(mixed), "corr_u": float(ru), "corr_v": float(rv),
                     "marginal_second_error": float(max(abs(w@(U*U)-1), abs(w@(W*W)-1))),
                     "marginal_fourth_error": float(max(abs(w@(U**4)-3), abs(w@(W**4)-3))),
                     "gaussian_objective": float(4*((sum(a)-sum(b))**2+2*sum((x-y)**2 for x,y in zip(a,b))))})
    return rows


def numerical_checks(output, make_plot=True):
    import numpy as np
    import scipy
    times = sorted(set([0., 1/7, .02, .05, .08, .1, .125, .15]))
    rows = []
    for order in [64, 128, 256, 384]:
        rows.extend(flow_quadrature(order, times))
    with (output / "flow_quadrature.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    coarse = [r for r in rows if r["order"] == 256]
    fine = [r for r in rows if r["order"] == 384]
    delta = max(abs(a["objective"]-b["objective"]) for a, b in zip(coarse, fine))
    explicit = min(fine, key=lambda r: abs(r["time"]-1/7))
    check(explicit["objective"] < 76/7,
          "numerical explicit flow does not agree with certified upper bound")
    check(delta < 1e-3, f"quadrature orders 256/384 disagree by {delta}")
    for r in fine:
        check(r["marginal_second_error"] < 1e-5 and r["marginal_fourth_error"] < 1e-3,
              f"marginal moment residual too large at t={r['time']}")
    # Independent derivative check from high-order quadrature at small positive t.
    small = flow_quadrature(64, [0., .0001, .0002])
    derivative = (small[-1]["objective"]-small[0]["objective"])/.0002
    check(abs(derivative + 16) < 1e-2, f"finite-difference slope {derivative} != -16")
    # A proportional control has exact Gaussian optimum; perturbation must not improve it.
    control = flow_quadrature(64, [0., .05, .1], (2., 1.), (4., 2.))
    check(min(r["objective"] for r in control) >= control[0]["gaussian_objective"]-1e-8,
          "proportional control contradicted")
    summary = {"numpy": np.__version__, "scipy": scipy.__version__,
               "ode_method": "DOP853", "rtol": 2e-11, "atol": 2e-12,
               "quadrature_orders": [64, 128, 256, 384],
               "max_objective_difference_orders_256_384": delta,
               "explicit_example": explicit,
               "finite_difference_at_zero": derivative,
               "exact_derivative_at_zero": -16,
               "proportional_control": control,
               "scope": "Floating-point quadrature and ODE illustrations, not certified global optima.",
               "status": "passed"}
    if make_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                             "axes.spines.top": False, "axes.spines.right": False})
        fig, ax = plt.subplots(figsize=(7.1, 4.2), layout="constrained")
        tx = np.array([r["time"] for r in fine]); q = np.array([r["objective"] for r in fine])
        ax.plot(tx, q, "o-", color="#8B5D82", lw=2, label="Feasible nonlinear coupling (quadrature)")
        ax.axhline(12., color="#AD6B35", ls="--", lw=1.5, label="Best Gaussian coupling: 12")
        ax.scatter([1/7], [76/7], marker="D", color="#4F7C78", s=48,
                   zorder=6, label="Proved upper bound at t = 1/7: 76/7")
        ax.set(xlabel="Flow time t", ylabel="Squared-distance GW objective Q",
               title="A standard Gaussian plane can be matched to a line\nmore accurately by a nonlinear coupling")
        ax.legend(frameon=False, loc="lower left", fontsize=8)
        ax.grid(axis="y", color="#D7D1CC", alpha=.6)
        fig.savefig(output / "nonlinear_improvement.png", dpi=220)
        fig.savefig(output / "nonlinear_improvement.pdf")
        plt.close(fig)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="results")
    parser.add_argument("--exact-only", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    result = {"python": platform.python_version(), "exact_checks": exact_checks(),
              "formal_proof_assistant": "not used"}
    if not args.exact_only:
        result["numerical"] = numerical_checks(out, not args.no_plot)
    result["elapsed_seconds"] = time.perf_counter()-start
    result["optimisation_flag"] = "-O (asserts stripped; explicit checks active)" if not __debug__ else "none"
    result["status"] = "PASS"  # reached only if every explicit check above passed
    (out / "verification_results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
