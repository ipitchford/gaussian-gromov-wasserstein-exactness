"""Independent exact polynomial audit; Python standard library only.

Checks algebraic identities used in the Gaussian GW flow proof. Gaussian
expectations integrate monomials by exact integer moments, not simulation.
This supplements, and does not replace, the analytic proof.
Checks raise SystemExit explicitly, so they remain active under ``python -O``.
"""

from fractions import Fraction
import json
from math import prod


def check(condition, message):
    """Raise SystemExit (not AssertionError) so -O cannot strip the test."""
    if not condition:
        raise SystemExit(f"CHECK FAILED: {message}")


class Poly:
    def __init__(self, terms, n):
        self.n = n
        self.terms = {tuple(k): Fraction(v) for k, v in terms.items() if v}

    @classmethod
    def const(cls, value, n):
        return cls({(0,) * n: value}, n)

    @classmethod
    def var(cls, index, n):
        powers = [0] * n
        powers[index] = 1
        return cls({tuple(powers): 1}, n)

    def __add__(self, other):
        if not isinstance(other, Poly):
            other = Poly.const(other, self.n)
        terms = dict(self.terms)
        for powers, coefficient in other.terms.items():
            terms[powers] = terms.get(powers, 0) + coefficient
        return Poly(terms, self.n)

    __radd__ = __add__

    def __neg__(self):
        return Poly({p: -v for p, v in self.terms.items()}, self.n)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if not isinstance(other, Poly):
            other = Poly.const(other, self.n)
        terms = {}
        for p, a in self.terms.items():
            for q, b in other.terms.items():
                r = tuple(i + j for i, j in zip(p, q))
                terms[r] = terms.get(r, 0) + a * b
        return Poly(terms, self.n)

    __rmul__ = __mul__

    def __pow__(self, degree):
        ans = Poly.const(1, self.n)
        for _ in range(degree):
            ans = ans * self
        return ans

    def diff(self, index):
        terms = {}
        for powers, coefficient in self.terms.items():
            if powers[index]:
                q = list(powers)
                q[index] -= 1
                terms[tuple(q)] = coefficient * powers[index]
        return Poly(terms, self.n)

    def gaussian_expectation(self, coordinates):
        terms = {}
        for powers, coefficient in self.terms.items():
            if any(powers[i] % 2 for i in coordinates):
                continue
            q = list(powers)
            for i in coordinates:
                coefficient *= prod(range(1, powers[i], 2))
                q[i] = 0
            terms[tuple(q)] = terms.get(tuple(q), 0) + coefficient
        return Poly(terms, self.n)

    def is_zero(self):
        return not self.terms


def audit():
    checked = []

    def equal(label, left, right=0):
        check((left - right).is_zero(), label)
        checked.append(label)

    n = 6
    u, v, ai, aj, bi, bj = [Poly.var(i, n) for i in range(n)]
    V = [u * (v**2 - 1), v * (1 - u**2)]

    def L(f):
        return f.diff(0) * V[0] + f.diff(1) * V[1]

    def E(f):
        return f.gaussian_expectation([0, 1])

    Pa = ai * u**2 + aj * v**2
    Pb = bi * u**2 + bj * v**2
    equal("Gaussian divergence", V[0].diff(0) + V[1].diff(1), u * V[0] + v * V[1])
    equal("Radial derivative", L(u**2 + v**2), 2 * (v**2 - u**2))
    equal("Mixed fourth moment derivative", E(Pa * L(Pb)), 4 * (bi * aj - bj * ai))
    equal("Generator norm", E(L(Pa)**2), 24 * (ai**2 + aj**2) - 32 * ai * aj)
    equal("First coordinate field norm", E(V[0]**2), 2)
    equal("Second coordinate field norm", E(V[1]**2), 2)
    equal("First covariance derivative", E(u * V[0]))
    equal("Second covariance derivative", E(v * V[1]))
    equal("First covariance second derivative", E(u * L(V[0])), -2)
    equal("Second covariance second derivative", E(v * L(V[1])), -2)
    equal("Skew-adjoint generator identity", E(Pa * L(L(Pb))), -E(L(Pa) * L(Pb)))

    d = 4
    n = 3 * d
    z = [Poly.var(i, n) for i in range(d)]
    a = [Poly.var(d+i, n) for i in range(d)]
    b = [Poly.var(2*d+i, n) for i in range(d)]
    Ha = sum(a[i] * z[i]**2 for i in range(d))
    Hb = sum(b[i] * z[i]**2 for i in range(d))
    Phi = Fraction(1, 2) * (Ha * Hb + sum(a) * Hb - sum(b) * Ha)
    W = [z[i] * sum((b[i]*a[j]-b[j]*a[i])*(z[j]**2-1)
                    for j in range(d) if j != i) for i in range(d)]
    for i in range(d):
        equal(f"Gradient decomposition coordinate {i}", 2*Ha*b[i]*z[i], Phi.diff(i)+W[i])
    equal("Full field Gaussian divergence", sum(W[i].diff(i)-z[i]*W[i] for i in range(d)))
    equal("Full field squared norm",
          sum(w**2 for w in W).gaussian_expectation(range(d)),
          4*sum((b[i]*a[j]-b[j]*a[i])**2 for i in range(d) for j in range(i+1, d)))
    equal("Lagrange spectral identity",
          sum((b[i]*a[j]-b[j]*a[i])**2 for i in range(d) for j in range(i+1, d)),
          sum(x**2 for x in a)*sum(x**2 for x in b)-sum(a[i]*b[i] for i in range(d))**2)
    equal("Full field radial cancellation", sum(z[i]*W[i] for i in range(d)),
          -sum(z[i]**2*sum(b[i]*a[j]-b[j]*a[i] for j in range(d)) for i in range(d)))
    C = [[b[i]*a[j]-a[i]*b[j] for j in range(d)] for i in range(d)]
    Ca = [sum(C[i][j]*a[j] for j in range(d)) for i in range(d)]
    Cb = [sum(C[i][j]*b[j] for j in range(d)) for i in range(d)]
    S = sum(C[i][j]**2 for i in range(d) for j in range(i+1, d))
    LHa = sum(Ha.diff(i)*W[i] for i in range(d))
    LHb = sum(Hb.diff(i)*W[i] for i in range(d))
    equal("Dimension-free source generator norm", (LHa**2).gaussian_expectation(range(d)),
          8*sum(v**2 for v in Ca)+16*sum(C[i][j]**2*(a[i]-a[j])**2
                                      for i in range(d) for j in range(i+1, d)))
    equal("Dimension-free target generator norm", (LHb**2).gaussian_expectation(range(d)),
          8*sum(v**2 for v in Cb)+16*sum(C[i][j]**2*(b[i]-b[j])**2
                                      for i in range(d) for j in range(i+1, d)))
    equal("Rank-two source contraction", sum(v**2 for v in Ca), sum(v**2 for v in a)*S)
    equal("Rank-two target contraction", sum(v**2 for v in Cb), sum(v**2 for v in b)*S)
    for i in range(d):
        equal(f"Full field coordinate norm {i}", (W[i]**2).gaussian_expectation(range(d)),
              2*sum(C[i][j]**2 for j in range(d)))
    equal("Full flow mixed-moment derivative", (Ha*LHb).gaussian_expectation(range(d)), 4*S)
    check(len(checked) == 28, f"expected 28 identities, ran {len(checked)}")
    # Status is assigned only here, after every explicit check has passed.
    return {"status": "PASS", "exact_identity_count": len(checked), "checks": checked,
            "method": "Exact Fraction polynomial arithmetic and Gaussian integer moments",
            "limitations": "No numerical approximation; does not verify analytic flow or completeness arguments.",
            "optimisation_flag": "-O (explicit checks active)" if not __debug__ else "none"}


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
