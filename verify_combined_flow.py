"""Independent finite checks of the combined Gaussian-preserving flow formulas.

Uses physicists' Gauss-Hermite nodes transformed to standard normals, and evaluates
the generator directly. It does not import the original polynomial checker. Checks raise SystemExit
explicitly (active under ``python -O``). Output JSON goes to --out DIR.
"""
import argparse
from pathlib import Path
from fractions import Fraction
import itertools
import json
import platform
import time
import numpy as np


def check(condition, message):
    """Raise SystemExit (not AssertionError) so -O cannot strip the test."""
    if not condition:
        raise SystemExit(f"CHECK FAILED: {message}")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default='.',help='directory for combined_flow_verification.json')
    args=parser.parse_args()
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter()
    nodes,weights=np.polynomial.hermite.hermgauss(5)
    nodes=np.sqrt(2)*nodes
    weights=weights/np.sqrt(np.pi)
    rng=np.random.default_rng(73331)
    max_resid={k:0. for k in ['generator_norm','score_derivative','field_norm',
                               'weighted_divergence','lagrange','Ca_norm','Cb_norm']}
    cases=0
    for d in [2,3,4]:
        index=np.array(list(itertools.product(range(5),repeat=d)))
        z=nodes[index]
        w=np.prod(weights[index],axis=1)
        for sample in range(40):
            a=np.sort(rng.integers(0,7,d))[::-1].astype(float)
            b=np.sort(rng.integers(0,7,d))[::-1].astype(float)
            C=np.outer(b,a)-np.outer(a,b)
            eta=z*z-1
            velocity=z*(eta@C.T)
            # Direct generator evaluation, without using its Hoeffding expansion.
            La=2*np.sum(a*z*velocity,axis=1)
            Lb=2*np.sum(b*z*velocity,axis=1)
            Ha=np.sum(a*z*z,axis=1)
            pairs=[(i,j) for i in range(d) for j in range(i+1,d)]
            P=lambda v:8*((C@v)@(C@v))+16*sum(C[i,j]**2*(v[i]-v[j])**2 for i,j in pairs)
            p_a,p_b=P(a),P(b)
            A2,B2=a@a,b@b
            A,B=np.sqrt(A2),np.sqrt(B2)
            dot=a@b
            S=A2*B2-dot*dot
            vals={
                'generator_norm':abs(w@(La*La)-p_a),
                'score_derivative':abs(w@(Ha*Lb)-4*S),
                'field_norm':np.max(abs(w@(velocity*velocity)-2*np.sum(C*C,axis=1))),
                'weighted_divergence':np.max(abs(np.sum(eta@C.T,axis=1)-np.sum(z*velocity,axis=1))),
                'lagrange':abs(sum(C[i,j]**2 for i,j in pairs)-S),
                'Ca_norm':abs((C@a)@(C@a)-A2*S),
                'Cb_norm':abs((C@b)@(C@b)-B2*S),
            }
            scale=max(1.,p_a,p_b,A2*S,B2*S,abs(S))
            for k,v in vals.items():
                max_resid[k]=max(max_resid[k],float(v/scale))
                check(v/scale<2e-12,f"{k} residual d={d} sample={sample}: {v} (scale {scale})")
            K=np.sqrt(p_a*p_b)+8*sum(C[i,j]**2*(a[i]*b[i]+a[j]*b[j]) for i,j in pairs)
            check(p_a <= 24*A2*S+1e-7,f"P_C(a) bound d={d} sample={sample}")
            check(p_b <= 24*B2*S+1e-7,f"P_C(b) bound d={d} sample={sample}")
            check(K <= 32*A*B*S+1e-7,f"K_C <= 32 A2 B2 S d={d} sample={sample}")
            cases+=1
    old_gain=4*Fraction(1,9)-Fraction(36,2)*Fraction(1,9)**2
    new_gain=4*Fraction(1,7)-Fraction(28,2)*Fraction(1,7)**2
    check(12-4*old_gain==Fraction(100,9),"old rational bound 100/9")
    check(12-4*new_gain==Fraction(76,7),"new rational bound 76/7")
    check(cases==120,f"expected 120 cases, ran {cases}")
    result={'python':platform.python_version(),'numpy':np.__version__,
            'cases':cases,'dimensions':[2,3,4],
            'quadrature_order_each_coordinate':5,
            'spectral_entries':'sorted integers 0 through 6, deterministic seed 73331',
            'max_scaled_residuals':max_resid,
            'old_rational_bound':str(12-4*old_gain),
            'new_rational_bound':str(12-4*new_gain),
            'elapsed_seconds':time.perf_counter()-start,
            'optimisation_flag':'-O (explicit checks active)' if not __debug__ else 'none'}
    # Status is assigned only after every explicit check above has passed.
    result['status']='All finite checks passed. General proofs remain analytic.'
    (out/'combined_flow_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
