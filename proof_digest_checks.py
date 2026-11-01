"""Independent exact boundary example and two-ODE flow-reduction checks.

Checks raise SystemExit explicitly (active under ``python -O``). Output JSON
goes to --out DIR (default: current directory).
"""
import argparse
from pathlib import Path
from fractions import Fraction as F
import hashlib
import itertools
import json
import platform
import time
import numpy as np
from scipy.integrate import solve_ivp


def check(condition, message):
    """Raise SystemExit (not AssertionError) so -O cannot strip the test."""
    if not condition:
        raise SystemExit(f"CHECK FAILED: {message}")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default='.',help='directory for proof_digest_checks.json')
    args=parser.parse_args()
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    square=list(itertools.product([-1,1],repeat=2))
    deterministic=[]
    for signs in itertools.product([-1,1],repeat=4):
        if sum(signs):
            continue
        kk=[sum(F(x[j]*y,4) for x,y in zip(square,signs)) for j in range(2)]
        score=F(2)+2*sum(k*k for k in kk)
        q=F(0)
        for i,xi in enumerate(square):
            for j,xj in enumerate(square):
                dx=sum((xi[k]-xj[k])**2 for k in range(2))
                dy=(signs[i]-signs[j])**2
                q+=F((dx-dy)**2,16)
        check(q==24-4*score,f"moment reduction Q=24-4G for signs {signs}")
        deterministic.append({'signs':signs,'covariance_squared':str(sum(k*k for k in kk)),
                              'Q':str(q)})
    check(len(deterministic)==6,"six balanced deterministic sign couplings")
    check(min(F(c['Q']) for c in deterministic)==8,"finite-square deterministic optimum equals 8")
    grid_count=0
    for vals in itertools.product([F(-1),F(-1,2),F(0),F(1,2),F(1)],repeat=4):
        if sum(vals): continue
        kk=[sum(F(x[j],4)*v for x,v in zip(square,vals)) for j in range(2)]
        check(sum(k*k for k in kk)<=sum(v*v for v in vals)/4<=1,f"Bessel bound for grid values {vals}")
        grid_count+=1

    rng=np.random.default_rng(44031)
    reduced_count=0
    max_error=0.
    max_rel=0.
    for dim in [2,3,7]:
        for sample in range(10):
            a=np.sort(rng.uniform(0,2,dim))[::-1]
            b=np.sort(rng.uniform(0,2,dim))[::-1]
            if sample==0: b=2*a
            if sample==1: a[-1]=0
            z=rng.normal(size=dim)
            if sample==2: z[0]=0
            if sample==3: z[:]=0
            C=np.outer(b,a)-np.outer(a,b)
            alpha,beta=a.sum(),b.sum()
            def full(_,zz): return zz*(C@(zz*zz-1))
            def reduced(_,uv):
                tt=z*np.exp(b*uv[0]-a*uv[1])
                return np.array([a@(tt*tt)-alpha,b@(tt*tt)-beta])
            for orientation in [-1,1]:
                end=orientation*.2/max(1,np.linalg.norm(a)*np.linalg.norm(b))
                times=np.linspace(0,end,41)
                fsol=solve_ivp(full,(0,end),z,t_eval=times,rtol=2e-12,atol=2e-13,method='DOP853')
                rsol=solve_ivp(reduced,(0,end),[0.,0.],t_eval=times,rtol=2e-12,atol=2e-13,method='DOP853')
                check(fsol.success and rsol.success,"ODE solver success")
                rebuilt=z[:,None]*np.exp(b[:,None]*rsol.y[0]-a[:,None]*rsol.y[1])
                error=float(np.max(abs(fsol.y-rebuilt)))
                relative=error/max(1,float(np.max(abs(fsol.y))))
                max_error=max(max_error,error)
                max_rel=max(max_rel,relative)
                check(relative<2e-10,f"flow reconstruction dim={dim} sample={sample} orientation={orientation}: {error}")
                reduced_count+=1
    result={'python':platform.python_version(),'numpy':np.__version__,
            'finite_square_deterministic_couplings':deterministic,
            'finite_square_conditional_mean_grid_cases':grid_count,
            'finite_square_exact_optimum':'8',
            'reduced_flow_cases':reduced_count,'dimensions':[2,3,7],
            'seed':44031,'rtol':2e-12,'atol':2e-13,
            'max_flow_reconstruction_abs_error':max_error,
            'max_flow_reconstruction_scaled_error':max_rel,
            'elapsed_seconds':time.perf_counter()-started,
            'optimisation_flag':'-O (explicit checks active)' if not __debug__ else 'none'}
    check(reduced_count==60,f"expected 60 reduced-flow cases, ran {reduced_count}")
    # Status is assigned only after every explicit check above has passed.
    result['status']='All checks passed. The finite-square optimality and flow reduction have independent analytic proofs.'
    (out/'proof_digest_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
