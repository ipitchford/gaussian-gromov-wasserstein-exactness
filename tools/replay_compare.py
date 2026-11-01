# SPDX-License-Identifier: MIT
"""Clean replay of the package (read-only on ROOT), run in a temporary copy.

Usage: python -I -B tools/replay_compare.py ROOT [--quick] [--keep DIR]

Runs every checker with the interpreter running this script and flags -I -O -B
(the same commands as build.sh, without the PDF build), each writing its retained
output file names in the copy. Requires every checker to exit 0, then compares each
regenerated JSON output (and JSON-valued .log output) with the retained one:
numbers within rel 1e-7 / abs 1e-10, everything else exactly; run-time,
interpreter and optimisation-flag fields are ignored. CSV numbers use the same
tolerance; figures are regenerated, not byte-compared. Fails if a checker writes an
unexpected file or if ROOT changes.

--quick omits the two slow numerical checkers (extension_checks.py, about six
minutes, and the quadrature/ODE part of verify_gw.py, for which only the exact
rational part is run, into a scratch directory that is not compared).
"""
import csv, hashlib, io, json, math, os, shutil, subprocess, sys, tempfile
from pathlib import Path

root = Path(sys.argv[1]).resolve()
quick = '--quick' in sys.argv
keep = Path(sys.argv[sys.argv.index('--keep') + 1]) if '--keep' in sys.argv else None
VOLATILE = {'python', 'numpy', 'scipy', 'mpmath', 'matplotlib', 'platform', 'interpreter', 'working_directory',
            'elapsed_seconds', 'seconds', 'timings_seconds', 'optimisation_flag'}

# (script, arguments, stdout file relative to the package root or None)
FULL = [('verify_gw.py', ['--out', 'results'], 'results/run_log.json'),
        ('independent_symbolic_audit.py', [], 'independent_symbolic_audit.json'),
        ('verify_combined_flow.py', ['--out', '.'], 'combined_flow_run.log'),
        ('proof_digest_checks.py', ['--out', '.'], 'proof_digest_checks_run.log'),
        ('gw_gaussian_bracket.py', [], 'gw_gaussian_bracket_demo.log'),
        ('extension_checks.py', ['--out', '.'], 'extension_checks_run.log'),
        ('negative_controls.py', ['--out', '.'], 'negative_controls_run.log')]


def tree(p):
    return {q.relative_to(p).as_posix(): hashlib.sha256(q.read_bytes()).hexdigest()
            for q in p.rglob('*') if q.is_file() and not q.relative_to(p).as_posix().startswith('.git/')}


def cmp(a, b, path, diffs, stats):
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            diffs.append(f'{path}: keys differ {sorted(set(a) ^ set(b))[:5]}')
        for k in set(a) & set(b):
            if k in VOLATILE:
                stats['volatile'] += 1
                continue
            cmp(a[k], b[k], f'{path}.{k}', diffs, stats)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            diffs.append(f'{path}: length {len(a)} vs {len(b)}')
            return
        for i, (x, y) in enumerate(zip(a, b)):
            cmp(x, y, f'{path}[{i}]', diffs, stats)
    elif isinstance(a, bool) or isinstance(b, bool) or not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
        stats['exact'] += 1
        if a != b:
            diffs.append(f'{path}: {str(a)[:60]!r} vs {str(b)[:60]!r}')
    else:
        stats['numeric'] += 1
        if a == b:
            return
        if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
            return
        if not math.isclose(a, b, rel_tol=1e-7, abs_tol=1e-10):
            diffs.append(f'{path}: {a!r} vs {b!r}')
        else:
            stats['last_place'] += 1


def parse(text):
    """JSON document, or a one-line header followed by JSON; None if neither."""
    try:
        return '', json.loads(text)
    except ValueError:
        pass
    if '\n' in text:
        header, body = text.split('\n', 1)
        try:
            return header, json.loads(body)
        except ValueError:
            pass
    return None


before = tree(root)
with tempfile.TemporaryDirectory(prefix='gw-replay-') as td:
    w = Path(td) / 'pkg'
    shutil.copytree(root, w, ignore=shutil.ignore_patterns('.git'))
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg', MPLCONFIGDIR=str(Path(td) / 'mpl'),
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    plan = [c for c in FULL if not (quick and c[0] == 'extension_checks.py')]
    if quick:
        plan[0] = ('verify_gw.py', ['--exact-only', '--out', str(Path(td) / 'exact-only')], None)
    for script, args, stdout in plan:
        r = subprocess.run([sys.executable, '-I', '-O', '-B', script] + args, cwd=w, capture_output=True,
                           text=True, env=env, timeout=3600)
        if r.returncode:
            sys.exit(f'FAIL {script} exit {r.returncode}: {r.stderr[-2000:]}')
        if stdout:
            (w / stdout).write_text(r.stdout)
        print(f'ok {script} {" ".join(args)} (-I -O -B)')
    after = tree(w)
    extra = sorted(set(after) - set(before))
    if extra:
        sys.exit(f'FAIL replay wrote unexpected files: {extra}')
    changed = sorted(k for k in before if after.get(k) != before[k])
    diffs, stats, compared = [], dict(exact=0, numeric=0, last_place=0, volatile=0), []
    for k in changed:
        if k.endswith(('.json', '.log')):
            pa, pb = parse((root / k).read_text()), parse((w / k).read_text())
            if pa is None or pb is None:
                diffs.append(f'{k}: non-JSON text output changed')
                continue
            if pa[0] != pb[0]:
                diffs.append(f'{k}: header line differs')
            cmp(pa[1], pb[1], k, diffs, stats)
            compared.append(k)
        elif k.endswith('.csv'):
            ra = list(csv.reader(io.StringIO((root / k).read_text())))
            rb = list(csv.reader(io.StringIO((w / k).read_text())))
            if [len(r) for r in ra] != [len(r) for r in rb]:
                diffs.append(f'{k}: CSV shape differs')
                continue
            for i, (x, y) in enumerate(zip([c for r in ra for c in r], [c for r in rb for c in r])):
                try:
                    fx, fy = float(x), float(y)
                except ValueError:
                    cmp(x, y, f'{k}[{i}]', diffs, stats)
                else:
                    cmp(fx, fy, f'{k}[{i}]', diffs, stats)
            compared.append(k)
        elif not k.endswith(('.png', '.pdf')):
            diffs.append(f'{k}: unexpected output changed')
    if keep:
        shutil.copytree(w, keep)
    print(f'changed outputs (compared): {compared}')
    print(f'leaves exact={stats["exact"]} numeric={stats["numeric"]} (within tolerance but not identical: '
          f'{stats["last_place"]}); volatile fields skipped {stats["volatile"]}')
    if diffs:
        print('\n'.join(diffs[:40]))
        sys.exit(f'FAIL {len(diffs)} differences beyond tolerance')
if tree(root) != before:
    sys.exit('FAIL root changed')
mode = 'quick (exact checkers, combined-flow, proof-digest, bracket demo, negative controls)' if quick else 'full (all seven checkers)'
print(f'PASS clean replay, {mode}: every checker passed under -I -O -B; regenerated outputs equal the retained ones '
      'within rel 1e-7 (run-time and interpreter fields ignored; figures regenerated, not byte-compared); '
      'no unexpected files; root unchanged')
