#!/usr/bin/env python3
"""
Bootstrap, leave-out and random-subsample tests for the BTFR outlier fractions.

Usage:  python bootstrap_analysis.py Data_For_Research - Dependece of Baryonic Tully Fisher Relation On Star Formation Rate.csv
Needs:  numpy, pandas, scipy
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

A = 50.0
THRESHOLD = 0.3
N_BOOT = 2000
N_SUB = 5000
SEED = 42

path = sys.argv[1] if len(sys.argv) > 1 else "Data_For_Research - Dependece of Baryonic Tully Fisher Relation On Star Formation Rate.csv"
d = pd.read_csv(path)
d.columns = [c.strip() for c in d.columns]

mbar = 10 ** d["logMstarTaylor"] + 1.33 * 10 ** d["loghimass"]
vobs = d["w50"] / (2 * np.sin(np.radians(d["Inc"])))
dv = np.log10(vobs / (mbar / A) ** 0.25).values
abs_dv = np.abs(dv)
inc = d["Inc"].values
sbin = d[[f"Bin {i}" for i in range(1, 8)]].values.argmax(axis=1)
n = len(d)


def analyse(idx):
    """Outlier fractions, chi-square p, U-shape flag and curvature for a set of galaxies."""
    b = sbin[idx]
    out = abs_dv[idx] >= THRESHOLD
    total = np.bincount(b, minlength=7)
    k = np.bincount(b, weights=out, minlength=7)
    if (total == 0).any() or k.sum() == 0:
        return None
    frac = k / total
    _, p, _, _ = stats.chi2_contingency(np.array([k, total - k]), correction=False)
    lowest_middle = frac[1:4].min()
    u_shape = frac[0] > lowest_middle and frac[5] > lowest_middle and frac[6] > lowest_middle
    x = np.arange(7) - 3
    design = np.vstack([x ** 2, x, np.ones(7)]).T
    w = np.sqrt(total)
    curvature = np.linalg.lstsq(design * w[:, None], frac * w, rcond=None)[0][0]
    return frac, p, u_shape, curvature


rng = np.random.default_rng(SEED)

# bootstrap
runs = []
while len(runs) < N_BOOT:
    r = analyse(rng.integers(0, n, n))
    if r is not None:
        runs.append(r)
frac = np.array([r[0] for r in runs])
p = np.array([r[1] for r in runs])
u = np.array([r[2] for r in runs])
curv = np.array([r[3] for r in runs])

print(f"Bootstrap: {N_BOOT} resamples of {n} galaxies, threshold {THRESHOLD} dex")
print(f"  U-shape found in {100 * u.mean():.1f} per cent of resamples")
print(f"  positive curvature in {100 * (curv > 0).mean():.1f} per cent")
print(f"  chi-square p < 0.05 in {100 * (p < 0.05).mean():.1f} per cent, p < 0.001 in {100 * (p < 0.001).mean():.1f} per cent")
print(f"  median p = {np.median(p):.2g}")
print("  bin   median   2.5 per cent   97.5 per cent")
for i in range(7):
    lo, med, hi = np.percentile(100 * frac[:, i], [2.5, 50, 97.5])
    print(f"   {i + 1}    {med:5.1f}    {lo:5.1f}          {hi:5.1f}")

# leave-out of the most extreme residuals
order = np.argsort(-abs_dv)
for m in (10, 20):
    keep = np.setdiff1d(np.arange(n), order[:m])
    f, pm, um, _ = analyse(keep)
    print(f"Leave-out of the {m} largest |dV|: p = {pm:.2g}, U-shape = {um}, fractions (per cent) = {np.round(100 * f, 1)}")

# restricted inclination sample against random subsamples of the same size
inside = np.where((inc > 35) & (inc < 70))[0]
f_in, p_in, _, _ = analyse(inside)
subs = []
while len(subs) < N_SUB:
    r = analyse(rng.choice(n, len(inside), replace=False))
    if r is not None:
        subs.append(r[1])
subs = np.array(subs)
print(f"35 < i < 70: {len(inside)} galaxies, p = {p_in:.3f}")
print(f"Random subsamples of {len(inside)} galaxies: median p = {np.median(subs):.2g}, "
      f"p < 0.05 in {100 * (subs < 0.05).mean():.1f} per cent, p >= {p_in:.3f} in {100 * (subs >= p_in).mean():.1f} per cent")
