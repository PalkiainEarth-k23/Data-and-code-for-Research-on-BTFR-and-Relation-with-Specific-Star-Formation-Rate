#!/usr/bin/env python3
"""
Reproduces the numbers in "Dependence of Baryonic Tully-Fisher Relation Outliers
on Specific Star Formation Rate" from the catalogue CSV.

Usage:  python reproduce_results.py Data_For_Research.csv
Needs:  numpy, pandas, scipy
All residuals are recomputed from the raw columns (W50, Inc, M*, M_HI), not read
from the stored result columns, so the CSV can be checked independently.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

path = sys.argv[1] if len(sys.argv) > 1 else "Data_For_Research.csv"
d = pd.read_csv(path)
d.columns = [c.strip() for c in d.columns]

# ---- core quantities (paper Eqs. 2, 5, 6, 7, 13) -------------------------
A_PRIMARY, A_LELLI = 50.0, 46.8
Mstar = 10 ** d["logMstarTaylor"]
MHI = 10 ** d["loghimass"]
Mbar = Mstar + 1.33 * MHI
Vobs = d["w50"] / (2 * np.sin(np.radians(d["Inc"])))
dV = np.log10(Vobs / (Mbar / A_PRIMARY) ** 0.25)          # signed residual, A = 50
d["dV"] = dV
d["lsSFR"] = d["logSFR22"] - d["logMstarTaylor"]
d["bin"] = d[[f"Bin {i}" for i in range(1, 8)]].values.argmax(axis=1) + 1
N = d.groupby("bin").size().values

def per_bin(flag):
    return pd.Series(np.asarray(flag)).groupby(d["bin"].values).sum().values

def chi2_2x7(k):
    """Pearson test of homogeneity on the 2x7 table (outlier / not outlier)."""
    c, p, _, _ = stats.chi2_contingency(np.array([k, N - k]), correction=False)
    return c, p

def line(label, k):
    c, p = chi2_2x7(k)
    print(f"{label:44s} total={int(k.sum()):4d}  k={k.tolist()}  chi2={c:6.2f}  p={p:.2e}")

print("=" * 100)
print(f"Galaxies: {len(d)}   bins N = {N.tolist()}")
edges = d.groupby("bin")["lsSFR"].agg(["min", "max"]).round(2)
print("log sSFR range per bin:\n", edges.T.to_string())

print("\n--- Main result (Table 2, Eq. 16) ---")
line("A=50, from stored is_outlier flag", per_bin(d["is_outlier"]))
line("A=50, recomputed |dV| >= 0.3", per_bin(np.abs(dV) >= 0.3))
n_flag_mismatch = int((d["is_outlier"] & (np.abs(dV) < 0.3)).sum())
print(f"note: {n_flag_mismatch} galaxy has is_outlier=True although |dV| < 0.3 "
      f"(dV = {dV[d['is_outlier'] & (np.abs(dV) < 0.3)].round(3).tolist()}); "
      "effect on chi2 is 31.70 -> 31.47")

print("\n--- Thresholds (Sec. 3.3) ---")
for thr in (0.25, 0.30, 0.40):
    line(f"A=50, |dV| >= {thr}", per_bin(np.abs(dV) >= thr))

print("\n--- Alternative calibration (Appendix, Lelli et al. A = 46.8) ---")
# A smaller A raises V_exp, so the SIGNED residual falls by (1/4) log10(50/46.8).
shift = 0.25 * np.log10(A_PRIMARY / A_LELLI)
print(f"signed-residual shift = -{shift:.5f} dex")
line("A=46.8, shift applied to SIGNED dV (correct)", per_bin(np.abs(dV - shift) >= 0.3))
line("A=46.8, shift applied to |dV| (not valid)", per_bin((np.abs(dV) - shift) >= 0.3))
print("(the second line reproduces 173 / chi2 = 30.10; it wrongly lowers |dV| for the "
      f"{int((dV <= -0.3).sum())} negative-residual outliers)")

print("\n--- Inclination checks ---")
inc = d["Inc"]
print(f"fraction in 35-70 deg: {((inc >= 35) & (inc <= 70)).mean():.3f}; <35: {(inc < 35).mean():.3f}; >70: {(inc > 70).mean():.3f}")
for lab, m in [("Inc >= 25", inc >= 25), ("35 <= Inc <= 70", (inc >= 35) & (inc <= 70))]:
    k = pd.Series((np.abs(dV[m]) >= 0.3)).groupby(d["bin"][m].values).sum().reindex(range(1, 8), fill_value=0).values
    n = d[m].groupby("bin").size().reindex(range(1, 8), fill_value=0).values
    c, p, _, _ = stats.chi2_contingency(np.array([k, n - k]), correction=False)
    print(f"{lab:20s} N={int(m.sum()):4d} f={np.round(k / np.maximum(n, 1), 3).tolist()} chi2={c:.1f} p={p:.1e}")

print("\n--- Mass and gas fraction (Sec. 4) ---")
d["logMbar"] = np.log10(Mbar)
d["fgas"] = 1.33 * MHI / Mbar
d["out"] = np.abs(dV) >= 0.3
for col, q in (("logMbar", 4), ("fgas", 4)):
    g = d.groupby(pd.qcut(d[col], q), observed=True)["out"].agg(["size", "sum"])
    c, p, _, _ = stats.chi2_contingency(np.array([g["sum"], g["size"] - g["sum"]]), correction=False)
    print(f"outlier fraction by {col} quartile: {(g['sum'] / g['size']).round(3).tolist()}  chi2={c:.1f} p={p:.1e}")
m = d["logMbar"] >= 9
k = pd.Series(d["out"][m].values).groupby(d["bin"][m].values).sum().values
n = d[m].groupby("bin").size().values
c, p, _, _ = stats.chi2_contingency(np.array([k, n - k]), correction=False)
print(f"log Mbar >= 9: N={int(m.sum())} f={np.round(k / n, 3).tolist()} chi2={c:.1f} p={p:.1e}")

print("\n--- Error propagation (Sec. 2.6) ---")
print(f"sigma_dV (flat floors) = {np.hypot(0.05, 0.25 * 0.10):.4f} dex ; 0.3/sigma = {0.3 / np.hypot(0.05, 0.025):.2f}")
rs = 1.4826 * np.median(np.abs(dV - np.median(dV)))
from math import erfc, sqrt
print(f"robust spread of dV = {rs:.3f} dex ; Gaussian expectation beyond 0.3 = {len(d) * erfc(0.3 / (rs * sqrt(2))):.0f} galaxies ; observed = {int((np.abs(dV) >= 0.3).sum())}")

print("\n--- Distance (Sec. 3.5) ---")
sd = 0.5 * (d["sigDist"] / d["Dist_1"]) / np.log(10)
print(f"median z = {d['z'].median():.4f}; fraction z<0.01 = {(d['z'] < 0.01).mean():.3f}; distance contribution to sigma_dV: median {sd.median():.3f}, max {sd.max():.3f} dex")
db = pd.cut(d["Dist_1"], [0, 15, 30, 50, 80, 300])
print("outlier fraction by distance bin (Mpc):", d.groupby(db, observed=True)["out"].mean().round(3).to_dict())
