# Data and code for "Dependence of Baryonic Tully–Fisher Relation Outliers on Specific Star Formation Rate"

This repository contains the catalogue and the script used to produce every number in the paper.

## Files

| File | Description |
|---|---|
| `Data_For_Research.csv` | Cross-matched ALFALFA–SDSS catalogue, 1640 galaxies |
| `reproduce_results.py` | Recomputes all results from the raw columns of the CSV |

## How to reproduce

```bash
pip install numpy pandas scipy
python reproduce_results.py Data_For_Research.csv
```

The script recomputes every residual from the raw columns (W50, inclination, stellar mass, HI mass). It does not rely on the stored result columns.

## Method in brief

- Baryonic mass: `Mbar = M_star + 1.33 * M_HI`, with `M_star` from the Taylor et al. (2011) stellar masses (`logMstarTaylor`) and `M_HI` from `loghimass`.
- Observed velocity: `Vobs = W50 / (2 sin i)`.
- Expected velocity: `Vexp = (Mbar / A)^(1/4)`, with A = 50 (McGaugh 2005, primary) and A = 46.8 (Lelli et al. 2016, alternative).
- Residual: `dV = log10(Vobs / Vexp)`. A galaxy is an outlier if `|dV| >= 0.3` dex.
- sSFR: `log sSFR = logSFR22 - logMstarTaylor`, divided into 7 bins (`Bin 1` to `Bin 7`, lowest to highest).
- Test: Pearson chi-square test of homogeneity on the 2 x 7 table (outlier / not outlier, 7 bins, 6 degrees of freedom).

## Main results (what the script prints)

| Calibration | Outliers | chi-square | p |
|---|---|---|---|
| A = 50, recomputed with `|dV| >= 0.3` (used in the paper) | 179 | 31.47 | 2.06e-05 |
| A = 46.8 (Lelli et al. 2016) | 183 | 28.50 | 7.6e-05 |

Outliers per bin (A = 50): 18, 10, 25, 50, 46, 15, 15, out of 111, 136, 318, 558, 364, 103 and 50 galaxies.

## Known issue in the catalogue (please read):

The stored `is_outlier` column contains a flagging error for one galaxy, **AGC 114598** (`|dV| = 0.135` dex, in Bin 5). It is marked as an outlier although it does not meet the 0.3 dex criterion. The paper applies `|dV| >= 0.3` directly, so it counts **179** outliers (chi-square = 31.47). Using the stored flag would give 180 outliers (chi-square = 31.70). The script prints both lines so the difference is visible. The same galaxy is also flagged in the `Bin 5 Outliers` column.

Also, galaxy AGC 7647 has both `Bin 1` and `Bin 4` set to True. The script and the paper assign it to the first bin (Bin 1), which gives the bin sizes used in the paper (Bin 4: N = 558).

## Main columns

| Column | Meaning |
|---|---|
| `AGC_1` | ALFALFA (AGC) catalogue number |
| `RA`, `DEC` | Position (degrees) |
| `Vhelio` | Heliocentric velocity (km/s) |
| `Dist_1` | Distance (Mpc) |
| `w50` | HI line width at 50% of peak (km/s) |
| `Inc` | Inclination (degrees) |
| `loghimass`, `logMH_err` | log10 HI mass (solar masses) and its uncertainty |
| `logMstarTaylor`, `logMstarTaylor_err` | log10 stellar mass (Taylor et al. 2011) and its uncertainty |
| `logSFR22`, `logSFR22_err` | log10 star formation rate from WISE 22 micron |
| `sSFR` | log10 specific star formation rate (per year) |
| `Mbar` | Baryonic mass (solar masses) |
| `DeltaV_abs` | Signed velocity residual dV for A = 50 (dex). Despite the name it is signed; the sign is used for the positive/negative split |
| `is_outlier` | Stored outlier flag (see known issue above) |
| `Bin 1` to `Bin 7` | sSFR bin membership |
| `z` | Redshift |

Other columns are intermediate quantities used during the analysis and are not needed to reproduce the results.

## Data sources and references

- ALFALFA survey: Haynes et al. 2011, *Astron. J.* 142, 170.
- SDSS: Abazajian et al. 2009, *Astrophys. J. Suppl. Ser.* 182, 543.
- Stellar masses: Taylor et al. 2011, *Mon. Not. R. Astron. Soc.* 418, 1587.
- BTFR normalisation: McGaugh 2005, *Astrophys. J.* 632, 859 (A = 50); Lelli, McGaugh and Schombert 2016, *Astron. J.* 152, 157 (A = 46.8).

## Citation

If you use these data, please cite the paper: [Krishna Reddy Mittapelly], "Dependence of Baryonic Tully–Fisher Relation Outliers on Specific Star Formation Rate" ([Hosted as a preprint on Research square]).

The underlying survey data belong to the ALFALFA and SDSS collaborations; please cite the original surveys as well.
