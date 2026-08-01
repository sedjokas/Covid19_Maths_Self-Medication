# Code

Python scripts and data used to produce the model results, calibration, and figures in the manuscript (Sections 4 and 5).

## Layout

```
codes/
  model.py                 SEIHM-R-D model: ODE right-hand side, R0, center-manifold coefficient a
  run_all.py                Section 4 scenarios (1-10), baseline parameter set
  calibrate_uganda.py        Uganda Delta-wave calibration (differential evolution)
  finalize_uganda.py         Goodness-of-fit stats, 95% CI ensemble, gamma sensitivity, Section 5 figures
  transfer_countries.py      Cross-country parameter transfer (Mozambique, Senegal, Cameroon)
  data/
    uganda_calib_data.csv      Cleaned Uganda daily-case / cumulative-death series used for calibration
    transfer_data.json         Cleaned per-country cumulative-death series for the transfer exercise
  results/
    uganda_fit_results.json    Output of calibrate_uganda.py
    uganda_final_results.json  Output of finalize_uganda.py (goodness-of-fit, CIs)
    transfer_results.json      Output of transfer_countries.py
```

Run order: `calibrate_uganda.py` -> `finalize_uganda.py` -> `transfer_countries.py`. Each reads the previous script's JSON output from `results/`. `run_all.py` is independent (Section 4 only) and uses `model.py` directly.

Requires: `numpy`, `scipy`, `matplotlib`, `scikit-learn`.

## Data provenance

Case and death counts are from Our World in Data's COVID-19 dataset (`https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv`), cross-checked against a Wayback Machine snapshot archived November 2021 to rule out later revision. `data/uganda_calib_data.csv` is the cleaned Uganda subset (1 June - 31 October 2021); seven days (23-29 August 2021) are flagged and excluded in `calibrate_uganda.py` as a reporting-backfill artifact (see the script and the manuscript, Section 5.2, for details). `data/transfer_data.json` holds the corresponding cumulative-death series for Mozambique (Jul-Nov 2020), Senegal (Jul-Oct 2020), and Cameroon (Jun-Sep 2020), plus each country's day-0 case count used to seed the initial condition. Country populations are World Bank 2020 estimates.

## Notes on parameters fitted vs. held fixed

For the Uganda calibration, `alpha` (a four-parameter logistic: `alpha_hi`, `alpha_lo`, `k_sw`, `t_sw`), `kappa`, the initial infectious seed `I(0)`, and the mortality-scaling factor `delta` are fitted by differential evolution. `gamma` is held at its literature value (0.20) rather than fitted: case-count data alone barely constrain it (daily-case RMSE changes by under 0.2% as `gamma` ranges from 0.95 to 0.998 in exploratory fits), so letting it float just drives it to whichever bound is given rather than reflecting a real signal in the data. The initial exposed compartment `E(0)` is set analytically from the observed day-0 case count rather than fitted. All other parameters are held at the literature-informed or assumed values in the manuscript's Table 2.

For the cross-country transfer, only `alpha`, `beta`, and `kappa` are re-fitted per country; the remaining 14 of 17 parameters (including Uganda's fitted `delta`) are held fixed. See the manuscript, Section 5.4, for the reasoning and caveats behind this exercise.
