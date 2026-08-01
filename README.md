# Modeling COVID-19 Spread in Low-Income Countries

## A SEIHM-R-D Framework with Bidirectional Hospital and Home-Care Transitions

**Authors:** Selain K. Kasereka · Ruffin-Benoît M. Ngoie · Emile Franc G. Doungmo · Eric M. Mafuta · Jean Chamberlain Chedjou · Emmanuel M. Kabengele · Kyandoghere Kyamakya

**Affiliations:** University of Klagenfurt (Austria) · ABIL-LAB Kinshasa (DRC) · Institut Supérieur Pédagogique de Mbanza-Ngungu (DRC) · University of South Africa · University of Kinshasa · University of Geneva

**Target journal:** Infectious Disease Modelling (KeAi / Elsevier)

**Status:** Manuscript prepared for submission. DOI: not yet assigned.

---

## Abstract

Most COVID-19 transmission models leave out people who manage illness at home rather than in a hospital. This pattern was common in low-income countries, driven by financial barriers, limited hospital capacity, and local care-seeking habits. We build a SEIHM-R-D compartmental model that separates hospitalized (H) from home-care (M) patients. Patients can move in both directions between the two care settings, and each setting carries its own level of infectiousness.

We establish the equilibria and the local and global stability of the six-dimensional living subsystem (S, E, I, H, M, R), with cumulative deaths D treated as a monotone output rather than as an equilibrium coordinate. We derive the basic reproduction number R₀ via the next-generation-matrix method, and establish the direction of bifurcation at R₀ = 1 analytically (center-manifold method) and numerically. The model is calibrated to Uganda's COVID-19 Delta-wave surveillance data (June–October 2021), and a subset of its parameters is transferred to Mozambique, Senegal, and Cameroon.

**Key results:**
- We prove the center-manifold coefficient is **strictly negative for every admissible parameter combination**: the model undergoes **only forward (transcritical) bifurcation** at R₀ = 1. **R₀ < 1 is both necessary and sufficient for elimination — there is no bistability regime.**
- Calibration to Uganda's real Delta-wave surveillance data: **R₀ = 1.97** at the start of the wave, daily-case RMSE of 44.6 cases/day, close agreement with the observed 2,853 cumulative deaths. Under an assumed 5% case-ascertainment fraction, roughly **95% of infections went uncaptured** by surveillance, with a peak home-care-to-hospital occupancy ratio of **9.3:1**.
- A three-parameter cross-country transfer (holding 14 of 17 Uganda-fitted parameters fixed) reproduces cumulative-mortality curves in Mozambique (R² = 0.990, MAPE 19.1%), Senegal (R² = 0.997, MAPE 6.0%), and Cameroon (R² = 0.981, MAPE 14.8%) — mean R² = 0.989.
- Community-surveillance and integrated-policy scenarios reduce simulated R₀ by up to **36.7%** and simulated mortality by up to **15.1%**.

---

## The SEIHM-R-D Model

### Compartments

| Symbol | Description |
|--------|-------------|
| **S** | Susceptible |
| **E** | Exposed (latent period) |
| **I** | Infectious (not yet in a care pathway) |
| **H** | Hospitalized |
| **M** | Home-care / self-medication |
| **R** | Recovered |
| **D** | Cumulative deaths (monotone output, not an equilibrium coordinate) |

The living population is N_L = S + E + I + H + M + R; D is tracked separately as a derived, non-decreasing output that has no influence on, and is not influenced by, the living subsystem.

### Force of Infection

$$\lambda = \frac{\alpha \, (I + \eta_M M + \xi_H H)}{N_L}$$

where η_M ∈ [0.50, 0.80] is the relative infectiousness of home-care patients and ξ_H ∈ [0.05, 0.20] for hospitalized patients, reflecting reduced but non-zero infection control in each setting. λ multiplies S separately in the S and E equations below.

### Basic Reproduction Number R₀

Derived via the next-generation-matrix method:

$$\mathcal{R}_0 = \frac{\alpha\beta\,[\mathcal{D} + \kappa\xi_H P_H + \kappa\eta_M P_M]}{(\mu_1+\beta)(\mu_1+\mu_2+\kappa)(\omega_1\omega_2 - \theta_1\theta_2)}$$

where P_H = γω₂ + (1−γ)θ₂ and P_M = γθ₁ + (1−γ)ω₁ capture transmission pathways through hospitalized and home-care patients respectively, and ω₁, ω₂ are the total exit rates from H and M (μ₁, ε₁, φ₁, π₁, θ₁ and the M-analogues).

### Baseline Parameters (Table 2 of the manuscript)

Used for the illustrative Section-4 scenarios below. Section 5's Uganda calibration uses a separately fitted parameter set (see that section).

| Parameter | Meaning | Value | Source |
|-----------|---------|-------|--------|
| Λ | Recruitment rate | 0.304 | Ndondo et al. 2021 |
| α | Contact rate | 0.035 | assumed |
| β | Progression rate E→I | 0.1961 | McAloon 2020; Lauer 2020 (incubation period) |
| γ | Proportion hospitalized | 0.20 | WHO-China Joint Mission 2020 |
| κ | Progression rate I→care states | 0.100 | calibrated |
| θ₁ | Hospital-to-home rate | 0.396 | assumed |
| θ₂ | Home-to-hospital rate | 0.0396 | assumed |
| π₁ / π₂ | COVID-19 mortality (H / M) | 0.0196 / 0.009 | assumed, order-of-magnitude consistent with CFR literature |
| φ₁ / φ₂ | Recovery rate (H / M) | 0.15 / 0.0714 | Alimohamadi 2022 (LOS); Tenforde 2020 (WHO recovery guidance) |
| μ₁ / μ₂ | Natural / disease-induced mortality | 0.00576 / 0.0008 | Sinan et al. 2021; Verity et al. 2020 |
| ε₁ / ε₂ | Return-to-susceptible rate (H / M) | 0.001 / 0.001 | assumed |
| η_M | Relative infectiousness of home-care | [0.50, 0.80] (baseline 0.50) | assumed |
| ξ_H | Relative infectiousness of hospitalized | [0.05, 0.20] (baseline 0.10) | assumed |

Under this baseline, R₀ = 0.4729.

Entries marked "assumed" have no direct literature analogue — the SEIHM-R-D care-setting split is, to our knowledge, novel — and are chosen to be broadly plausible rather than fitted. See the manuscript's Table 2 footnotes and Limitations section for the full discussion of which parameters are literature-derived versus assumed, and why some literature values were not used (e.g., where they broke the model's intended qualitative behavior in a model-specific, non-portable way).

---

## Mathematical Analysis

### Direction of Bifurcation at R₀ = 1 (Proposition, Section 3 of the manuscript)

The center-manifold calculation fully accounts for the state-dependence of the total-population denominator N_L in the frequency-dependent force of infection. The resulting quadratic center-manifold coefficient **a** factors as the product of two strictly positive quantities for every admissible (non-negative) parameter combination, so **a is strictly negative unconditionally** — independent of ε₁, ε₂, φ₁, φ₂, θ₁, θ₂. This means:

- The model undergoes **only forward (transcritical) bifurcation** at R₀ = 1.
- **R₀ < 1 is both necessary and sufficient for elimination.** There is no bistability window and no dependence on initial outbreak size.

This is verified three independent ways: (i) a closed-form, parameter-free sign argument; (ii) direct symbolic differentiation with no simplifying assumptions; (iii) numerical eigendecomposition of the Jacobian across a sweep of φ₁ ∈ [0.02, 0.95] and ε₁ = ε₂ ∈ {0, 0.001, 0.02, 0.10, 0.30}. It is additionally confirmed by direct (non-linearized) integration of the full nonlinear system from five initial conditions, from a tiny seed to a very large initial outbreak, at R₀ ≈ 0.988 — all five converge to the disease-free equilibrium. See `codes/bifurcation_verification.py`.

**A numerical caveat that matters for interpreting simulations near R₀ ≈ 1:** the relaxation eigenvalue near the bifurcation point is extremely small (≈ −6.7×10⁻⁴ / day at φ₁=0.15, α=0.0731, a ~1,480-day relaxation time). Simulations truncated at a few hundred days can show trajectories that have simply not yet reached their unique common equilibrium, which is easy to mistake for convergence to two distinct equilibria — `codes/bifurcation_verification.py` therefore integrates Scenario 5 out to 60,000 days rather than a few hundred.

### Global Stability

- **Theorem (GAS of the DFE):** If R₀ ≤ 1, the disease-free equilibrium is globally asymptotically stable in Ω (Lyapunov function + LaSalle's invariance principle).
- **Theorem (GAS of the endemic equilibrium):** If R₀ > 1, a unique endemic equilibrium exists and is globally asymptotically stable in Ω \ {X₀} (Volterra–Goh-type Lyapunov function).

Both are proved on the six-dimensional living subsystem (S, E, I, H, M, R); D is excluded from the equilibrium notion and reported only as a derived, path-dependent output.

---

## Simulation Scenarios

The paper analyzes 10 scenarios (`codes/model.py` + `codes/run_all.py`), all under the Table 2 baseline (R₀ = 0.4729) unless a scenario explicitly varies a parameter:

| # | Scenario | Key variable | Main finding |
|---|----------|-------------|-------------|
| 1 | Baseline dynamics | All compartments | Transient outbreak with a 9.3:1 peak M:H ratio; cumulative mortality 0.97/1,000 by day 150 |
| 2 | Care-seeking behavior | γ (hospitalization proportion) | Mortality ranked High Trust (0.93/1,000) < Moderate (0.95) < Self-medication (0.98) |
| 3 | Contact-rate sweep | α | Raising α 8-fold (0.037→0.296) crosses R₀=1 at α≈0.074; mortality rises 81.5-fold |
| 4 | Bifurcation verification | φ₁ ∈ [0.02, 0.95] | Coefficient a stays strictly negative throughout |
| 5 | No-bistability demonstration | 5 initial conditions, 60,000-day horizon | All trajectories converge to the same disease-free equilibrium (R₀≈0.988) |
| 6 | Hospital discharge policy | θ₁ (discharge rate under capacity pressure) | Aggressive discharge raises mortality ~4.6% |
| 7 | Joint sensitivity of mortality | γ, η_M | Mortality more sensitive to η_M than to γ under Table 2 |
| 8 | Community surveillance | η_M, θ₂ | Intensive vs. none: R₀ falls 36.7%, deaths fall 30.5% |
| 9 | Integrated policy comparison | γ, η_M, θ₂ combined | Integrated approach: 15.1% mortality reduction, 22.7% R₀ reduction |
| 10 | Convergence to the endemic equilibrium | S-I phase plane, R₀=2.5 | Unique endemic attractor from all tested initial conditions |

---

## Calibration and Validation

### Uganda Calibration (Delta wave, June–October 2021)

- Population: 47,123,531 (UBOS 2021 projection)
- Real Delta-wave surveillance totals (OWID, cross-checked against an independently archived Nov-2021 snapshot): 78,410 confirmed cases, 2,853 confirmed deaths
- Fitted contact rate: logistic decline from α_hi=0.323 to α_lo=0.100, switch centered at day ≈17.5 (17-18 June 2021 — closely matching Uganda's actual nationwide lockdown, announced 18 June 2021)
- Calibrated R₀ = **1.97** at the start of the wave, declining to R_eff = **0.61**
- Cumulative deaths reproduced: observed 2,853 (by construction of the fitted mortality-scaling factor δ ≈ 0.0222)
- RMSE = 44.6 cases/day on daily incidence (in-sample; 7 days in late August excluded as a reporting-backfill artifact, see below)
- Under an assumed, literature-informed (not independently estimated) 5% case-ascertainment fraction: **95% of infections implied undetected**
- Peak home-care-to-hospital occupancy ratio: **9.3:1**

These case-ascertainment and occupancy-ratio figures are model-implied quantities conditional on the assumed ρ = 0.05, not independently measured facts — see the manuscript's Limitations section.

**Data note:** OWID's 7-day smoothed case series for Uganda shows an isolated ≈3,100 cases/day plateau for exactly 23–29 August 2021, collapsing to ≈150/day the next day. This is a reporting-backfill artifact (zero reported cases 16–21 August, followed by a single-day dump of 21,653 cases on 22 August, smeared by the rolling average), not a genuine second epidemic peak. These seven days are excluded from the calibration objective; see `codes/calibrate_uganda.py`.

### Cross-Country Parameter Transfer

Three of seventeen parameters (α, β, κ) are re-fitted per country against real OWID cumulative-mortality curves; the remaining fourteen (including Uganda's fitted δ) are held at their Uganda values.

| Country | Window | R² | MAPE |
|---------|--------|----|------|
| Mozambique | Jul–Nov 2020 | 0.990 | 19.1% |
| Senegal | Jul–Oct 2020 | 0.997 | 6.0% |
| Cameroon | Jun–Sep 2020 | 0.981 | 14.8% |
| **Mean** | | **0.989** | |

MAPE is elevated, especially for Mozambique, because cumulative deaths are small (single/low double digits) early in each window, so a handful of deaths' absolute error is a large percentage error. This is reported as evidence against gross incompatibility of the fitted disease-progression parameters across settings — a necessary but not sufficient condition for transferability — not as independent external validation. See the manuscript's Limitations section for the full discussion.

---

## Repository Structure

```
./
├── README.md
├── LICENSE
├── codes/                                  # All Python code and the data it uses
│   ├── README.md                             # Code-specific documentation
│   ├── model.py                               # Core SEIHM-R-D model (RHS, R0, center-manifold a)
│   ├── run_all.py                              # Section 4 scenarios (1-10)
│   ├── bifurcation_verification.py             # Independent check of Scenarios 4 & 5
│   ├── calibrate_uganda.py                      # Uganda Delta-wave calibration
│   ├── finalize_uganda.py                        # Goodness-of-fit, 95% CIs, gamma sensitivity, figures
│   ├── transfer_countries.py                      # Cross-country parameter transfer
│   ├── data/                                       # Cleaned data used directly by the scripts above
│   └── results/                                     # JSON outputs (fitted parameters, goodness-of-fit)
│
└── figures/                                # The manuscript's actual figures (as included in the PDF)
```

`codes/data/` holds the cleaned Uganda and transfer-country case/death series actually used for fitting, plus the bifurcation-verification script's own output data. It does not include the full multi-country OWID dump those series were extracted from (tens of megabytes); `codes/README.md` documents how to obtain it.

---

## Requirements

```bash
pip install numpy scipy matplotlib scikit-learn
```

Python ≥ 3.8. Simulations use `scipy.integrate.odeint` (an adaptive-step LSODA-based solver); the manuscript itself reports results from a fixed-step 4th-order Runge–Kutta integration with step-size convergence checked by halving h. Calibration uses `scipy.optimize.differential_evolution`.

## Usage

```bash
cd codes

# Section 4: baseline scenarios (1-10)
python run_all.py

# Section 3: independent bifurcation check (Scenarios 4 & 5)
python bifurcation_verification.py

# Section 5: Uganda calibration, then goodness-of-fit/CIs/figures, then cross-country transfer
python calibrate_uganda.py
python finalize_uganda.py
python transfer_countries.py
```

See `codes/README.md` for the full run order, data provenance, and notes on which parameters are fitted versus held fixed at each stage.

---

## Key Policy Implications

1. **R₀ < 1 is both necessary and sufficient for elimination in this model.** There is no bistable "trap" to avoid: driving R₀ below 1, by whatever mix of care-seeking, infection-control, or surveillance measures is feasible, suffices.
2. **The home-care compartment is the dominant transmission driver** under baseline and Uganda-calibrated parameters alike, both because more infections are routed there and because it clears more slowly.
3. **Community surveillance** (reducing η_M, increasing θ₂) is the most effective single intervention tested, reducing simulated R₀ by up to 36.7%.
4. **Integrated policies** (surveillance + awareness) yield the largest simulated mortality reduction (15.1%) among the scenarios tested.
5. **Parameter transferability** across Sub-Saharan Africa (mean R² = 0.989 with 14 of 17 parameters fixed from Uganda) is consistent with, though does not prove, shared disease-progression structure across the four countries studied.

All of the above are scenario-conditional, assumption-driven projections from a model that has not yet been benchmarked against simpler baseline models or validated on temporally held-out data — see the manuscript's Limitations section (Section 6) before treating any specific percentage as a validated causal intervention effect.

---

## Citation

> Kasereka S.K., Ngoie R.-B.M., Doungmo E.F.G., Mafuta E.M., Chedjou J.C., Kabengele E.M., Kyamakya K. Modeling COVID-19 Spread in Low-Income Countries: A SEIHM-R-D Framework with Bidirectional Hospital and Home-Care Transitions. Prepared for submission to *Infectious Disease Modelling*.

**Corresponding author:** Selain K. Kasereka — selain.kasereka@aau.at

---

## License

This project is licensed under the [MIT License](LICENSE).
