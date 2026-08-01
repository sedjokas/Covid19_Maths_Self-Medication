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
- Calibration to Uganda: **R₀ = 2.23**, daily-case RMSE of 272 cases/day, close agreement with the observed 1,448 cumulative deaths. Under an assumed 5% case-ascertainment fraction, roughly **95% of infections went uncaptured** by surveillance, with a peak home-care-to-hospital occupancy ratio of **27.7:1**.
- A three-parameter cross-country transfer (holding 14 of 17 Uganda-fitted parameters fixed) reproduces cumulative-mortality curves in Mozambique (R² = 0.964, MAPE 4.5%), Senegal (R² = 0.958, MAPE 3.8%), and Cameroon (R² = 0.966, MAPE 3.5%) — mean R² = 0.963.
- Community-surveillance and integrated-policy scenarios reduce simulated R₀ by up to **70.4%** and simulated mortality by up to **22.9%**.

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

where P_H = γω₂ + (1−γ)θ₂ and P_M = γθ₁ + (1−γ)ω₁ capture transmission pathways through hospitalized and home-care patients respectively, and ω₁, ω₂ are the total exit rates from H and M (Λ, μ₁, ε₁, φ₁, π₁, θ₁ and the M-analogues).

### Baseline Parameters (Table 2 of the manuscript)

| Parameter | Meaning | Value |
|-----------|---------|-------|
| Λ | Recruitment rate | 0.304 |
| α | Contact rate | 0.035 |
| β | Progression rate E→I | 0.095 |
| γ | Proportion of infected hospitalized | 0.29 |
| κ | Progression rate I→care states | 0.100 |
| θ₁ | Hospital-to-home discharge rate | 0.396 |
| θ₂ | Home-to-hospital transfer rate | 0.0396 |
| π₁ / π₂ | COVID-19 mortality (H / M) | 0.0196 / 0.009 |
| φ₁ / φ₂ | Recovery rate (H / M) | 0.900 / 6.25×10⁻⁵ |
| μ₁ / μ₂ | Natural / disease-induced mortality | 0.00576 / 0.196 |
| ε₁ / ε₂ | Return-to-susceptible rate (H / M) | 0.001 / 0.001 |
| η_M | Relative infectiousness of home-care | [0.50, 0.80] (baseline 0.50) |
| ξ_H | Relative infectiousness of hospitalized | [0.05, 0.20] (baseline 0.10) |

Under this baseline, R₀ = 0.2102.

---

## Mathematical Analysis

### Direction of Bifurcation at R₀ = 1 (Proposition, Section 3 of the manuscript)

The center-manifold calculation fully accounts for the state-dependence of the total-population denominator N_L in the frequency-dependent force of infection. The resulting quadratic center-manifold coefficient **a** factors as the product of two strictly positive quantities for every admissible (non-negative) parameter combination, so **a is strictly negative unconditionally** — independent of ε₁, ε₂, φ₁, φ₂, θ₁, θ₂. This means:

- The model undergoes **only forward (transcritical) bifurcation** at R₀ = 1.
- **R₀ < 1 is both necessary and sufficient for elimination.** There is no bistability window and no dependence on initial outbreak size.

This is verified three independent ways: (i) a closed-form, parameter-free sign argument; (ii) direct symbolic differentiation with no simplifying assumptions; (iii) numerical eigendecomposition of the Jacobian across a sweep of φ₁ ∈ [0.02, 0.95] and ε₁ = ε₂ ∈ {0, 0.001, 0.02, 0.10, 0.30}. It is additionally confirmed by direct (non-linearized) integration of the full nonlinear system from five initial conditions, from a tiny seed to a very large initial outbreak, at R₀ ≈ 0.988 — all five converge to the disease-free equilibrium. See `bifurcation_verification.py`.

**A numerical caveat that matters for interpreting simulations near R₀ ≈ 1:** the relaxation eigenvalue near the bifurcation point is extremely small (≈ −2.7×10⁻⁴ / day here, a ~3,700-day relaxation time). Simulations truncated at a few hundred days can show trajectories that have simply not yet reached their unique common equilibrium, which is easy to mistake for convergence to two distinct equilibria — `bifurcation_verification.py` therefore integrates Scenario 5 out to 60,000 days rather than a few hundred.

### Global Stability

- **Theorem (GAS of the DFE):** If R₀ ≤ 1, the disease-free equilibrium is globally asymptotically stable in Ω (Lyapunov function + LaSalle's invariance principle).
- **Theorem (GAS of the endemic equilibrium):** If R₀ > 1, a unique endemic equilibrium exists and is globally asymptotically stable in Ω \ {X₀} (Volterra–Goh-type Lyapunov function).

Both are proved on the six-dimensional living subsystem (S, E, I, H, M, R); D is excluded from the equilibrium notion and reported only as a derived, path-dependent output.

---

## Simulation Scenarios

The paper analyzes 10 scenarios using 4th-order Runge–Kutta numerical integration:

| # | Scenario | Key variable | Main finding |
|---|----------|-------------|-------------|
| 1 | Baseline dynamics | All compartments | 18:1 home-care/hospital ratio at peak; R₀ = 0.2102 |
| 2 | Care-seeking behavior | γ (hospitalization proportion) | Mortality ranked High Trust < Moderate < Self-medication |
| 3 | Contact-rate sweep | α | Raising α 7.1-fold crosses R₀ = 1; mortality rises 26.9-fold |
| 4 | Bifurcation verification | φ₁ ∈ [0.02, 0.95] | Coefficient a stays strictly negative throughout |
| 5 | No-bistability demonstration | 5 initial conditions, 60,000-day horizon | All trajectories converge to the same disease-free equilibrium |
| 6 | Hospital discharge policy | θ₁ (early discharge under capacity pressure) | Aggressive discharge raises mortality ~7.3% |
| 7 | Joint sensitivity of mortality | γ, η_M | Mortality contour spans ~3.8–5.0 deaths/1,000 |
| 8 | Community surveillance | η_M, θ₂ | Intensive vs. none: R₀ falls 70.4%, deaths fall 33.3% |
| 9 | Integrated policy comparison | γ, η_M, θ₂ combined | Integrated approach: 22.9% mortality reduction |
| 10 | Convergence to the endemic equilibrium | Phase-plane, R₀ > 1 | Unique endemic attractor from all tested initial conditions |

---

## Calibration and Validation

### Uganda Calibration (Delta wave, June–October 2021)
- Population: 47,123,531
- Calibrated R₀ = **2.23**
- Cumulative deaths reproduced: observed 1,448 (by construction of the fitted mortality-scaling factor δ ≈ 1.30)
- RMSE = 272 cases/day on daily incidence (in-sample)
- Under an assumed, literature-informed (not independently estimated) 5% case-ascertainment fraction: **95% of infections implied undetected**
- Peak home-care-to-hospital occupancy ratio: **27.7:1**

These case-ascertainment and occupancy-ratio figures are model-implied quantities conditional on the assumed ρ = 0.05, not independently measured facts — see the manuscript's Limitations section.

### Cross-Country Parameter Transfer

Three of seventeen parameters (α, β, κ) are re-fitted per country; the remaining fourteen are held at their Uganda-fitted values.

| Country | Window | R² | MAPE |
|---------|--------|----|------|
| Mozambique | Jul–Nov 2020 | 0.9639 | 4.5% |
| Senegal | Jul–Oct 2020 | 0.9582 | 3.8% |
| Cameroon | Jun–Sep 2020 | 0.9660 | 3.5% |
| **Mean** | | **0.963** | |

This is reported as evidence against gross incompatibility of the fitted disease-progression parameters across settings — a necessary but not sufficient condition for transferability — not as independent external validation. See the manuscript's Limitations section for the full discussion.

---

## Repository Structure

```
./
├── README.md
├── LICENSE
├── covid_simulation_white_bg.py       # Main SEIHM-R-D simulation (Scenarios 1-3, 6-9)
├── bifurcation_verification.py        # Bifurcation verification (Scenarios 4 & 5): center-manifold
│                                       #   coefficient sweep + long-horizon no-bistability demo
├── export_gamma_panels.py             # Uganda calibration panels & gamma sensitivity
│
├── data/
│   ├── baseline_simulation_data.csv       # Time series - all 7 compartments, baseline run
│   ├── endemic_comparison_data.csv        # Baseline vs. endemic equilibrium comparison
│   ├── summary_table.csv                  # R0 and mortality summary across scenarios
│   ├── bifurcation_coefficient_a.csv      # Coefficient a across the (phi_1, epsilon) sweep
│   └── no_bistability_trajectories.csv    # 5 initial conditions x 60,000-day integration
│
└── figures/
    ├── fig_bifurcation_coefficient_a.png      # Scenario 4: coefficient a is negative throughout
    ├── fig_no_bistability.png         # Scenario 5: convergence to a common DFE
    └── ... (per-scenario figures, see filenames)
```

---

## Requirements

```bash
pip install numpy scipy matplotlib pandas
```

Python ≥ 3.8. Simulations use `scipy.integrate.odeint` (an adaptive-step LSODA-based solver); the manuscript itself reports results from a fixed-step 4th-order Runge–Kutta integration with step-size convergence checked by halving h.

## Usage

```bash
# Main simulation - policy scenarios (1-3, 6-9)
python covid_simulation_white_bg.py

# Bifurcation verification - coefficient-a sweep and long-horizon no-bistability demo (4 & 5)
python bifurcation_verification.py

# Uganda calibration panels (gamma sensitivity, validation figures)
python export_gamma_panels.py
```

Each script generates publication-ready figures at 300 DPI (white background).

---

## Key Policy Implications

1. **R₀ < 1 is both necessary and sufficient for elimination in this model.** There is no bistable "trap" to avoid: driving R₀ below 1, by whatever mix of care-seeking, infection-control, or surveillance measures is feasible, suffices.
2. **The home-care compartment is the dominant transmission driver** under baseline and Uganda-calibrated parameters alike, both because more infections are routed there and because it clears more slowly.
3. **Community surveillance** (reducing η_M, increasing θ₂) is the most effective single intervention tested, reducing simulated R₀ by up to 70.4%.
4. **Integrated policies** (surveillance + awareness) yield the largest simulated mortality reduction (22.9%) among the scenarios tested.
5. **Parameter transferability** across Sub-Saharan Africa (mean R² = 0.963 with 14 of 17 parameters fixed from Uganda) is consistent with, though does not prove, shared disease-progression structure across the four countries studied.

All of the above are scenario-conditional, assumption-driven projections from a model that has not yet been benchmarked against simpler baseline models or validated on temporally held-out data — see the manuscript's Limitations section (Section 6) before treating any specific percentage as a validated causal intervention effect.

---

## Citation

> Kasereka S.K., Ngoie R.-B.M., Doungmo E.F.G., Mafuta E.M., Chedjou J.C., Kabengele E.M., Kyamakya K. Modeling COVID-19 Spread in Low-Income Countries: A SEIHM-R-D Framework with Bidirectional Hospital and Home-Care Transitions. Prepared for submission to *Infectious Disease Modelling*.

**Corresponding author:** Selain K. Kasereka — selain.kasereka@aau.at

---

## License

This project is licensed under the [MIT License](LICENSE).
