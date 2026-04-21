# Understanding COVID-19 Spread in Low-Income Countries
## A SEIHM-R-D Modeling Approach Including Hospitalized and Home-Care Populations

**Authors:** Selain K. Kasereka · Ruffin-Benoît M. Ngoie · Emile Franc G. Doungmo · Eric M. Mafuta · Jean Chamberlain Chedjou · Emmanuel M. Kabengele · Kyandoghere Kyamakya

**Affiliations:** University of Klagenfurt (Austria) · ABIL-LAB Kinshasa (DRC) · Institut Supérieur Pédagogique de Mbanza-Ngungu (DRC) · University of South Africa · University of Kinshasa · University of Geneva

**Journal:** Journal of Infection and Public Health (Elsevier)

**Status:** Preprint submitted — April 2026 · DOI: not yet assigned

---

## Abstract

In low-income countries, many COVID-19 patients avoided hospitalisation due to financial constraints and limited hospital capacity, opting for home-care instead. Standard epidemic models rarely capture this behaviour, leading to systematic underestimation of transmission. This work develops a **SEIHM-R-D** model explicitly distinguishing hospitalised (*H*) from home-care (*M*) patients, with differential infectiousness and bidirectional transitions between care settings.

Key results:
- Mathematical analysis establishes **backward bifurcation**: reducing R₀ below 1 may be insufficient for elimination when hospital recovery is poor
- Model calibrated on Uganda's COVID-19 Delta-wave data (June–October 2021): **R₀ = 2.23**, MAPE < 0.1%, with **95% of infections remaining undetected** and a 27.7:1 peak home-care-to-hospital ratio
- External validation on Mozambique, Senegal, and Cameroon: mean **R² = 0.963** with 82% of parameters fixed from Uganda calibration
- Community surveillance reduces R₀ by **70.4%** and prevents **33% of deaths**; integrated policies yield **22.9% mortality reduction**

---

## The SEIHM-R-D Model

### Compartments

| Symbol | Description |
|--------|-------------|
| **S** | Susceptible |
| **E** | Exposed (latent period) |
| **I** | Infectious (unmanaged, community) |
| **H** | Hospitalised |
| **M** | Home-care / self-medication |
| **R** | Recovered |
| **D** | Cumulative deaths |

The living population is N_L = S + E + I + H + M + R; D is tracked separately.

### Force of Infection

$$\lambda = \frac{\alpha \, S \, (I + \eta_M M + \xi_H H)}{N_L}$$

where η_M ∈ [0.50, 0.80] is the relative infectiousness of home-care patients and ξ_H ∈ [0.05, 0.20] for hospitalised patients — reflecting reduced but non-zero infection control in household settings.

### Basic Reproduction Number R₀

Derived via the next-generation matrix method (Eq. 14):

$$\mathcal{R}_0 = \frac{\alpha\beta\,[\mathcal{D} + \kappa\xi_H P_H + \kappa\eta_M P_M]}{(\mu_1+\beta)(\mu_1+\mu_2+\kappa)(\omega_1\omega_2 - \theta_1\theta_2)}$$

where P_H = γω₂ + (1−γ)θ₂ and P_M = γθ₁ + (1−γ)ω₁ capture transmission pathways through hospitalised and home-care patients respectively.

### Baseline Parameters (Table 2 of the manuscript)

| Parameter | Meaning | Value |
|-----------|---------|-------|
| Λ | Recruitment rate | 0.304 |
| α | Contact rate | 0.035 |
| β | Progression rate E→I | 0.095 |
| γ | Proportion of infected hospitalised | 0.29 |
| κ | Progression rate I→care states | 0.100 |
| θ₁ | Hospital-to-home discharge rate | 0.396 |
| θ₂ | Home-to-hospital transfer rate | 0.0396 |
| π₁ / π₂ | COVID-19 mortality (H / M) | 0.0196 / 0.009 |
| φ₁ / φ₂ | Recovery rate (H / M) | 0.900 / 6.25×10⁻⁵ |
| μ₁ / μ₂ | Natural / disease-induced mortality | 0.00576 / 0.196 |
| η_M | Relative infectiousness of home-care | [0.50, 0.80] |
| ξ_H | Relative infectiousness of hospitalised | [0.05, 0.20] |

---

## Mathematical Analysis

### Backward Bifurcation (Theorem 3.1)

The model undergoes **backward (subcritical) bifurcation** at R₀ = 1 when the hospital recovery rate φ₁ falls below a critical threshold. In this regime, two endemic equilibria coexist for R̃₀ < R₀ < 1, while the disease-free equilibrium remains locally stable — **bistability**. Consequently, reducing R₀ below 1 is necessary but not sufficient for disease elimination; the actual threshold is R̃₀ < 1.

The bifurcation condition depends on α̃:
- If α < α̃ → **backward bifurcation** (bistability, initial conditions determine outcome)
- If α ≥ α̃ → **forward bifurcation** (classical threshold behaviour)

This is a critical policy implication: in settings with poor hospital recovery, large outbreaks can persist even when R₀ < 1.

### Global Stability

- **Theorem 3.3** (GAS of DFE): If R₀ ≤ 1 under the forward-bifurcation condition, the disease-free equilibrium is globally asymptotically stable (Lyapunov function, LaSalle's principle).
- **Theorem 3.4** (GAS of endemic equilibrium): If R₀ > 1, the endemic equilibrium E* is globally asymptotically stable in Ω \ {E₀} (Volterra–Goh function).

---

## Simulation Scenarios

The paper analyses 10 scenarios using RK4 numerical integration:

| # | Scenario | Key variable | Main finding |
|---|----------|-------------|-------------|
| 1 | **Baseline epidemic dynamics** | All compartments | 27.7:1 home-care/hospital ratio at peak |
| 2 | **Impact of self-medication** | γ (hospitalisation rate) | Higher γ → lower transmission and mortality |
| 3 | **R₀ < 1 to R₀ > 1 transition** | α (contact rate) | Loss of epidemic control and endemicity |
| 4 | **Hospital recovery bifurcation** | φ₁ (hospital recovery) | φ₁ governs forward vs. backward bifurcation |
| 5 | **Backward bifurcation & bistability** | Initial conditions | Identical R₀, deaths ranging 0.5–351.2 per 1,000 |
| 6 | **Hospital capacity constraints** | θ₁ (early discharge) | Overcrowding forces patients home, raises mortality |
| 7 | **Epidemic mortality dynamics** | π₁, π₂, η_M | Home-care dominates transmission and deaths |
| 8 | **Community surveillance** | η_M, θ₂ | R₀ reduced by 70.4%; 33% fewer deaths |
| 9 | **Integrated policy packages** | Combined parameters | 22.9% mortality reduction vs. status quo |
| 10 | **Phase-plane bifurcation** | All 7 state variables | Hysteresis and basin of attraction in bistability region |

---

## Calibration and Validation

### Uganda Calibration (June–October 2021)
- Population: 47,123,531
- Calibrated R₀ = **2.23**
- Cumulative deaths reproduced: observed 1,448 — **MAPE < 0.1%**
- RMSE = 272 cases/day on daily incidence
- **95% of infections remained undetected** (reporting rate ρ = 5%)
- Peak home-care-to-hospital ratio: **27.7:1**

### External Validation
| Country | R² | Parameters fixed from Uganda |
|---------|----|------------------------------|
| Mozambique | ~0.963 (mean) | 82% (14 of 17) |
| Senegal | ~0.963 (mean) | 82% (14 of 17) |
| Cameroon | ~0.963 (mean) | 82% (14 of 17) |

---

## Repository Structure

```
./
├── README.md
├── LICENSE
├── covid_simulation_white_bg.py       # Main SEIHM-R-D simulation (Scenarios 1–9)
├── backward_bifurcation_demo.py       # Backward bifurcation analysis (Scenarios 4 & 5)
├── export_gamma_panels.py             # Uganda calibration panels & γ sensitivity
│
├── data/
│   ├── baseline_simulation_data.csv   # Time series – all 7 compartments, baseline run
│   ├── bistability_data.csv           # Bistability – 5 initial conditions (Scenario 5)
│   ├── endemic_comparison_data.csv    # Baseline vs. endemic equilibrium comparison
│   └── summary_table.csv             # R₀ and mortality summary across scenarios
│
└── figures/
    ├── Scenario 1  – Plot_Covid1.png / Plot_Covid3.png / Plot_Covid5.png / Plot_Covid7.png
    │                 IandD1.jpeg / IandD2.jpeg / IandD3.jpeg
    │                 SimGen1_All.png / SimGen1_Infected.png / SimGen1_Deaths.png
    │                 SimGen1_Living.png / SimGen1_PeakMH.png
    ├── Scenario 2  – gamma1.png / gamma2.png / gamma3.png / gamma4.png
    │                 all150_alpha025.png / inf150_alpha025.png
    ├── Scenario 3  – S_endemic.png / I_endemic.png / H_endemic.png
    │                 M_endemic.png / D_endemic.png / Peak_endemic.png
    ├── Scenarios 4 & 5 – forward1.png / backward1.png / forward_backward1.png
    │                     Sbistab.png / Ebistab.png / Ebistab1.png / Ibistab.png
    │                     HCbistab.png / Dbistab.png
    ├── Scenario 6  – capa1.png / Capa2.png / Capa3.png / Capa4.png
    ├── Scenario 7  – Morts.png / MortH.png / PeakMH.png
    │                 Death_policy.png / Death_policy2.png
    ├── Scenario 8  – EIMPlot01.png / EIMPlot03.png / EIMPlot06.png / EIMPlot09.png
    ├── Scenario 9  – Inf_policy.png / CumulDeathsCI.png / Dealycase95CI.png
    └── Scenario 10 – scenario10_panel_a_SI.png / scenario10_panel_b_SM.png
                      scenario10_panel_c_SH.png / scenario10_panel_d_IM.png
                      scenario10_panel_e_IH.png / scenario10_panel_f_HM.png
                      Bassin.png
```

---

## Requirements

```bash
pip install numpy scipy matplotlib pandas
```

Python ≥ 3.8. Simulations use the classical 4th-order Runge-Kutta method (via `scipy.integrate.odeint`).

## Usage

```bash
# Main simulation — all 9 policy scenarios
python covid_simulation_white_bg.py

# Backward bifurcation analysis — forward vs. backward comparison
python backward_bifurcation_demo.py

# Uganda calibration panels (γ sensitivity, validation figures)
python export_gamma_panels.py
```

Each script generates publication-ready figures at 300 DPI (white background, dark color palette).

---

## Key Policy Implications

1. **R₀ < 1 is not sufficient for elimination** when hospital recovery is poor — backward bifurcation creates a bistability zone where large outbreaks persist even with R₀ < 1.
2. **The home-care compartment is the dominant transmission driver** — near-zero home-care recovery (φ₂ ≈ 6.25×10⁻⁵) creates a persistent endemic reservoir.
3. **Community surveillance** (reducing η_M and increasing θ₂) is the most effective single intervention, reducing R₀ by 70.4%.
4. **Integrated policies** (surveillance + awareness + capacity) yield the greatest mortality reduction (22.9%).
5. **Parameter transferability** across Sub-Saharan Africa (R² = 0.963 with 82% parameters fixed) confirms the model's applicability for regional public health planning.

---

## Citation

> Kasereka S.K., Ngoie R.-B.M., Doungmo E.F.G., Mafuta E.M., Chedjou J.C., Kabengele E.M., Kyamakya K. (2026). *Understanding COVID-19 Spread in Low-Income Countries: A Modeling Approach Including Hospitalized and Home-Care Populations.* Preprint submitted to Elsevier.

**Corresponding author:** Selain K. Kasereka — selain.kasereka@aau.at

---

## License

This project is licensed under the [MIT License](LICENSE).
