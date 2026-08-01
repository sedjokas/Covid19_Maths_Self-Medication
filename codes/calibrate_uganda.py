import numpy as np
import csv
from pathlib import Path
from scipy.integrate import odeint
from scipy.optimize import differential_evolution

HERE = Path(__file__).resolve().parent
DATA_PATH = HERE / 'data' / 'uganda_calib_data.csv'
RESULTS_PATH = HERE / 'results' / 'uganda_fit_results.json'
RESULTS_PATH.parent.mkdir(exist_ok=True)

# ---- Load real Uganda calibration data (OWID, archived Nov 2021 snapshot) ----
days, dates, cases_smoothed, cum_deaths = [], [], [], []
with open(DATA_PATH) as f:
    reader = csv.DictReader(f)
    for r in reader:
        days.append(int(r['day']))
        dates.append(r['date'])
        cases_smoothed.append(float(r['new_cases_smoothed']))
        cum_deaths.append(float(r['cum_deaths']))

days = np.array(days)
cases_smoothed = np.array(cases_smoothed)
cum_deaths = np.array(cum_deaths)
deaths_baseline = cum_deaths[0]
deaths_net = cum_deaths - deaths_baseline  # model D(0)=0, compare to net increase

# Days 83-89 (2021-08-23 to 2021-08-29) carry a reporting-backfill artifact:
# raw counts were 0 for 2021-08-16 to 2021-08-21, then a 21,653-case dump landed
# on 2021-08-22. The 7-day rolling average smears that single-day dump across
# the following week, producing a flat ~3,100/day plateau (days 83-89) that
# collapses back to ~150/day the instant the dump exits the averaging window
# (day 90). No genuine epidemic curve is flat for exactly 7 days then falls off
# a cliff in one day. These days are excluded from the fit; deaths reporting
# does not show a comparable artifact (cum_deaths rises smoothly throughout)
# so no exclusion is applied there.
ARTIFACT_DAYS = set(range(83, 90))
fit_mask = np.array([d not in ARTIFACT_DAYS for d in days])

POP = 47123531  # UBOS population projection (real people)

# ---- Fixed parameters from the reconciled Table 2 ----
FIXED = {
    'Lambda': 0.304, 'beta': 1/5.1,
    'theta_1': 0.396, 'theta_2': 0.0396, 'pi_1': 0.0196, 'pi_2': 0.009,
    'phi_1': 0.15, 'phi_2': 1/14, 'mu_1': 0.00576, 'mu_2': 0.0008,
    'epsilon_1': 0.001, 'epsilon_2': 0.001, 'eta_M': 0.50, 'xi_H': 0.10,
    'gamma': 0.20,  # fixed at the literature baseline (WHO-China Joint Mission
                     # severity split); case-count data alone barely identifies
                     # gamma (RMSE flat from gamma=0.95 to 0.998), so leaving it
                     # free just drives it to whichever bound is given
}
RHO = 0.05  # assumed case-ascertainment fraction (unchanged from manuscript)

# The model's own demographic disease-free equilibrium is N* = Lambda/mu_1,
# not an arbitrary round number. The simulation must be seeded there so total
# living population stays near steady state absent disease (otherwise a slow
# demographic relaxation, timescale 1/mu_1 ~ 174 days, would be conflated
# with the ~152-day epidemic signal). SCALE converts model units to real
# people, consistently, for both cases and deaths.
N_STAR = FIXED['Lambda'] / FIXED['mu_1']
SCALE = POP / N_STAR

# Anchor the initial exposed compartment to the observed day-0 case count
# (rather than an arbitrary guess), using the same case-reporting relation
# used throughout the fit: reported_cases = RHO * beta * E * SCALE.
E0 = cases_smoothed[0] / (RHO * FIXED['beta'] * SCALE)


def rhs(y, t, p):
    S, E, I, H, M, R, D = y
    N_L = S + E + I + H + M + R
    alpha_t = p['alpha_lo'] + (p['alpha_hi'] - p['alpha_lo']) / (1 + np.exp(p['k_sw'] * (t - p['t_sw'])))
    lam = alpha_t * (I + p['eta_M'] * M + p['xi_H'] * H) / N_L
    omega_1 = p['mu_1'] + p['epsilon_1'] + p['phi_1'] + p['pi_1'] + p['theta_1']
    omega_2 = p['mu_1'] + p['epsilon_2'] + p['phi_2'] + p['pi_2'] + p['theta_2']
    dS = p['Lambda'] - p['mu_1'] * S - lam * S + p['epsilon_1'] * H + p['epsilon_2'] * M
    dE = lam * S - (p['mu_1'] + p['beta']) * E
    dI = p['beta'] * E - (p['mu_1'] + p['mu_2'] + p['kappa']) * I
    dH = p['kappa'] * p['gamma'] * I + p['theta_2'] * M - omega_1 * H
    dM = p['kappa'] * (1 - p['gamma']) * I + p['theta_1'] * H - omega_2 * M
    dR = p['phi_1'] * H + p['phi_2'] * M - p['mu_1'] * R
    dDmod = p['pi_1'] * H + p['pi_2'] * M + p['mu_2'] * I
    return [dS, dE, dI, dH, dM, dR, dDmod]


def simulate(params_vec, t):
    alpha_hi, alpha_lo, k_sw, t_sw, kappa, I0 = params_vec
    p = dict(FIXED, alpha_hi=alpha_hi, alpha_lo=alpha_lo, k_sw=k_sw, t_sw=t_sw, kappa=kappa)
    y0 = [N_STAR - E0 - I0, E0, I0, 0, 0, 0, 0]
    sol = odeint(rhs, y0, t, args=(p,))
    return sol


def objective(params_vec):
    sol = simulate(params_vec, days.astype(float))
    E = sol[:, 1]
    modeled_cases = RHO * FIXED['beta'] * E * SCALE
    rmse_cases = np.sqrt(np.mean((modeled_cases[fit_mask] - cases_smoothed[fit_mask]) ** 2))
    return rmse_cases


bounds = [
    (0.05, 3.0),    # alpha_hi
    (0.005, 0.5),   # alpha_lo
    (0.02, 1.0),    # k_sw
    (1, 140),       # t_sw
    (0.02, 1.0),    # kappa
    (0.001, 0.5),   # I0
]


def R0(p):
    K = p['mu_1'] + p['beta']
    L = p['mu_1'] + p['mu_2'] + p['kappa']
    omega_1 = p['mu_1'] + p['epsilon_1'] + p['phi_1'] + p['pi_1'] + p['theta_1']
    omega_2 = p['mu_1'] + p['epsilon_2'] + p['phi_2'] + p['pi_2'] + p['theta_2']
    Dcal = omega_1 * omega_2 - p['theta_1'] * p['theta_2']
    P_H = p['gamma'] * omega_2 + (1 - p['gamma']) * p['theta_2']
    P_M = p['gamma'] * p['theta_1'] + (1 - p['gamma']) * omega_1
    return (p['alpha'] * p['beta'] * (Dcal + p['kappa'] * p['xi_H'] * P_H + p['kappa'] * p['eta_M'] * P_M)) / (K * L * Dcal)


if __name__ == '__main__':
    print("Running differential evolution (this may take a bit)...")
    result = differential_evolution(objective, bounds, seed=42, maxiter=250, popsize=25, tol=1e-8, polish=True, updating='deferred', workers=1)
    print("Best RMSE (daily cases, artifact days excluded):", result.fun)
    alpha_hi, alpha_lo, k_sw, t_sw, kappa, I0 = result.x
    gamma = FIXED['gamma']
    print(f"alpha_hi={alpha_hi:.5f} alpha_lo={alpha_lo:.5f} k_sw={k_sw:.5f} t_sw={t_sw:.2f} gamma={gamma:.5f} (fixed) kappa={kappa:.5f} I0={I0:.5f}")
    print(f"E0 (anchored to day-0 case count) = {E0:.5f}")

    # Now compute delta (residual CFR-adjustment factor, applied on top of the
    # population SCALE already used for cases) to match final cumulative deaths.
    # delta should be O(1) if pi_1/pi_2/mu_2 are realistic; a value far from 1
    # signals a real mismatch, not just noise.
    sol = simulate(result.x, days.astype(float))
    Dmod = sol[:, 6]  # cumulative death-flow integral, model units
    Dmod_people = Dmod * SCALE
    target_deaths = deaths_net[-1]
    delta = target_deaths / Dmod_people[-1] if Dmod_people[-1] > 0 else np.nan
    print(f"\nModel D at day 152 (real people, pre-delta): {Dmod_people[-1]:.2f}")
    print(f"Target net deaths: {target_deaths}")
    print(f"delta (residual CFR-adjustment factor): {delta:.4f}")

    p_start = dict(FIXED, alpha=alpha_hi, kappa=kappa)
    p_end = dict(FIXED, alpha=alpha_lo, kappa=kappa)
    print(f"\nR0 at wave start (alpha_hi): {R0(p_start):.4f}")
    print(f"R_eff at wave end (alpha_lo): {R0(p_end):.4f}")

    modeled_cases = RHO * FIXED['beta'] * sol[:, 1] * SCALE
    rmse_all = np.sqrt(np.mean((modeled_cases - cases_smoothed) ** 2))
    rmse_fit = np.sqrt(np.mean((modeled_cases[fit_mask] - cases_smoothed[fit_mask]) ** 2))
    print(f"\nDaily-case RMSE (excluding artifact days): {rmse_fit:.2f} cases/day")
    print(f"Daily-case RMSE (all days, incl. artifact): {rmse_all:.2f} cases/day")

    import json
    with open(RESULTS_PATH, 'w') as f:
        json.dump({
            'alpha_hi': alpha_hi, 'alpha_lo': alpha_lo, 'k_sw': k_sw, 't_sw': t_sw,
            'gamma': gamma, 'kappa': kappa, 'I0': I0, 'E0': E0,
            'delta': delta, 'N_STAR': N_STAR, 'SCALE': SCALE,
            'R0_start': R0(p_start), 'R0_end': R0(p_end),
            'rmse_cases_fit': rmse_fit, 'rmse_cases_all': rmse_all,
            'target_deaths_net': target_deaths, 'model_D_people_final': Dmod_people[-1],
        }, f, indent=2)
    print("\nSaved to uganda_fit_results.json")
