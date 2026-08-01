"""
Final Uganda Delta-wave calibration analysis: locks in the fitted parameters,
computes goodness-of-fit statistics, a perturbation-ensemble 95% CI, the
delta-sensitivity table, and generates all Section 5 figures.
"""
import json
import numpy as np
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint

HERE = Path(__file__).resolve().parent
DATA_PATH = HERE / 'data' / 'uganda_calib_data.csv'
FIGDIR = HERE / 'results' / 'figures'
FIGDIR.mkdir(parents=True, exist_ok=True)
FIGDIR = str(FIGDIR) + '/'

with open(HERE / 'results' / 'uganda_fit_results.json') as f:
    FIT = json.load(f)

days, cases_smoothed, cum_deaths = [], [], []
with open(DATA_PATH) as f:
    for r in csv.DictReader(f):
        days.append(int(r['day']))
        cases_smoothed.append(float(r['new_cases_smoothed']))
        cum_deaths.append(float(r['cum_deaths']))
days = np.array(days)
cases_smoothed = np.array(cases_smoothed)
cum_deaths = np.array(cum_deaths)
deaths_net = cum_deaths - cum_deaths[0]

ARTIFACT_DAYS = set(range(83, 90))
fit_mask = np.array([d not in ARTIFACT_DAYS for d in days])

POP = 47123531
FIXED = {
    'Lambda': 0.304, 'beta': 1/5.1,
    'theta_1': 0.396, 'theta_2': 0.0396, 'pi_1': 0.0196, 'pi_2': 0.009,
    'phi_1': 0.15, 'phi_2': 1/14, 'mu_1': 0.00576, 'mu_2': 0.0008,
    'epsilon_1': 0.001, 'epsilon_2': 0.001, 'eta_M': 0.50, 'xi_H': 0.10,
    'gamma': 0.20,
}
RHO = 0.05
N_STAR = FIXED['Lambda'] / FIXED['mu_1']
SCALE = POP / N_STAR
E0 = FIT['E0']

# Real observed Delta-wave totals (net over the 153-day window), from OWID
# archived Nov-2021 snapshot (see calibrate_uganda.py for provenance).
REAL_NET_CASES = 78410.0
REAL_NET_DEATHS = 2853.0
REAL_CFR = REAL_NET_DEATHS / REAL_NET_CASES

# True epidemic peak (artifact days 83-89 excluded)
peak_idx_real = np.argmax(np.where(fit_mask, cases_smoothed, -1))
PEAK_REAL_VALUE = cases_smoothed[peak_idx_real]
PEAK_REAL_DAY = int(days[peak_idx_real])


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


def simulate(p, t, I0):
    y0 = [N_STAR - E0 - I0, E0, I0, 0, 0, 0, 0]
    return odeint(rhs, y0, t, args=(p,))


t_fine = np.linspace(0, 152, 610)
p_fit = dict(FIXED, alpha_hi=FIT['alpha_hi'], alpha_lo=FIT['alpha_lo'],
             k_sw=FIT['k_sw'], t_sw=FIT['t_sw'], kappa=FIT['kappa'])
sol_fine = simulate(p_fit, t_fine, FIT['I0'])
modeled_cases_fine = RHO * FIXED['beta'] * sol_fine[:, 1] * SCALE
peak_idx_model = np.argmax(modeled_cases_fine)
PEAK_MODEL_VALUE = modeled_cases_fine[peak_idx_model]
PEAK_MODEL_DAY = t_fine[peak_idx_model]

sol_days = simulate(p_fit, days.astype(float), FIT['I0'])
modeled_cases_days = RHO * FIXED['beta'] * sol_days[:, 1] * SCALE
Dmod_people = sol_days[:, 6] * SCALE
delta = FIT['delta']
modeled_deaths_days = Dmod_people * delta

rmse_fit = np.sqrt(np.mean((modeled_cases_days[fit_mask] - cases_smoothed[fit_mask]) ** 2))

print(f"Real peak (artifact excluded): {PEAK_REAL_VALUE:.1f} cases/day at day {PEAK_REAL_DAY}")
print(f"Model peak: {PEAK_MODEL_VALUE:.1f} cases/day at day {PEAK_MODEL_DAY:.1f}")
print(f"Relative error at peak: {(PEAK_MODEL_VALUE-PEAK_REAL_VALUE)/PEAK_REAL_VALUE*100:.1f}%")
print(f"RMSE (fit days): {rmse_fit:.2f} cases/day")
print(f"Real CFR (Delta wave, net): {REAL_CFR*100:.3f}%")
print(f"delta (fitted): {delta:.5f}")

# ---- Perturbation ensemble for 95% CIs ----
# Perturb 5 parameters (alpha_hi, alpha_lo, t_sw, kappa, delta) by +/-15%
# (uniform), re-simulate (no re-fit), collect peak cases, day of peak,
# and cumulative deaths at day 152.
rng = np.random.default_rng(42)
N_ENS = 500
peaks, peak_days, deaths_final = [], [], []
for _ in range(N_ENS):
    a_hi = FIT['alpha_hi'] * rng.uniform(0.85, 1.15)
    a_lo = FIT['alpha_lo'] * rng.uniform(0.85, 1.15)
    tsw = FIT['t_sw'] * rng.uniform(0.85, 1.15)
    kap = FIT['kappa'] * rng.uniform(0.85, 1.15)
    delt = delta * rng.uniform(0.85, 1.15)
    p_pert = dict(FIXED, alpha_hi=a_hi, alpha_lo=a_lo, k_sw=FIT['k_sw'], t_sw=tsw, kappa=kap)
    sol_p = simulate(p_pert, t_fine, FIT['I0'])
    cases_p = RHO * FIXED['beta'] * sol_p[:, 1] * SCALE
    idx = np.argmax(cases_p)
    peaks.append(cases_p[idx])
    peak_days.append(t_fine[idx])
    deaths_final.append(sol_p[-1, 6] * SCALE * delt)

peaks = np.array(peaks); peak_days = np.array(peak_days); deaths_final = np.array(deaths_final)
peak_ci = np.percentile(peaks, [2.5, 97.5])
peak_day_ci = np.percentile(peak_days, [2.5, 97.5])
deaths_ci = np.percentile(deaths_final, [2.5, 97.5])
print(f"\nPeak cases 95% CI: [{peak_ci[0]:.0f}, {peak_ci[1]:.0f}]")
print(f"Peak day 95% CI: [{peak_day_ci[0]:.0f}, {peak_day_ci[1]:.0f}]")
print(f"Cumulative deaths 95% CI: [{deaths_ci[0]:.0f}, {deaths_ci[1]:.0f}]")

# ---- Save all results ----
results = {
    'alpha_hi': FIT['alpha_hi'], 'alpha_lo': FIT['alpha_lo'], 'k_sw': FIT['k_sw'],
    't_sw': FIT['t_sw'], 'gamma': FIXED['gamma'], 'kappa': FIT['kappa'], 'I0': FIT['I0'], 'E0': E0,
    'delta': delta, 'N_STAR': N_STAR, 'SCALE': SCALE,
    'R0_start': FIT['R0_start'], 'R0_end': FIT['R0_end'],
    'rmse_fit': rmse_fit,
    'real_net_cases': REAL_NET_CASES, 'real_net_deaths': REAL_NET_DEATHS, 'real_cfr': REAL_CFR,
    'peak_real_value': float(PEAK_REAL_VALUE), 'peak_real_day': PEAK_REAL_DAY,
    'peak_model_value': float(PEAK_MODEL_VALUE), 'peak_model_day': float(PEAK_MODEL_DAY),
    'peak_rel_error_pct': float((PEAK_MODEL_VALUE-PEAK_REAL_VALUE)/PEAK_REAL_VALUE*100),
    'peak_ci': peak_ci.tolist(), 'peak_day_ci': peak_day_ci.tolist(), 'deaths_ci': deaths_ci.tolist(),
    'modeled_deaths_final': float(modeled_deaths_days[-1]),
}
with open(HERE / 'results' / 'uganda_final_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved uganda_final_results.json")

# ---- Figures ----
# 1. Daily cases with CI band (using ensemble at each day-of-t_fine, approx via
#    same perturbation set evaluated on t_fine grid already collected as `peaks`
#    only at argmax; build a proper band by re-collecting full curves for a
#    smaller ensemble for speed)
N_BAND = 150
curves = np.zeros((N_BAND, len(t_fine)))
for i in range(N_BAND):
    a_hi = FIT['alpha_hi'] * rng.uniform(0.85, 1.15)
    a_lo = FIT['alpha_lo'] * rng.uniform(0.85, 1.15)
    tsw = FIT['t_sw'] * rng.uniform(0.85, 1.15)
    kap = FIT['kappa'] * rng.uniform(0.85, 1.15)
    p_pert = dict(FIXED, alpha_hi=a_hi, alpha_lo=a_lo, k_sw=FIT['k_sw'], t_sw=tsw, kappa=kap)
    sol_p = simulate(p_pert, t_fine, FIT['I0'])
    curves[i] = RHO * FIXED['beta'] * sol_p[:, 1] * SCALE
band_lo = np.percentile(curves, 2.5, axis=0)
band_hi = np.percentile(curves, 97.5, axis=0)

fig, ax = plt.subplots(figsize=(8, 5))
ax.fill_between(t_fine, band_lo, band_hi, color='steelblue', alpha=0.2, label='95% CI (perturbation ensemble)')
ax.plot(t_fine, modeled_cases_fine, 'r-', lw=1.8, label='Model (fitted)')
ax.plot(days[fit_mask], cases_smoothed[fit_mask], 'k.', ms=4, label='Observed (OWID, smoothed)')
ax.plot(days[~fit_mask], cases_smoothed[~fit_mask], 'x', color='gray', ms=5, label='Reporting artifact (excluded)')
ax.set_xlabel('Day (from 1 June 2021)')
ax.set_ylabel('Daily new cases')
ax.set_title('Uganda Delta wave: simulated vs. observed daily cases')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(FIGDIR + 'Dealycase95CI.png', dpi=150)
plt.close(fig)

# 2. Cumulative deaths with CI band
deaths_curves = np.zeros((N_BAND, len(t_fine)))
for i in range(N_BAND):
    a_hi = FIT['alpha_hi'] * rng.uniform(0.85, 1.15)
    a_lo = FIT['alpha_lo'] * rng.uniform(0.85, 1.15)
    tsw = FIT['t_sw'] * rng.uniform(0.85, 1.15)
    kap = FIT['kappa'] * rng.uniform(0.85, 1.15)
    delt = delta * rng.uniform(0.85, 1.15)
    p_pert = dict(FIXED, alpha_hi=a_hi, alpha_lo=a_lo, k_sw=FIT['k_sw'], t_sw=tsw, kappa=kap)
    sol_p = simulate(p_pert, t_fine, FIT['I0'])
    deaths_curves[i] = sol_p[:, 6] * SCALE * delt
dband_lo = np.percentile(deaths_curves, 2.5, axis=0)
dband_hi = np.percentile(deaths_curves, 97.5, axis=0)
modeled_deaths_fine = sol_fine[:, 6] * SCALE * delta

fig, ax = plt.subplots(figsize=(8, 5))
ax.fill_between(t_fine, dband_lo, dband_hi, color='indianred', alpha=0.2, label='95% CI')
ax.plot(t_fine, modeled_deaths_fine, 'r-', lw=1.8, label='Model (fitted)')
ax.plot(days, deaths_net, 'k.', ms=4, label='Observed (OWID)')
ax.set_xlabel('Day (from 1 June 2021)')
ax.set_ylabel('Cumulative deaths')
ax.set_title('Uganda Delta wave: simulated vs. observed cumulative deaths')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(FIGDIR + 'CumulDeathsCI.png', dpi=150)
plt.close(fig)

# 3. Compartment trajectories (per 1000: rescale model units by 1000/N_STAR)
per1000 = 1000.0 / N_STAR
S, E, I, H, M, R, D = [sol_fine[:, k] * per1000 for k in range(7)]

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(t_fine, S, 'b-')
ax.set_xlabel('Day'); ax.set_ylabel('S per 1,000'); ax.set_title('Susceptible')
plt.tight_layout(); plt.savefig(FIGDIR + 'fig_comp_susceptible.pdf'); plt.close(fig)

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(t_fine, E, 'orange', label='E')
ax.plot(t_fine, I, 'g-', label='I')
ax.set_xlabel('Day'); ax.set_ylabel('per 1,000'); ax.set_title('Exposed & Infectious'); ax.legend()
plt.tight_layout(); plt.savefig(FIGDIR + 'fig_comp_exposed_infected.pdf'); plt.close(fig)

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(t_fine, H, 'purple', label='H (hospitalized)')
ax.plot(t_fine, M, 'brown', label='M (home care)')
ax.set_xlabel('Day'); ax.set_ylabel('per 1,000'); ax.set_title('Hospitalized & Home care'); ax.legend()
plt.tight_layout(); plt.savefig(FIGDIR + 'fig_comp_hosp_homecare.pdf'); plt.close(fig)

D_people = sol_fine[:, 6] * SCALE * delta  # cumulative deaths, real people
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(t_fine, R, 'teal', label='R')
ax.set_xlabel('Day'); ax.set_ylabel('R per 1,000', color='teal')
ax2 = ax.twinx()
ax2.plot(t_fine, D_people, 'k--', label='D (cumulative deaths, real people)')
ax2.set_ylabel('Cumulative deaths (real people)')
ax.set_title('Recovered & Deaths')
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc='upper left')
plt.tight_layout(); plt.savefig(FIGDIR + 'fig_comp_recovered_deaths.pdf'); plt.close(fig)

peak_M_idx = np.argmax(M)
peak_H_at_that_time = H[peak_M_idx]
mh_ratio = M[peak_M_idx] / peak_H_at_that_time if peak_H_at_that_time > 0 else np.nan
print(f"\nPeak M (per 1000): {M[peak_M_idx]:.3f} at day {t_fine[peak_M_idx]:.1f}, H at that time: {peak_H_at_that_time:.3f}, M:H ratio = {mh_ratio:.2f}:1")
print(f"S declines from {S[0]:.1f} to {S[-1]:.1f} per 1000")

results['peak_M_per1000'] = float(M[peak_M_idx])
results['peak_H_at_Mpeak_per1000'] = float(peak_H_at_that_time)
results['MH_ratio'] = float(mh_ratio)
results['S_start_per1000'] = float(S[0])
results['S_end_per1000'] = float(S[-1])
with open(HERE / 'results' / 'uganda_final_results.json', 'w') as f:
    json.dump(results, f, indent=2)

# 4. Gamma sensitivity panels
gammas = np.linspace(0.05, 0.5, 20)
peakM_vs_gamma, deaths_vs_gamma = [], []
for g in gammas:
    p_g = dict(p_fit, gamma=g)
    sol_g = simulate(p_g, t_fine, FIT['I0'])
    Mg = sol_g[:, 4] * per1000
    peakM_vs_gamma.append(Mg.max())
    deaths_vs_gamma.append(sol_g[-1, 6] * SCALE * delta)

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(gammas, peakM_vs_gamma, 'o-', color='brown')
ax.axvline(FIXED['gamma'], color='gray', ls='--', lw=1)
ax.set_xlabel(r'$\gamma$ (hospitalization proportion)')
ax.set_ylabel('Peak M per 1,000')
ax.set_title('Peak home-care prevalence vs. $\\gamma$')
plt.tight_layout(); plt.savefig(FIGDIR + 'fig_gamma_homecare.pdf'); plt.close(fig)

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(gammas, deaths_vs_gamma, 'o-', color='indianred')
ax.axvline(FIXED['gamma'], color='gray', ls='--', lw=1)
ax.set_xlabel(r'$\gamma$ (hospitalization proportion)')
ax.set_ylabel('Cumulative deaths')
ax.set_title('Cumulative mortality vs. $\\gamma$')
plt.tight_layout(); plt.savefig(FIGDIR + 'fig_gamma_mortality.pdf'); plt.close(fig)

# report specific gamma=0.05 and gamma=0.25 values (bracketing fitted 0.20)
def eval_gamma(g):
    p_g = dict(p_fit, gamma=g)
    sol_g = simulate(p_g, t_fine, FIT['I0'])
    Mg = sol_g[:, 4] * per1000
    return Mg.max(), sol_g[-1, 6] * SCALE * delta

m05, d05 = eval_gamma(0.05)
m25, d25 = eval_gamma(0.25)
m_base, d_base = eval_gamma(0.20)
print(f"\ngamma=0.05: peak M={m05:.2f}/1000, deaths={d05:.0f} ({(d05-d_base)/d_base*100:+.1f}%)")
print(f"gamma=0.20 (fitted): peak M={m_base:.2f}/1000, deaths={d_base:.0f}")
print(f"gamma=0.25: peak M={m25:.2f}/1000, deaths={d25:.0f} ({(d25-d_base)/d_base*100:+.1f}%)")

results.update({
    'gamma_sens_005_peakM': float(m05), 'gamma_sens_005_deaths': float(d05),
    'gamma_sens_025_peakM': float(m25), 'gamma_sens_025_deaths': float(d25),
    'gamma_sens_base_peakM': float(m_base), 'gamma_sens_base_deaths': float(d_base),
})
with open(HERE / 'results' / 'uganda_final_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nAll figures saved to", FIGDIR)
