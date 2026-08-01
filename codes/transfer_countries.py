"""
Cross-country parameter transfer: refit only (alpha, beta, kappa) per country
to match real cumulative-mortality curves (OWID), holding the other 14 of 17
Table 2 parameters at Uganda-calibrated/literature values (including Uganda's
fitted delta and gamma=0.20).
"""
import json
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import differential_evolution
from sklearn.metrics import r2_score

HERE = Path(__file__).resolve().parent
FIGDIR_PATH = HERE / 'results' / 'figures'
FIGDIR_PATH.mkdir(parents=True, exist_ok=True)
FIGDIR = str(FIGDIR_PATH) + '/'

with open(HERE / 'results' / 'uganda_final_results.json') as f:
    UGA = json.load(f)
with open(HERE / 'data' / 'transfer_data.json') as f:
    TDATA = json.load(f)

# World Bank 2020 total-population estimates (standard, widely cited vintage).
POPULATIONS = {'MOZ': 31178238, 'SEN': 16743930, 'CMR': 26545864}

FIXED_BASE = {
    'Lambda': 0.304,
    'theta_1': 0.396, 'theta_2': 0.0396, 'pi_1': 0.0196, 'pi_2': 0.009,
    'phi_1': 0.15, 'phi_2': 1/14, 'mu_1': 0.00576, 'mu_2': 0.0008,
    'epsilon_1': 0.001, 'epsilon_2': 0.001, 'eta_M': 0.50, 'xi_H': 0.10,
    'gamma': 0.20,
}
RHO = 0.05
N_STAR = FIXED_BASE['Lambda'] / FIXED_BASE['mu_1']
DELTA = UGA['delta']  # held fixed at Uganda-fitted value


def rhs(y, t, p):
    S, E, I, H, M, R, D = y
    N_L = S + E + I + H + M + R
    lam = p['alpha'] * (I + p['eta_M'] * M + p['xi_H'] * H) / N_L
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


results = {}
for code, meta in TDATA.items():
    name = meta['name']
    cum_deaths = np.array(meta['cum_deaths'])
    deaths_net = cum_deaths - cum_deaths[0]
    n = len(cum_deaths)
    t = np.arange(n).astype(float)
    pop = POPULATIONS[code]
    scale = pop / N_STAR
    new_cases_day0 = float(meta['new_cases_day0'])

    def simulate(params_vec):
        alpha, beta, kappa, I0 = params_vec
        E0 = new_cases_day0 / (RHO * beta * scale)
        p = dict(FIXED_BASE, alpha=alpha, beta=beta, kappa=kappa)
        y0 = [N_STAR - E0 - I0, E0, I0, 0, 0, 0, 0]
        return odeint(rhs, y0, t, args=(p,)), E0

    def objective(params_vec):
        sol, _ = simulate(params_vec)
        Dmod_people = sol[:, 6] * scale * DELTA
        return np.sqrt(np.mean((Dmod_people - deaths_net) ** 2))

    bounds = [(0.02, 1.5), (1/20, 1/2), (0.02, 1.0), (0.0001, 0.5)]
    res = differential_evolution(objective, bounds, seed=7, maxiter=300, popsize=25,
                                  tol=1e-10, polish=True, updating='deferred', workers=1)
    alpha_f, beta_f, kappa_f, I0_f = res.x
    sol, E0_f = simulate(res.x)
    Dmod_people = sol[:, 6] * scale * DELTA
    r2 = r2_score(deaths_net, Dmod_people)
    mape = np.nanmean(np.abs((Dmod_people - deaths_net) / np.where(deaths_net == 0, np.nan, deaths_net))) * 100

    print(f"{name}: alpha={alpha_f:.4f} beta={beta_f:.4f} (1/beta={1/beta_f:.2f}d) kappa={kappa_f:.4f} I0={I0_f:.4f}")
    print(f"  R2={r2:.4f}  MAPE={mape:.2f}%  RMSE={res.fun:.2f} deaths")

    results[code] = {
        'name': name, 'alpha': alpha_f, 'beta': beta_f, 'kappa': kappa_f, 'I0': I0_f, 'E0': E0_f,
        'R2': r2, 'MAPE': mape, 'RMSE_deaths': res.fun, 'population': pop,
        'final_deaths_obs': float(deaths_net[-1]), 'final_deaths_model': float(Dmod_people[-1]),
        'n_days': n,
    }

    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    ax.plot(t, deaths_net, 'k.', ms=3, label='Observed')
    ax.plot(t, Dmod_people, 'r-', lw=1.6, label='Model (3-param transfer)')
    ax.set_xlabel('Day'); ax.set_ylabel('Cumulative deaths')
    ax.set_title(f'{name} ($R^2$={r2:.3f})')
    ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(FIGDIR + f'{name}_Validation_Clean.pdf')
    plt.close(fig)

mean_r2 = np.mean([results[c]['R2'] for c in results])
print(f"\nMean R2 across 3 countries: {mean_r2:.4f}")
results['mean_R2'] = float(mean_r2)

with open(HERE / 'results' / 'transfer_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved transfer_results.json")
