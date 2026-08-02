"""
Generates the manuscript's actual Figure 2 (fig_corrected_coeff_a_v2.png):
center-manifold coefficient a, swept over the hospital recovery rate phi_1
and the return-to-susceptibility rates epsilon_1=epsilon_2, at the Table 2
baseline. Uses the same center-manifold formula as model.py / Section 3.

This is the real manuscript figure, distinct from
bifurcation_verification.py's own independent verification_* plot of the
same underlying quantity (see that script's docstring).
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import BASELINE, derived

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results' / 'figures'
OUT.mkdir(parents=True, exist_ok=True)


def center_manifold_a(p):
    K, L, omega_1, omega_2, Dcal, P_H, P_M = derived(p)
    alpha_c = (p['mu_1'] + p['beta']) * L * Dcal / (
        p['beta'] * (Dcal + p['kappa'] * p['xi_H'] * P_H + p['kappa'] * p['eta_M'] * P_M))
    w_E = L / p['beta']
    w_H = p['kappa'] * P_H / Dcal
    w_M = p['kappa'] * P_M / Dcal
    w_R = (p['phi_1'] * w_H + p['phi_2'] * w_M) / p['mu_1']
    v_E = p['beta'] / (p['mu_1'] + p['beta'])
    B = 1 + p['xi_H'] * w_H + p['eta_M'] * w_M
    A_prime = w_E + 1 + w_H + w_M + w_R
    S0 = p['Lambda'] / p['mu_1']
    a = -(2 * alpha_c * v_E / S0) * B * A_prime
    return a


phi1_range = np.linspace(0.02, 0.95, 60)
eps_values = [0.0, 0.001, 0.02, 0.10, 0.30]

fig, ax = plt.subplots(figsize=(7, 5))
a_max_overall = -np.inf
for eps in eps_values:
    a_vals = []
    for phi1 in phi1_range:
        p = dict(BASELINE, phi_1=phi1, epsilon_1=eps, epsilon_2=eps)
        a = center_manifold_a(p)
        a_vals.append(a)
        a_max_overall = max(a_max_overall, a)
    ax.plot(phi1_range, a_vals, lw=2.4, label=f"eps1=eps2={eps}")

ax.axhline(0, color='black', linewidth=1.5)
ax.axvline(0.15, color='gray', ls=':', lw=1.5, label='baseline phi1=0.15')
ax.set_xlabel('Hospital recovery rate phi_1')
ax.set_ylabel('center-manifold coefficient a')
ax.set_title('a is negative throughout (forward bifurcation only)')
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(OUT / 'fig_corrected_coeff_a_v2.png', dpi=300)
plt.close(fig)

print(f"max(a) over full sweep = {a_max_overall:.6e} (should be < 0)")
assert a_max_overall < 0, "a was not strictly negative everywhere -- investigate before trusting the figure"
print(f"Saved to {OUT / 'fig_corrected_coeff_a_v2.png'}")
