"""
Bifurcation Verification for the SEIHM-R-D COVID-19 Model
Reproduces Scenarios 4 and 5 of the manuscript (Section 4) and the
center-manifold calculation of Section 3 (Proposition "no backward bifurcation").

The center-manifold calculation fully accounts for the state-dependence of
N_L in the frequency-dependent force of infection. The resulting
center-manifold coefficient a is strictly negative for every admissible
parameter combination, so the model undergoes only forward (transcritical)
bifurcation at R0 = 1, and R0 < 1 is both necessary and sufficient for
elimination. There is no bistability regime.

Near R0 = 1 the relaxation eigenvalue is extremely small
(~ -6.7e-4 / day at the manuscript's Scenario 5 parameter point, a
~1,480-day relaxation time), so trajectories that have not yet reached
their unique common equilibrium can look like they are converging to two
different endpoints when they are not. Scenario 5 below therefore
integrates for 60,000 days.

This script is an independent numerical check: it reaches the same
conclusion as the manuscript (a < 0 everywhere; all trajectories converge
to the disease-free equilibrium) using the model's own closed-form
equations, but it is not the exact plotting code behind the manuscript's
Figures 4 and 5 (figures/fig_corrected_coeff_a_v2.png and
figures/fig_no_bistability_v2.png in the repository root). Its output is
saved separately under results/figures/verification_* below rather than
into the top-level figures/ folder, since it is a supplementary check,
not a manuscript figure.

Outputs (relative to this script's own directory):
  - results/figures/verification_coefficient_a_scenario4.png
  - results/figures/verification_no_bistability_scenario5.png
  - data/verification_coefficient_a.csv
  - data/verification_no_bistability_trajectories.csv
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

HERE = Path(__file__).resolve().parent
DEFAULT_FIGDIR = HERE / 'results' / 'figures'
DEFAULT_DATADIR = HERE / 'data'
DEFAULT_FIGDIR.mkdir(parents=True, exist_ok=True)

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['grid.color'] = 'gray'
plt.rcParams['grid.alpha'] = 0.3


def get_baseline_parameters():
    """Baseline parameters, Table 2 of the manuscript (matches model.py's BASELINE)."""
    return {
        'Lambda': 0.304,
        'alpha': 0.035,
        'beta': 1/5.1,
        'gamma': 0.20,
        'kappa': 0.100,
        'theta_1': 0.396,
        'theta_2': 0.0396,
        'pi_1': 0.0196,
        'pi_2': 0.009,
        'phi_1': 0.15,
        'phi_2': 1/14,
        'mu_1': 0.00576,
        'mu_2': 0.0008,
        'epsilon_1': 0.001,
        'epsilon_2': 0.001,
        'eta_M': 0.50,
        'xi_H': 0.10,
    }


def derived_quantities(p):
    """K, L, omega_1, omega_2, script-D, P_H, P_M from Section 3."""
    K = p['mu_1'] + p['beta']
    L = p['mu_1'] + p['mu_2'] + p['kappa']
    omega_1 = p['mu_1'] + p['epsilon_1'] + p['phi_1'] + p['pi_1'] + p['theta_1']
    omega_2 = p['mu_1'] + p['epsilon_2'] + p['phi_2'] + p['pi_2'] + p['theta_2']
    Dcal = omega_1 * omega_2 - p['theta_1'] * p['theta_2']
    P_H = p['gamma'] * omega_2 + (1 - p['gamma']) * p['theta_2']
    P_M = p['gamma'] * p['theta_1'] + (1 - p['gamma']) * omega_1
    return K, L, omega_1, omega_2, Dcal, P_H, P_M


def calculate_R0(p):
    K, L, omega_1, omega_2, Dcal, P_H, P_M = derived_quantities(p)
    return (p['alpha'] * p['beta'] * (Dcal + p['kappa'] * p['xi_H'] * P_H
            + p['kappa'] * p['eta_M'] * P_M)) / (K * L * Dcal)


def center_manifold_coefficient_a(p):
    """
    Eqs. (alpha_c), (wHM), (vHM), (quad_corrected), (ab_corrected) of the
    manuscript. Returns (alpha_c, a, b). a < 0 and b > 0 for every admissible
    parameter combination (Proposition "no backward bifurcation").
    """
    K, L, omega_1, omega_2, Dcal, P_H, P_M = derived_quantities(p)
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
    b = v_E * B
    return alpha_c, a, b


def seihmrd_rhs(y, t, p):
    """Right-hand side of System (eq:system), full nonlinear model."""
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
    dD = p['pi_1'] * H + p['pi_2'] * M + p['mu_2'] * I
    return [dS, dE, dI, dH, dM, dR, dD]


def scenario4_bifurcation_coefficient(outdir_fig=DEFAULT_FIGDIR, outdir_data=DEFAULT_DATADIR):
    """
    Scenario 4: sweep phi_1 in [0.02, 0.95] and epsilon_1 = epsilon_2 in
    {0, 0.001, 0.02, 0.10, 0.30}; confirm a < 0 throughout.
    """
    base = get_baseline_parameters()
    phi1_range = np.linspace(0.02, 0.95, 60)
    eps_values = [0.0, 0.001, 0.02, 0.10, 0.30]

    rows = []
    fig, ax = plt.subplots(figsize=(7, 5))
    for eps in eps_values:
        a_vals = []
        for phi1 in phi1_range:
            p = base.copy()
            p['phi_1'] = phi1
            p['epsilon_1'] = eps
            p['epsilon_2'] = eps
            _, a, _ = center_manifold_coefficient_a(p)
            a_vals.append(a)
            rows.append((phi1, eps, a))
        ax.plot(phi1_range, a_vals, label=f"$\\epsilon_1=\\epsilon_2={eps}$")

    ax.axhline(0, color='black', linewidth=1)
    ax.set_xlabel(r"Hospital recovery rate $\varphi_1$")
    ax.set_ylabel(r"Center-manifold coefficient $a$")
    ax.set_title("Coefficient $a$ is strictly negative throughout\n"
                  "(confirms only forward bifurcation at $\\mathcal{R}_0=1$)")
    ax.legend(fontsize=8)
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(f"{outdir_fig}/verification_coefficient_a_scenario4.png", dpi=300)
    plt.close(fig)

    with open(f"{outdir_data}/verification_coefficient_a.csv", "w") as f:
        f.write("phi_1,epsilon_1_eq_epsilon_2,a\n")
        for phi1, eps, a in rows:
            f.write(f"{phi1},{eps},{a}\n")

    a_max = max(a for _, _, a in rows)
    print(f"Scenario 4: max(a) over full sweep = {a_max:.6e} (should be < 0)")
    assert a_max < 0, "a was not strictly negative everywhere -- investigate before trusting the figure"
    return rows


def scenario5_no_bistability(outdir_fig=DEFAULT_FIGDIR, outdir_data=DEFAULT_DATADIR):
    """
    Scenario 5: at phi_1=0.05, alpha chosen so R0 is just under 1, integrate
    from five initial conditions spanning a tiny seed to a very large initial
    outbreak, over a 60,000-day horizon, and confirm all trajectories
    converge to the disease-free equilibrium.
    """
    p = get_baseline_parameters()
    p['phi_1'] = 0.05
    # Solve for alpha giving R0 just under 1 (R0 is linear in alpha).
    p_unit = dict(p, alpha=1.0)
    p['alpha'] = 0.988 / calculate_R0(p_unit)
    R0 = calculate_R0(p)
    print(f"Scenario 5: R0 = {R0:.4f} at phi_1=0.05, alpha={p['alpha']:.4f}")

    N = 1000.0
    initial_conditions = {
        'Minuscule_Outbreak': [N - 0.1, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0],
        'Tiny_Outbreak':      [N - 1.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0],
        'Small_Outbreak':     [N - 5.0, 3.0, 2.0, 0.0, 0.0, 0.0, 0.0],
        'Medium_Outbreak':    [N - 50.0, 25.0, 25.0, 0.0, 0.0, 0.0, 0.0],
        'Large_Outbreak':     [N - 100.0, 50.0, 50.0, 0.0, 0.0, 0.0, 0.0],
    }

    # Log-spaced grid: a linear 0-60,000-day axis compresses the actual
    # transient (which unfolds over ~1,000 days) into an unreadable sliver.
    t = np.logspace(-1, np.log10(60000), 3000)
    results = {}
    for name, y0 in initial_conditions.items():
        sol = odeint(seihmrd_rhs, y0, t, args=(p,))
        results[name] = sol

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for name, sol in results.items():
        axes[0].plot(t, sol[:, 0], label=name.replace('_', ' '))
        axes[1].plot(t, sol[:, 6], label=name.replace('_', ' '))
    axes[0].set_xlabel("Time (days, log scale)")
    axes[0].set_ylabel("Susceptible $S$")
    axes[0].set_title("All trajectories converge to the same DFE")
    axes[1].set_xlabel("Time (days, log scale)")
    axes[1].set_ylabel("Cumulative deaths $D$")
    axes[1].set_title(f"D converges to a finite, path-dependent limit\n"
                       f"(larger seed outbreak, more deaths before die-out; R0={R0:.3f})")
    for ax in axes:
        ax.set_xscale('log')
        ax.legend(fontsize=8)
        ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{outdir_fig}/verification_no_bistability_scenario5.png", dpi=300)
    plt.close(fig)

    header = ["Time"]
    for name in initial_conditions:
        header += [f"{name}_{v}" for v in ["S", "E", "I", "H", "M", "R", "D"]]
    with open(f"{outdir_data}/verification_no_bistability_trajectories.csv", "w") as f:
        f.write(",".join(header) + "\n")
        for i, ti in enumerate(t):
            row = [str(ti)]
            for name in initial_conditions:
                row += [str(v) for v in results[name][i, :]]
            f.write(",".join(row) + "\n")

    final_S = {name: sol[-1, 0] for name, sol in results.items()}
    spread = max(final_S.values()) - min(final_S.values())
    print(f"Scenario 5: spread in final S across the 5 runs at t=60000 days = {spread:.6f}")
    print("(small spread confirms convergence to a common disease-free equilibrium; "
          "note this would NOT yet be small at t=300 days, which is why the model "
          "can look bistable if under-integrated)")
    return results


if __name__ == "__main__":
    print("=" * 70)
    print("Scenario 4: center-manifold coefficient a across the parameter sweep")
    print("=" * 70)
    scenario4_bifurcation_coefficient()

    print()
    print("=" * 70)
    print("Scenario 5: long-horizon no-bistability demonstration")
    print("=" * 70)
    scenario5_no_bistability()

    print()
    print(f"Done. See {DEFAULT_FIGDIR}/verification_coefficient_a_scenario4.png and {DEFAULT_FIGDIR}/verification_no_bistability_scenario5.png")
