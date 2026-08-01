"""
SEIHM-R-D core model, exactly matching manuscript.tex Eq. (eq:system), (eq:R0), (eq:ab_corrected).
Reconciled, literature-anchored baseline parameters (see param_sources.md for full citations).
"""
import numpy as np
from scipy.integrate import odeint

BASELINE = {
    'Lambda': 0.304,      # verified real: Ndondo et al. 2021, Results in Physics 24:104096
    'alpha': 0.035,        # kept original value; no portable literature value exists for a
                            #   model-specific contact-rate coefficient; relabeled "assumed" (as kappa already was)
    'beta': 1/5.1,          # updated: real incubation period, McAloon et al. 2020 BMJ Open;
                            #   Lauer et al. 2020 Ann Intern Med (corroborating)
    'gamma': 0.20,          # updated: WHO-China Joint Mission 2020 severity split (13.8%+6.1%)
    'kappa': 0.100,         # unchanged, calibrated (as in the original table)
    'theta_1': 0.396,       # kept original value; relabeled "assumed" (no clean literature analog)
    'theta_2': 0.0396,      # kept original value; relabeled "assumed"
    'pi_1': 0.0196,         # kept original value; relabeled "assumed, order-of-magnitude
                            #   consistent with hospitalized-case CFR literature"
    'pi_2': 0.009,          # kept original value; relabeled "assumed, order-of-magnitude
                            #   consistent with community-case CFR literature"
    'phi_1': 0.15,          # updated (was 0.900, implying a physically impossible ~1.1-day
                            #   stay): ~6.7-day rate, near the Africa-specific subgroup of
                            #   Alimohamadi et al. 2022's LOS meta-analysis (8.56d, wide CI
                            #   reflecting capacity-driven variability)
    'phi_2': 1/14,           # updated (was 0.0000625, physically impossible): WHO recovery
                            #   guidance; Tenforde et al. 2020 MMWR (corroborating)
    'mu_1': 0.00576,        # kept value; citation corrected to Sinan et al. 2021 (exact match there)
    'mu_2': 0.0008,          # updated (was 0.196, physically impossible -- exceeded its own
                            #   cited source's stated upper bound): Verity et al. 2020, Lancet
                            #   Infect Dis (IFR + time-to-death)
    'epsilon_1': 0.001,     # kept; no literature analog exists; assumed
    'epsilon_2': 0.001,     # kept; matches epsilon_1; assumed
    'eta_M': 0.50,           # kept; relabeled "assumed, informed by household-transmission literature"
    'xi_H': 0.10,            # kept; relabeled "assumed, informed by nosocomial-transmission literature"
}


def derived(p):
    K = p['mu_1'] + p['beta']
    L = p['mu_1'] + p['mu_2'] + p['kappa']
    omega_1 = p['mu_1'] + p['epsilon_1'] + p['phi_1'] + p['pi_1'] + p['theta_1']
    omega_2 = p['mu_1'] + p['epsilon_2'] + p['phi_2'] + p['pi_2'] + p['theta_2']
    Dcal = omega_1 * omega_2 - p['theta_1'] * p['theta_2']
    P_H = p['gamma'] * omega_2 + (1 - p['gamma']) * p['theta_2']
    P_M = p['gamma'] * p['theta_1'] + (1 - p['gamma']) * omega_1
    return K, L, omega_1, omega_2, Dcal, P_H, P_M


def R0(p):
    K, L, omega_1, omega_2, Dcal, P_H, P_M = derived(p)
    return (p['alpha'] * p['beta'] * (Dcal + p['kappa'] * p['xi_H'] * P_H
            + p['kappa'] * p['eta_M'] * P_M)) / (K * L * Dcal)


def alpha_for_R0(p, target_R0):
    """Solve for alpha giving a target R0 (R0 is exactly linear in alpha)."""
    p1 = dict(p, alpha=1.0)
    return target_R0 / R0(p1)


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
    return alpha_c, a


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
    dD = p['pi_1'] * H + p['pi_2'] * M + p['mu_2'] * I
    return [dS, dE, dI, dH, dM, dR, dD]


def simulate(p, t, y0):
    return odeint(rhs, y0, t, args=(p,))


def default_ic(N=1000.0):
    return [N - 5, 3, 2, 0, 0, 0, 0]
