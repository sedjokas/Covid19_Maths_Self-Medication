"""
export_gamma_panels.py
======================
Exports the two panels of fig_val_gamma_sensitivity.pdf as separate PDFs,
without any title, suptitle, or caption text.

  fig_gamma_homecare.pdf  – panel (a): Home-Care Burden M(t)
  fig_gamma_mortality.pdf – panel (b): Cumulative Mortality
"""

import os
import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

OUT = os.path.dirname(os.path.abspath(__file__))

# ── Parameters (identical to simulate_uganda.py) ─────────────
N_pop    = 47_123_531
alpha_hi = 0.351
alpha_lo = 0.0752
t_sw     = 52
k_sw     = 5.0

beta   = 0.095
gamma  = 0.15
kappa  = 0.100
theta1 = 0.396
theta2 = 0.0396
pi1    = 0.0196
pi2    = 0.009
phi1   = 0.900
phi2   = 6.25e-5
mu1    = 0.00576
mu2    = 0.196
eps1   = 0.001
eps2   = 0.001
eta_M  = 0.50
xi_H   = 0.10
Lambda = 5219.0 / N_pop
rho    = 0.05
scale  = N_pop / 1000.0

omega1 = mu1 + eps1 + phi1 + pi1 + theta1
omega2 = mu1 + eps2 + phi2 + pi2 + theta2
K      = mu1 + beta
L      = mu1 + mu2 + kappa
D_val  = omega1 * omega2 - theta1 * theta2

# ── Dynamic QSS ICs ──────────────────────────────────────────
r_dyn = 0.042
I_E   = beta / (L + r_dyn)
D_r   = (omega1 + r_dyn) * (omega2 + r_dyn) - theta1 * theta2
PH_r  = gamma  * (omega2 + r_dyn) + (1 - gamma)  * theta2
PM_r  = (1 - gamma) * (omega1 + r_dyn) + gamma * theta1
H_E   = kappa * I_E * PH_r / D_r
M_E   = kappa * I_E * PM_r / D_r

E0  = 200.0 / (rho * beta * scale)
I0  = I_E * E0
H0  = H_E * E0
M0  = M_E * E0
S0  = 968.0
R0v = 1000.0 - S0 - E0 - I0 - H0 - M0
D0  = 0.0

# ── Mortality scale (from baseline run) ──────────────────────
def alpha_t(t):
    return alpha_lo + (alpha_hi - alpha_lo) / (1.0 + np.exp(k_sw * (t - t_sw)))

def seimhrd(y, t):
    S, E, I, H, M, R, Dc = y
    a   = alpha_t(t)
    NL  = max(S + E + I + H + M + R, 1e-9)
    lam = a * S * (I + eta_M * M + xi_H * H) / NL
    dS  = Lambda * NL - mu1 * S - lam + eps1 * H + eps2 * M
    dE  = lam - K * E
    dI  = beta * E - L * I
    dH  = kappa * gamma       * I + theta2 * M - omega1 * H
    dM  = kappa * (1 - gamma) * I + theta1 * H - omega2 * M
    dR  = phi1 * H + phi2 * M - mu1 * R
    dDc = pi1  * H + pi2  * M + mu2 * I
    return [dS, dE, dI, dH, dM, dR, dDc]

t = np.linspace(0, 152, 1521)
y0 = [S0, E0, I0, H0, M0, R0v, D0]
sol = odeint(seimhrd, y0, t, rtol=1e-9, atol=1e-11)
D_t = sol[:, 6]
cum_d_raw  = D_t * scale
mort_scale = 1448.0 / max(cum_d_raw[-1], 1.0)

# ── Scenarios ────────────────────────────────────────────────
scenarios = [
    (0.05, "#d62728", r"$\gamma=0.05$  (very low access, 95% home-care)"),
    (0.15, "#1f77b4", r"$\gamma=0.15$  (calibrated, Uganda baseline)"),
    (0.25, "#2ca02c", r"$\gamma=0.25$  (improved access)"),
]

M_results  = []
D_results  = []

for gv, col, lbl in scenarios:
    PH_r_g = gv * (omega2 + r_dyn) + (1 - gv) * theta2
    PM_r_g = (1 - gv) * (omega1 + r_dyn) + gv * theta1
    H_E_g  = kappa * I_E * PH_r_g / D_r
    M_E_g  = kappa * I_E * PM_r_g / D_r
    y0_g   = [S0, E0, I_E*E0, H_E_g*E0, M_E_g*E0,
              1000 - S0 - E0 - I_E*E0 - H_E_g*E0 - M_E_g*E0, D0]

    def sys_g(y, t_, gv_=gv):
        S_, E_, I_, H_, M_, R_, Dc_ = y
        a_  = alpha_t(t_)
        NL_ = max(S_ + E_ + I_ + H_ + M_ + R_, 1e-9)
        lm_ = a_ * S_ * (I_ + eta_M*M_ + xi_H*H_) / NL_
        return [
            Lambda*NL_ - mu1*S_ - lm_ + eps1*H_ + eps2*M_,
            lm_          - K * E_,
            beta * E_    - L * I_,
            kappa * gv_        * I_ + theta2*M_ - omega1*H_,
            kappa * (1 - gv_)  * I_ + theta1*H_ - omega2*M_,
            phi1*H_ + phi2*M_  - mu1*R_,
            pi1*H_  + pi2*M_   + mu2*I_,
        ]

    sg = odeint(sys_g, y0_g, t, rtol=1e-9, atol=1e-11)
    M_results.append((col, lbl, sg[:, 4]))
    D_results.append((col, lbl, sg[:, 6] * scale * mort_scale))

# ── Figure style ─────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "serif",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.grid":         True,
    "grid.alpha":        0.3,
    "figure.dpi":        150,
})

# ── Panel (a): Home-Care Burden ───────────────────────────────
fig, ax = plt.subplots(figsize=(6, 4.5))
for col, lbl, M_g in M_results:
    ax.plot(t, M_g, color=col, lw=2.5, label=lbl)
ax.set_xlabel("Days since 1 June 2021", fontsize=11)
ax.set_ylabel("Home-care $M(t)$ per 1,000", fontsize=11)
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig_gamma_homecare.pdf"), bbox_inches="tight")
plt.close()
print("Saved fig_gamma_homecare.pdf")

# ── Panel (b): Cumulative Mortality ───────────────────────────
fig, ax = plt.subplots(figsize=(6, 4.5))
for col, lbl, D_g_cal in D_results:
    ax.plot(t, D_g_cal, color=col, lw=2.5, label=lbl)
    ax.annotate(f"{int(D_g_cal[-1]):,}", xy=(152, D_g_cal[-1]),
                fontsize=8, color=col, xytext=(145, D_g_cal[-1] + 20))
ax.set_xlabel("Days since 1 June 2021", fontsize=11)
ax.set_ylabel("Cumulative deaths (calibrated)", fontsize=11)
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig_gamma_mortality.pdf"), bbox_inches="tight")
plt.close()
print("Saved fig_gamma_mortality.pdf")

# ════════════════════════════════════════════════════════════
# Figure 3 – 4 individual compartment panels (no title)
# ════════════════════════════════════════════════════════════

# Run baseline simulation to get all compartments
sol_base = odeint(seimhrd, y0, t, rtol=1e-9, atol=1e-11)
S_t, E_t, I_t, H_t, M_t, R_t, D_t2 = sol_base.T
cum_d_cal = D_t2 * scale * mort_scale
daily_rep  = rho * beta * E_t * scale
peak_day   = t[np.argmax(daily_rep)]

# Panel (a): Susceptible S(t)
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(t, S_t, color="#1f77b4", lw=2)
ax.set_xlabel("Days since 1 June 2021", fontsize=11)
ax.set_ylabel("Per 1,000 individuals", fontsize=11)
ax.set_ylim(bottom=0)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig_comp_susceptible.pdf"), bbox_inches="tight")
plt.close()
print("Saved fig_comp_susceptible.pdf")

# Panel (b): Exposed E(t) and Infected I(t)
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(t, E_t, color="#ff7f0e", lw=2, ls="--", label="Exposed $E$")
ax.plot(t, I_t, color="#d62728", lw=2, label="Infected $I$")
ax.axvline(peak_day, color="gray", ls=":", lw=1)
ax.set_xlabel("Days since 1 June 2021", fontsize=11)
ax.set_ylabel("Per 1,000 individuals", fontsize=11)
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig_comp_exposed_infected.pdf"), bbox_inches="tight")
plt.close()
print("Saved fig_comp_exposed_infected.pdf")

# Panel (c): Hospitalised H(t) vs Home-care M(t)
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(t, H_t, color="#17becf", lw=2, label="Hospitalised $H$")
ax.plot(t, M_t, color="#9467bd", lw=2, label="Home-care $M$")
peak_M = np.max(M_t); peak_H = np.max(H_t)
ax.annotate(
    f"Peak $M/H \\approx {peak_M/peak_H:.1f}$",
    xy=(t[np.argmax(M_t)], peak_M),
    xytext=(80, peak_M * 0.80),
    arrowprops=dict(arrowstyle="->"), fontsize=9)
ax.set_xlabel("Days since 1 June 2021", fontsize=11)
ax.set_ylabel("Per 1,000 individuals", fontsize=11)
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig_comp_hosp_homecare.pdf"), bbox_inches="tight")
plt.close()
print("Saved fig_comp_hosp_homecare.pdf")

# Panel (d): Recovered R(t) and Cumulative Deaths D(t)
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(t, R_t, color="#2ca02c", lw=2, label="Recovered $R$")
ax2 = ax.twinx()
ax2.plot(t, cum_d_cal, color="#d62728", lw=2, ls="--",
         label="Cumul. deaths $D$ (calibrated)")
ax2.set_ylabel("Cumulative deaths (absolute)", fontsize=9)
ax2.spines["right"].set_visible(True)
ax.set_xlabel("Days since 1 June 2021", fontsize=11)
ax.set_ylabel("Per 1,000 individuals", fontsize=11)
l1, lb1 = ax.get_legend_handles_labels()
l2, lb2 = ax2.get_legend_handles_labels()
ax.legend(l1 + l2, lb1 + lb2, fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig_comp_recovered_deaths.pdf"), bbox_inches="tight")
plt.close()
print("Saved fig_comp_recovered_deaths.pdf")

# ════════════════════════════════════════════════════════════
# fig_val_daily_cases.pdf – without title
# ════════════════════════════════════════════════════════════
obs_days = np.array([  0,  7, 14, 21, 28, 35, 42, 49,
                       52, 56, 63, 70, 77, 84, 91, 98,
                      105, 112, 119, 126, 133, 140, 147, 152])
obs_cases = np.array([ 200,  290,  430,  680, 1010, 1310, 1590, 1760,
                      1800, 1740, 1510, 1230,  960,  740,  570,  430,
                       330,  250,  200,  160,  130,  110,   90,   80])

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(t, daily_rep, color="#1f77b4", lw=2.5,
        label=r"Model  ($\rho=5\%$, $\mathcal{R}_0=2.23$)")
ax.plot(obs_days, obs_cases, "ro", ms=5, zorder=5,
        label="Observed 7-day avg. (OWID)")
ax.axvline(peak_day, color="#1f77b4", ls="--", lw=1.2,
           label=f"Model peak – day {peak_day:.0f}")
ax.axvline(52, color="salmon", ls=":", lw=1.2,
           label="Observed peak – day 51–57")
ax.set_xlabel("Days since 1 June 2021", fontsize=11)
ax.set_ylabel("Daily confirmed cases", fontsize=11)
ax.legend(fontsize=9)
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig_daily_cases_notitle.pdf"), bbox_inches="tight")
plt.close()
print("Saved fig_daily_cases_notitle.pdf")

# ════════════════════════════════════════════════════════════
# fig_val_deaths.pdf – without title
# ════════════════════════════════════════════════════════════
obs_cum_deaths = np.array([
       0,   12,   32,   70,  130,  210,  310,  430,
     480,  540,  640,  750,  840,  930, 1010, 1080,
    1140, 1200, 1260, 1310, 1360, 1400, 1430, 1448])

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.fill_between(t, cum_d_cal, alpha=0.15, color="#1f77b4")
ax.plot(t, cum_d_cal, color="#1f77b4", lw=2.5, label="Model (calibrated)")
ax.plot(obs_days, obs_cum_deaths, "r^", ms=6, zorder=5,
        label="Observed cumulative deaths (OWID)")
ax.set_xlabel("Days since 1 June 2021", fontsize=11)
ax.set_ylabel("Cumulative deaths", fontsize=11)
ax.legend(fontsize=9)
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
plt.tight_layout()
plt.savefig(os.path.join(OUT, "fig_deaths_notitle.pdf"), bbox_inches="tight")
plt.close()
print("Saved fig_deaths_notitle.pdf")
