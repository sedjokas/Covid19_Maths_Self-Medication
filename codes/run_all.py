import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import BASELINE, R0, rhs, simulate, default_ic, center_manifold_a, alpha_for_R0, derived

plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'scenario_figures')
os.makedirs(OUT, exist_ok=True)

def report(title):
    print("\n" + "="*70)
    print(title)
    print("="*70)

# ---------- Scenario 1: Baseline ----------
report("SCENARIO 1: Baseline")
t = np.linspace(0, 150, 1500)
y0 = default_ic(1000)
sol = simulate(BASELINE, t, y0)
S,E,I,H,M,R,D = sol.T
N_L = S+E+I+H+M+R
pk = np.argmax(M)
print(f"R0={R0(BASELINE):.4f}")
print(f"peak M={M.max():.2f}/1000 ({M[pk]/N_L[pk]*100:.1f}%), peak H at same day={H[pk]:.2f}/1000 ({H[pk]/N_L[pk]*100:.1f}%), ratio={M[pk]/H[pk]:.2f}:1")
print(f"cumulative deaths day150={D[-1]:.2f}/1000")
fig, axes = plt.subplots(1,2, figsize=(11,4.5))
for name,arr in [('S',S),('E',E),('I',I),('H',H),('M',M),('R',R)]:
    axes[0].plot(t, arr, label=name)
axes[0].legend(fontsize=8); axes[0].set_xlabel('Time (days)'); axes[0].set_ylabel('Population per 1,000'); axes[0].set_title(f'Baseline dynamics, R0={R0(BASELINE):.2f}')
axes[1].plot(t, D, color='black'); axes[1].set_xlabel('Time (days)'); axes[1].set_ylabel('Cumulative deaths per 1,000'); axes[1].set_title('Cumulative mortality')
fig.tight_layout(); fig.savefig(f"{OUT}/scenario1_baseline.png", dpi=200); plt.close(fig)

# ---------- Scenario 2: Care-seeking (vary gamma) ----------
report("SCENARIO 2: Care-seeking (gamma)")
fig, ax = plt.subplots(figsize=(6,4.5))
results2 = {}
for label, g in [('High Trust',0.9), ('Moderate',0.5), ('Self-medication',0.1)]:
    p = dict(BASELINE, gamma=g)
    sol = simulate(p, t, y0)
    S,E,I,H,M,R,D = sol.T
    results2[label] = (R0(p), D[-1])
    ax.plot(t, D, label=f"{label} (γ={g})")
    print(f"{label}: gamma={g} R0={R0(p):.4f} deaths_day150={D[-1]:.2f}/1000")
ax.legend(fontsize=8); ax.set_xlabel('Time (days)'); ax.set_ylabel('Cumulative deaths per 1,000'); ax.set_title('Effect of hospitalization proportion gamma')
fig.tight_layout(); fig.savefig(f"{OUT}/scenario2_gamma.png", dpi=200); plt.close(fig)

# ---------- Scenario 3: Contact-rate sweep ----------
report("SCENARIO 3: Contact-rate sweep")
alpha_lo = alpha_for_R0(BASELINE, 0.5)
alpha_hi = alpha_for_R0(BASELINE, 4.0)
print(f"alpha range for R0 0.5->4.0: {alpha_lo:.5f} -> {alpha_hi:.5f} (fold={alpha_hi/alpha_lo:.2f})")
fig, axes = plt.subplots(1,2, figsize=(11,4.5))
peak_I=[]; deaths150=[]; alphas=np.linspace(alpha_lo, alpha_hi, 8)
for a in alphas:
    p = dict(BASELINE, alpha=a)
    sol = simulate(p, t, y0)
    S,E,I,H,M,R,D = sol.T
    peak_I.append(I.max()); deaths150.append(D[-1])
    print(f"alpha={a:.5f} R0={R0(p):.3f} peakI={I.max():.3f} deaths150={D[-1]:.2f}")
axes[0].plot([R0(dict(BASELINE,alpha=a)) for a in alphas], peak_I, 'o-')
axes[0].axvline(1, color='red', ls='--', lw=1)
axes[0].set_xlabel('R0'); axes[0].set_ylabel('Peak infected I per 1,000'); axes[0].set_title('Peak I vs R0')
axes[1].plot([R0(dict(BASELINE,alpha=a)) for a in alphas], deaths150, 'o-', color='darkred')
axes[1].axvline(1, color='red', ls='--', lw=1)
axes[1].set_xlabel('R0'); axes[1].set_ylabel('Cumulative deaths, day 150 per 1,000'); axes[1].set_title('Mortality vs R0')
fig.tight_layout(); fig.savefig(f"{OUT}/scenario3_contact_sweep.png", dpi=200); plt.close(fig)
print(f"lowest alpha deaths={deaths150[0]:.3f} highest alpha deaths={deaths150[-1]:.3f} ratio={deaths150[-1]/deaths150[0]:.2f}")
print(f"lowest alpha peakI={peak_I[0]:.3f} highest alpha peakI={peak_I[-1]:.3f} ratio={peak_I[-1]/peak_I[0]:.2f}")

# ---------- Scenario 4: Bifurcation verification ----------
report("SCENARIO 4: Bifurcation verification (sweep phi_1)")
phi1_range = np.linspace(0.002, 0.095, 60)
eps_values = [0.0, 0.001, 0.02, 0.10, 0.30]
fig, ax = plt.subplots(figsize=(7,5))
amax_overall = -np.inf
for eps in eps_values:
    avals=[]
    for phi1 in phi1_range:
        p = dict(BASELINE, phi_1=phi1, epsilon_1=eps, epsilon_2=eps)
        _, a = center_manifold_a(p)
        avals.append(a)
    amax_overall = max(amax_overall, max(avals))
    ax.plot(phi1_range, avals, label=f"eps1=eps2={eps}")
ax.axhline(0, color='black', lw=1)
ax.set_xlabel('Hospital recovery rate phi_1'); ax.set_ylabel('center-manifold coefficient a')
ax.legend(fontsize=8); ax.set_title('a is negative throughout (forward bifurcation only)')
fig.tight_layout(); fig.savefig(f"{OUT}/scenario4_bifurcation_a.png", dpi=200); plt.close(fig)
print(f"max(a) over sweep = {amax_overall:.6e} (should be <0)")

# ---------- Scenario 5: No-bistability, long horizon near R0=1 ----------
report("SCENARIO 5: No-bistability demonstration")
alpha_near1 = alpha_for_R0(BASELINE, 0.988)
p5 = dict(BASELINE, alpha=alpha_near1)
print(f"alpha for R0~0.988: {alpha_near1:.6f}  R0={R0(p5):.4f}")
tlong = np.linspace(0, 60000, 3000)
ics = {
    'Minuscule': [999.9,0.05,0.05,0,0,0,0],
    'Tiny':      [999.0,0.5,0.5,0,0,0,0],
    'Small':     [995.0,3,2,0,0,0,0],
    'Medium':    [950,25,25,0,0,0,0],
    'Large':     [900,50,50,0,0,0,0],
}
fig, ax = plt.subplots(figsize=(7,5))
finalS=[]
for name,y0i in ics.items():
    sol = simulate(p5, tlong, y0i)
    ax.plot(tlong, sol[:,0], label=name)
    finalS.append(sol[-1,0])
ax.legend(fontsize=8); ax.set_xlabel('Time (days)'); ax.set_ylabel('Susceptible S')
ax.set_title(f'All trajectories converge to the same DFE (R0={R0(p5):.3f})')
fig.tight_layout(); fig.savefig(f"{OUT}/scenario5_no_bistability.png", dpi=200); plt.close(fig)
print(f"spread in final S at t=60000: {max(finalS)-min(finalS):.6f}")

# ---------- Scenario 6: Hospital discharge policy (theta_1) ----------
report("SCENARIO 6: Hospital discharge policy (theta_1)")
base_theta1 = BASELINE['theta_1']
fig, ax = plt.subplots(figsize=(6,4.5))
d6=[]
for label, th1 in [('Sufficient (low discharge)',base_theta1*0.5), ('Baseline',base_theta1), ('Overcrowded',base_theta1*3), ('Critical-overflow',base_theta1*6)]:
    p = dict(BASELINE, theta_1=th1)
    sol = simulate(p, t, y0)
    D = sol[:,6]
    d6.append((label, th1, D[-1]))
    ax.plot(t, D, label=f"{label} (θ1={th1:.4f})")
    print(f"{label}: theta1={th1:.4f} deaths150={D[-1]:.3f}")
ax.legend(fontsize=7); ax.set_xlabel('Time (days)'); ax.set_ylabel('Cumulative deaths per 1,000'); ax.set_title('Effect of hospital discharge rate theta_1')
fig.tight_layout(); fig.savefig(f"{OUT}/scenario6_discharge.png", dpi=200); plt.close(fig)
pct_increase = (d6[3][2]-d6[0][2])/d6[0][2]*100
print(f"pct increase critical-overflow vs sufficient: {pct_increase:.1f}%")

# ---------- Scenario 7: Joint sensitivity gamma, eta_M ----------
report("SCENARIO 7: Joint sensitivity of mortality to gamma, eta_M")
gamma_range = np.linspace(0.1,0.9,9)
etaM_range = np.linspace(0.1,0.9,9)
mort = np.zeros((len(gamma_range), len(etaM_range)))
for i,g in enumerate(gamma_range):
    for j,em in enumerate(etaM_range):
        p = dict(BASELINE, gamma=g, eta_M=em)
        sol = simulate(p, t, y0)
        mort[i,j] = sol[-1,6]
fig, ax = plt.subplots(figsize=(6,5))
c = ax.contourf(etaM_range, gamma_range, mort, levels=20, cmap='viridis')
fig.colorbar(c, label='Cumulative deaths per 1,000, day150')
ax.set_xlabel('eta_M'); ax.set_ylabel('gamma'); ax.set_title('Joint sensitivity of mortality')
fig.tight_layout(); fig.savefig(f"{OUT}/scenario7_joint_sensitivity.png", dpi=200); plt.close(fig)
print(f"mortality range: {mort.min():.2f} to {mort.max():.2f} per 1000")
print(f"raising eta_M 0.1->0.9 at gamma=0.5: {mort[4,0]:.2f} -> {mort[4,-1]:.2f}")
print(f"raising gamma 0.1->0.9 at eta_M=0.5: {mort[0,4]:.2f} -> {mort[-1,4]:.2f}")

# ---------- Scenario 8: Community surveillance ----------
report("SCENARIO 8: Community surveillance")
configs8 = [('No surveillance', 0.8, 0.005), ('Moderate', 0.5, BASELINE['theta_2']), ('Intensive', 0.2, BASELINE['theta_2']*3)]
res8=[]
for label, em, th2 in configs8:
    p = dict(BASELINE, eta_M=em, theta_2=th2)
    sol = simulate(p, t, y0)
    S,E,I,H,M,R,D = sol.T
    res8.append((label, em, th2, R0(p), D[-1]))
    print(f"{label}: eta_M={em} theta2={th2:.4f} R0={R0(p):.4f} deaths150={D[-1]:.2f}")
r0_reduction = (res8[0][3]-res8[2][3])/res8[0][3]*100
deaths_reduction = (res8[0][4]-res8[2][4])/res8[0][4]*100
print(f"R0 reduction intensive vs none: {r0_reduction:.1f}%  deaths reduction: {deaths_reduction:.1f}%")

# ---------- Scenario 9: Integrated policy ----------
report("SCENARIO 9: Integrated policy comparison")
p_base9 = dict(BASELINE, gamma=0.2, eta_M=0.5, theta_2=0.02)
configs9 = [
    ('Status quo', dict(p_base9)),
    ('Awareness campaign', dict(p_base9, gamma=0.6)),
    ('Community surveillance', dict(p_base9, eta_M=0.2, theta_2=0.08)),
    ('Integrated', dict(p_base9, gamma=0.6, eta_M=0.2, theta_2=0.08)),
]
res9=[]
for label, p in configs9:
    sol = simulate(p, t, y0)
    D = sol[:,6]
    res9.append((label, R0(p), D[-1]))
    print(f"{label}: R0={R0(p):.4f} deaths150={D[-1]:.3f}")
base_d = res9[0][2]
for label, r0v, dv in res9[1:]:
    print(f"  {label}: deaths averted={base_d-dv:.3f} ({(base_d-dv)/base_d*100:.1f}%), R0 reduction={(res9[0][1]-r0v)/res9[0][1]*100:.1f}%")

# ---------- Scenario 10: Convergence to endemic equilibrium (phase plane) ----------
report("SCENARIO 10: Convergence to endemic equilibrium")
tlong2 = np.linspace(0, 1500, 3000)
ics10 = [
    [999,0.5,0.5,0,0,0,0],
    [900,50,50,0,0,0,0],
    [700,150,100,50,0,0,0],
    [500,200,150,100,50,0,0],
]
fig, ax = plt.subplots(figsize=(6,5))
for y0i in ics10:
    sol = simulate(BASELINE, tlong2, y0i)
    ax.plot(sol[:,0], sol[:,2], alpha=0.7)
    ax.plot(sol[0,0], sol[0,2], 'o', color='green', ms=4)
sol_eq = simulate(BASELINE, np.linspace(0,3000,2000), y0)
ax.plot(sol_eq[-1,0], sol_eq[-1,2], '*', color='red', ms=15, label='Endemic equilibrium X*')
ax.set_xlabel('S'); ax.set_ylabel('I'); ax.legend(); ax.set_title(f'Phase plane, R0={R0(BASELINE):.2f}: convergence to X*')
fig.tight_layout(); fig.savefig(f"{OUT}/scenario10_phase.png", dpi=200); plt.close(fig)
print(f"Endemic equilibrium (S*,I*): ({sol_eq[-1,0]:.3f}, {sol_eq[-1,2]:.3f})")

print("\nALL DONE")
