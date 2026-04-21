"""
COVID-19 SEIHM-R-D Model Simulation - White Background Version
Understanding COVID-19 Spread in Low-Income Countries
Model with Hospitalized and Home-Care Populations

Each scenario generates separate figure files with white backgrounds
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd

# Set style for white background and professional appearance
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['grid.color'] = 'gray'
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.linestyle'] = '--'

# Dark color palette for better visibility
DARK_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']


class COVIDModel:
    """
    SEIHM-R-D Compartmental Model for COVID-19
    
    Compartments:
    S - Susceptible
    E - Exposed
    I - Infected
    H - Hospitalized
    M - Home-care (self-medication)
    R - Recovered
    D - Deaths (cumulative)
    """
    
    def __init__(self, params):
        """Initialize model with baseline parameters from Table 1"""
        self.params = params
        
    def force_of_infection(self, S, I, M, H, N_L):
        """Calculate force of infection λ"""
        alpha = self.params['alpha']
        eta_M = self.params['eta_M']
        xi_H = self.params['xi_H']
        
        lambda_t = alpha * S * (I + eta_M * M + xi_H * H) / N_L
        return lambda_t
    
    def derivatives(self, y, t):
        """
        System of differential equations for SEIHM-R-D model
        """
        S, E, I, H, M, R, D = y
        
        # Living population
        N_L = S + E + I + H + M + R
        
        # Parameters
        Lambda = self.params['Lambda']
        alpha = self.params['alpha']
        beta = self.params['beta']
        gamma = self.params['gamma']
        kappa = self.params['kappa']
        theta_1 = self.params['theta_1']
        theta_2 = self.params['theta_2']
        pi_1 = self.params['pi_1']
        pi_2 = self.params['pi_2']
        phi_1 = self.params['phi_1']
        phi_2 = self.params['phi_2']
        mu_1 = self.params['mu_1']
        mu_2 = self.params['mu_2']
        epsilon_1 = self.params['epsilon_1']
        epsilon_2 = self.params['epsilon_2']
        
        # Composite rates
        omega_1 = mu_1 + epsilon_1 + phi_1 + pi_1 + theta_1
        omega_2 = mu_1 + epsilon_2 + phi_2 + pi_2 + theta_2
        
        # Force of infection
        lambda_t = self.force_of_infection(S, I, M, H, N_L)
        
        # Differential equations
        dS = Lambda - mu_1*S - lambda_t + epsilon_1*H + epsilon_2*M
        dE = lambda_t - (mu_1 + beta)*E
        dI = beta*E - (mu_1 + mu_2 + kappa)*I
        dH = kappa*gamma*I + theta_2*M - omega_1*H
        dM = kappa*(1-gamma)*I + theta_1*H - omega_2*M
        dR = phi_1*H + phi_2*M - mu_1*R
        dD = pi_1*H + pi_2*M + mu_2*I
        
        return [dS, dE, dI, dH, dM, dR, dD]
    
    def simulate(self, y0, t):
        """Simulate the model"""
        solution = odeint(self.derivatives, y0, t)
        return solution
    
    def calculate_R0(self):
        """Calculate basic reproduction number R0"""
        alpha = self.params['alpha']
        beta = self.params['beta']
        gamma = self.params['gamma']
        kappa = self.params['kappa']
        mu_1 = self.params['mu_1']
        mu_2 = self.params['mu_2']
        xi_H = self.params['xi_H']
        eta_M = self.params['eta_M']
        theta_1 = self.params['theta_1']
        theta_2 = self.params['theta_2']
        epsilon_1 = self.params['epsilon_1']
        epsilon_2 = self.params['epsilon_2']
        phi_1 = self.params['phi_1']
        phi_2 = self.params['phi_2']
        pi_1 = self.params['pi_1']
        pi_2 = self.params['pi_2']
        
        # Composite rates
        K = mu_1 + beta
        L = mu_1 + mu_2 + kappa
        omega_1 = mu_1 + epsilon_1 + phi_1 + pi_1 + theta_1
        omega_2 = mu_1 + epsilon_2 + phi_2 + pi_2 + theta_2
        D = omega_1 * omega_2 - theta_1 * theta_2
        
        # Transmission weights
        P_H = gamma * omega_2 + (1 - gamma) * theta_2
        P_M = gamma * theta_1 + (1 - gamma) * omega_1
        
        # R0 formula from Eq. (14)
        R0 = (alpha * beta * (D + kappa * xi_H * P_H + kappa * eta_M * P_M)) / (K * L * D)
        
        return R0


def get_baseline_parameters():
    """
    Baseline parameters from Table 1 in the paper
    """
    params = {
        'Lambda': 0.304,          # Recruitment rate
        'alpha': 0.035,           # Contact rate
        'beta': 0.095,            # Transmission rate
        'gamma': 0.29,            # Rate of infected people hospitalized
        'kappa': 0.100,           # Progression rate from I to care states
        'theta_1': 0.396,         # Rate from hospital to home
        'theta_2': 0.0396,        # Rate from home to hospital
        'pi_1': 0.0196,           # Mortality rate for hospitalized
        'pi_2': 0.009,            # Mortality rate for home-care
        'phi_1': 0.900,           # Recovery rate for hospitalized
        'phi_2': 0.0000625,       # Recovery rate for home-care
        'mu_1': 0.00576,          # Natural mortality rate
        'mu_2': 0.196,            # Disease-induced mortality from I
        'epsilon_1': 0.001,       # Rate from H to S
        'epsilon_2': 0.001,       # Rate from M to S
        'eta_M': 0.50,            # Relative infectiousness of home-care
        'xi_H': 0.10,             # Relative infectiousness of hospitalized
    }
    return params


def get_initial_conditions(N=1000):
    """
    Initial conditions for the model
    """
    S0 = N - 5
    E0 = 3
    I0 = 2
    H0 = 0
    M0 = 0
    R0 = 0
    D0 = 0
    
    return [S0, E0, I0, H0, M0, R0, D0]


def save_single_plot(fig, filename):
    """Helper function to save a single plot"""
    fig.savefig(f'/home/claude/{filename}', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  ✓ Saved: {filename}")


def scenario_1_self_medication():
    """
    Scenario 1: Impact of Self-Medication (varying gamma)
    Creates 4 separate figures
    """
    print("\n" + "="*60)
    print("SCENARIO 1: Impact of Self-Medication (varying γ)")
    print("="*60)
    
    t = np.linspace(0, 150, 1500)
    y0 = get_initial_conditions()
    
    # Three scenarios with different hospitalization rates
    scenarios = {
        'High Trust (γ=0.9)': 0.9,
        'Moderate (γ=0.5)': 0.5,
        'High Self-Med (γ=0.1)': 0.1
    }
    
    results = {}
    
    for scenario_name, gamma_val in scenarios.items():
        params = get_baseline_parameters()
        params['gamma'] = gamma_val
        
        model = COVIDModel(params)
        R0 = model.calculate_R0()
        sol = model.simulate(y0, t)
        
        results[scenario_name] = {
            'solution': sol,
            'R0': R0,
            'total_deaths': sol[-1, 6],
            'peak_infected': np.max(sol[:, 2]),
            'peak_homecare': np.max(sol[:, 4])
        }
        
        print(f"\n{scenario_name}:")
        print(f"  R0 = {R0:.3f}")
        print(f"  Total deaths (day 150) = {sol[-1, 6]:.1f}")
    
    # Figure 1a: Infected (I)
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 2], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Infected (I)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 1a: Infected Population by Hospitalization Rate', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario1a_infected.png')
    
    # Figure 1b: Home-care (M)
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 4], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Home-care (M)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 1b: Home-Care Patients (Self-Medication)', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario1b_homecare.png')
    
    # Figure 1c: Cumulative Deaths (D)
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 6], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Cumulative Deaths (D)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 1c: Cumulative Mortality by Hospitalization Rate', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario1c_deaths.png')
    
    # Figure 1d: Summary metrics
    fig, ax = plt.subplots(figsize=(10, 6))
    metrics_data = []
    for scenario_name, data in results.items():
        metrics_data.append([
            scenario_name,
            data['total_deaths'],
            data['peak_infected'],
            data['peak_homecare']
        ])
    
    df = pd.DataFrame(metrics_data, 
                     columns=['Scenario', 'Total Deaths', 'Peak Infected', 'Peak Home-Care'])
    
    x = np.arange(len(df))
    width = 0.25
    
    ax.bar(x - width, df['Total Deaths'], width, label='Total Deaths', 
           alpha=0.8, color='#d62728', edgecolor='black', linewidth=1.2)
    ax.bar(x, df['Peak Infected'], width, label='Peak Infected', 
           alpha=0.8, color='#1f77b4', edgecolor='black', linewidth=1.2)
    ax.bar(x + width, df['Peak Home-Care'], width, label='Peak Home-Care', 
           alpha=0.8, color='#2ca02c', edgecolor='black', linewidth=1.2)
    
    ax.set_xlabel('Scenario', fontsize=13, fontweight='bold')
    ax.set_ylabel('Number of People', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 1d: Key Metrics Comparison', fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(['High Trust\n(γ=0.9)', 'Moderate\n(γ=0.5)', 'High Self-Med\n(γ=0.1)'], 
                       fontsize=11)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3, axis='y')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario1d_metrics.png')
    
    return results


def scenario_2_community_surveillance():
    """
    Scenario 2: Community Surveillance Impact
    Creates 4 separate figures
    """
    print("\n" + "="*60)
    print("SCENARIO 2: Community Surveillance (varying ηM and θ2)")
    print("="*60)
    
    t = np.linspace(0, 150, 1500)
    y0 = get_initial_conditions()
    
    scenarios = {
        'No Surveillance': {'eta_M': 0.8, 'theta_2': 0.01},
        'Moderate Surveillance': {'eta_M': 0.5, 'theta_2': 0.04},
        'Intensive Surveillance': {'eta_M': 0.2, 'theta_2': 0.10}
    }
    
    results = {}
    
    for scenario_name, scenario_params in scenarios.items():
        params = get_baseline_parameters()
        params.update(scenario_params)
        
        model = COVIDModel(params)
        R0 = model.calculate_R0()
        sol = model.simulate(y0, t)
        
        results[scenario_name] = {
            'solution': sol,
            'R0': R0,
            'total_deaths': sol[-1, 6],
        }
        
        print(f"\n{scenario_name}:")
        print(f"  R0 = {R0:.3f}, Total deaths = {sol[-1, 6]:.1f}")
    
    # Figure 2a: Home-care (M)
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 4], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Home-care (M)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 2a: Home-Care Patients by Surveillance Level', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario2a_homecare.png')
    
    # Figure 2b: Hospitalized (H)
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 3], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Hospitalized (H)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 2b: Hospitalized Patients by Surveillance Level', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario2b_hospitalized.png')
    
    # Figure 2c: Infected (I)
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 2], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Infected (I)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 2c: Active Infections by Surveillance Level', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario2c_infected.png')
    
    # Figure 2d: Deaths (D)
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 6], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Cumulative Deaths (D)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 2d: Cumulative Mortality by Surveillance Level', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario2d_deaths.png')
    
    return results


def scenario_3_hospital_capacity():
    """
    Scenario 3: Hospital Capacity and Early Discharge
    Creates 4 separate figures
    """
    print("\n" + "="*60)
    print("SCENARIO 3: Hospital Capacity (varying θ1)")
    print("="*60)
    
    t = np.linspace(0, 150, 1500)
    y0 = get_initial_conditions()
    
    scenarios = {
        'Sufficient Capacity': 0.05,
        'Overcrowded': 0.30,
        'Critical Overflow': 0.50
    }
    
    results = {}
    
    for scenario_name, theta_1_val in scenarios.items():
        params = get_baseline_parameters()
        params['theta_1'] = theta_1_val
        
        model = COVIDModel(params)
        R0 = model.calculate_R0()
        sol = model.simulate(y0, t)
        
        results[scenario_name] = {
            'solution': sol,
            'R0': R0,
            'total_deaths': sol[-1, 6],
        }
        
        print(f"\n{scenario_name}: R0 = {R0:.3f}, Deaths = {sol[-1, 6]:.1f}")
    
    # Figure 3a: Hospitalized
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 3], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Hospitalized (H)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 3a: Hospitalized Patients by Capacity Level', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario3a_hospitalized.png')
    
    # Figure 3b: Home-care
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 4], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Home-care (M)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 3b: Home-Care Patients by Capacity Level', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario3b_homecare.png')
    
    # Figure 3c: Cumulative Deaths
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 6], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Cumulative Deaths (D)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 3c: Cumulative Mortality by Capacity Level', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario3c_deaths.png')
    
    # Figure 3d: Ratio H/M
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        H = data['solution'][:, 3]
        M = data['solution'][:, 4]
        ratio = np.divide(H, M + 1e-10)
        ax.plot(t, ratio, label=scenario_name, linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('H/M Ratio', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 3d: Hospital to Home-Care Ratio', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario3d_ratio.png')
    
    return results


def scenario_4_integrated_policy():
    """
    Scenario 4: Integrated Policy Package
    Creates 4 separate figures
    """
    print("\n" + "="*60)
    print("SCENARIO 4: Integrated Policy Packages")
    print("="*60)
    
    t = np.linspace(0, 150, 1500)
    y0 = get_initial_conditions()
    
    baseline_params = get_baseline_parameters()
    baseline_params.update({
        'gamma': 0.2,
        'eta_M': 0.5,
        'theta_1': 0.4,
        'theta_2': 0.02
    })
    
    scenarios = {
        'Status Quo': baseline_params.copy(),
        'Awareness Campaign': baseline_params.copy(),
        'Community Surveillance': baseline_params.copy(),
        'Integrated Approach': baseline_params.copy()
    }
    
    scenarios['Awareness Campaign']['gamma'] = 0.6
    scenarios['Community Surveillance']['eta_M'] = 0.2
    scenarios['Community Surveillance']['theta_2'] = 0.08
    scenarios['Integrated Approach']['gamma'] = 0.6
    scenarios['Integrated Approach']['eta_M'] = 0.2
    scenarios['Integrated Approach']['theta_1'] = 0.15
    scenarios['Integrated Approach']['theta_2'] = 0.08
    
    results = {}
    
    for scenario_name, params in scenarios.items():
        model = COVIDModel(params)
        R0 = model.calculate_R0()
        sol = model.simulate(y0, t)
        
        results[scenario_name] = {
            'solution': sol,
            'R0': R0,
            'total_deaths': sol[-1, 6],
            'peak_infected': np.max(sol[:, 2]),
        }
        
        print(f"\n{scenario_name}: R0 = {R0:.3f}, Deaths = {sol[-1, 6]:.1f}")
    
    # Figure 4a: All compartments for Status Quo
    fig, ax = plt.subplots(figsize=(10, 6))
    sol = results['Status Quo']['solution']
    compartments = [
        ('S', 0, '#1f77b4'),
        ('E', 1, '#ff7f0e'),
        ('I', 2, '#2ca02c'),
        ('H', 3, '#d62728'),
        ('M', 4, '#9467bd'),
        ('R', 5, '#8c564b')
    ]
    for label, idx, color in compartments:
        ax.plot(t, sol[:, idx], label=label, linewidth=2.5, color=color)
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Number of People', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 4a: Status Quo - All Compartments', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario4a_status_quo.png')
    
    # Figure 4b: Infected comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 2], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Infected (I)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 4b: Active Infections by Policy', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario4b_infected.png')
    
    # Figure 4c: Cumulative deaths
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (scenario_name, data) in enumerate(results.items()):
        ax.plot(t, data['solution'][:, 6], label=scenario_name, 
                linewidth=2.5, color=DARK_COLORS[idx])
    ax.set_xlabel('Time (days)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Cumulative Deaths (D)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 4c: Cumulative Mortality by Policy', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=11, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_single_plot(fig, 'scenario4c_deaths.png')
    
    # Figure 4d: Summary bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    scenario_names = list(results.keys())
    deaths = [results[s]['total_deaths'] for s in scenario_names]
    
    bars = ax.bar(range(len(scenario_names)), deaths, 
                   color=DARK_COLORS[:len(scenario_names)], 
                   alpha=0.8, edgecolor='black', linewidth=1.5)
    ax.set_xlabel('Policy Scenario', fontsize=13, fontweight='bold')
    ax.set_ylabel('Total Deaths (Day 150)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 4d: Mortality by Policy Scenario', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(range(len(scenario_names)))
    ax.set_xticklabels(['Status\nQuo', 'Awareness\nCampaign', 
                        'Community\nSurveillance', 'Integrated\nApproach'], 
                       fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    for bar, death_count in zip(bars, deaths):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{death_count:.1f}',
                ha='center', va='bottom', fontweight='bold', fontsize=11)
    
    save_single_plot(fig, 'scenario4d_comparison.png')
    
    return results


def scenario_5_sensitivity_heatmap():
    """
    Scenario 5: Sensitivity Analysis Heatmap
    Creates 2 separate heatmap figures
    """
    print("\n" + "="*60)
    print("SCENARIO 5: Sensitivity Analysis Heatmap")
    print("="*60)
    
    t = np.linspace(0, 150, 1500)
    y0 = get_initial_conditions()
    
    gamma_range = np.linspace(0.1, 0.9, 15)
    eta_M_range = np.linspace(0.1, 0.9, 15)
    
    mortality_matrix = np.zeros((len(gamma_range), len(eta_M_range)))
    peak_infected_matrix = np.zeros((len(gamma_range), len(eta_M_range)))
    
    print(f"\nRunning {len(gamma_range) * len(eta_M_range)} simulations...")
    
    for i, gamma_val in enumerate(gamma_range):
        for j, eta_M_val in enumerate(eta_M_range):
            params = get_baseline_parameters()
            params['gamma'] = gamma_val
            params['eta_M'] = eta_M_val
            
            model = COVIDModel(params)
            sol = model.simulate(y0, t)
            
            mortality_matrix[i, j] = sol[-1, 6]
            peak_infected_matrix[i, j] = np.max(sol[:, 2])
    
    print("✓ Simulations complete!")
    
    # Figure 5a: Mortality heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(mortality_matrix, cmap='YlOrRd', aspect='auto', origin='lower',
                   extent=[eta_M_range[0], eta_M_range[-1], gamma_range[0], gamma_range[-1]])
    ax.set_xlabel('Home-care Infectiousness (ηM)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Hospitalization Rate (γ)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 5a: Total Mortality Sensitivity (Day 150)', 
                 fontsize=14, fontweight='bold', pad=20)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Total Deaths', fontsize=12, fontweight='bold')
    
    CS = ax.contour(eta_M_range, gamma_range, mortality_matrix, 
                    levels=8, colors='black', alpha=0.4, linewidths=1.5)
    ax.clabel(CS, inline=True, fontsize=9)
    
    save_single_plot(fig, 'scenario5a_mortality_heatmap.png')
    
    # Figure 5b: Peak infections heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(peak_infected_matrix, cmap='Blues', aspect='auto', origin='lower',
                   extent=[eta_M_range[0], eta_M_range[-1], gamma_range[0], gamma_range[-1]])
    ax.set_xlabel('Home-care Infectiousness (ηM)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Hospitalization Rate (γ)', fontsize=13, fontweight='bold')
    ax.set_title('Scenario 5b: Peak Infections Sensitivity', 
                 fontsize=14, fontweight='bold', pad=20)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Peak Infected', fontsize=12, fontweight='bold')
    
    CS = ax.contour(eta_M_range, gamma_range, peak_infected_matrix, 
                    levels=8, colors='black', alpha=0.4, linewidths=1.5)
    ax.clabel(CS, inline=True, fontsize=9)
    
    save_single_plot(fig, 'scenario5b_peak_infections_heatmap.png')
    
    min_idx = np.unravel_index(np.argmin(mortality_matrix), mortality_matrix.shape)
    max_idx = np.unravel_index(np.argmax(mortality_matrix), mortality_matrix.shape)
    
    print(f"\nMortality range: {mortality_matrix[min_idx]:.1f} to {mortality_matrix[max_idx]:.1f}")
    
    return mortality_matrix, peak_infected_matrix


def create_summary_table(all_results):
    """Create summary table"""
    print("\n" + "="*60)
    print("SUMMARY TABLE")
    print("="*60)
    
    summary_data = []
    
    for scenario_name, data in all_results['scenario1'].items():
        summary_data.append({
            'Category': 'Self-Medication',
            'Scenario': scenario_name,
            'R0': data['R0'],
            'Total Deaths': data['total_deaths'],
        })
    
    for scenario_name, data in all_results['scenario4'].items():
        summary_data.append({
            'Category': 'Policy Package',
            'Scenario': scenario_name,
            'R0': data['R0'],
            'Total Deaths': data['total_deaths'],
        })
    
    df = pd.DataFrame(summary_data)
    df.to_csv('/home/claude/summary_table.csv', index=False)
    print("\n✓ Summary table saved: summary_table.csv\n")
    print(df.to_string(index=False))
    
    return df


def main():
    """Main function"""
    print("\n" + "="*70)
    print("COVID-19 SEIHM-R-D MODEL - WHITE BACKGROUND VERSION")
    print("="*70)
    
    all_results = {}
    
    all_results['scenario1'] = scenario_1_self_medication()
    all_results['scenario2'] = scenario_2_community_surveillance()
    all_results['scenario3'] = scenario_3_hospital_capacity()
    all_results['scenario4'] = scenario_4_integrated_policy()
    
    scenario_5_sensitivity_heatmap()
    
    summary_df = create_summary_table(all_results)
    
    print("\n" + "="*70)
    print("COMPLETE! 18 figures generated (white background, dark colors)")
    print("="*70)
    
    return all_results, summary_df


if __name__ == "__main__":
    results, summary = main()
