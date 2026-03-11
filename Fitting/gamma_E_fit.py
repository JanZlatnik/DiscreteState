import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from scipy.optimize import curve_fit
import scipy.special as sp 
import gc
import os
import sys
import argparse

script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.append(script_dir)
import settings_fit

parser = argparse.ArgumentParser(description="Script for fitting complex F(E).")
parser.add_argument('--model', type=str, default='2DModel2', help='Select the calculation mode defined in settings_fit.py')
args = parser.parse_args()

if args.model not in settings_fit.MODELS:
    print(f"\n[ERROR] Mode '{args.model}' is not defined in settings_fit.py!")
    sys.exit(1)

cfg = settings_fit.MODELS[args.model]
print(f"--> Running {os.path.basename(__file__)} in mode: {args.model}")

R0_target = cfg.get('Fit_R0', 2.0)

base_path = cfg.get('base_path', r"./DATA/DSState")
output_path = cfg.get('output_path', r"./Graphs/Fitting")
os.makedirs(output_path, exist_ok=True)
aspect_ratio = cfg.get('aspect_ratio', 'square')
SIZE = 6
SMALL_SIZE = 14; MEDIUM_SIZE = 18; BIGGER_SIZE = 24; LEGEND_SIZE = 12
X_axis = cfg.get('X_axis_E', r'Energy $(\mathrm{eV})$')
Y_axis_gamma = cfg.get('Y_axis_gamma', r'$\Gamma(R_0, E)\, (\mathrm{eV})$')
Y_axis_delta = cfg.get('Y_axis_delta', r'$\Delta(R_0, E)\, (\mathrm{eV})$')
color_data = cfg.get('color_data', 'black')
color_fit = cfg.get('color_fit', 'red')
y_lim_delta = cfg.get('y_lim_delta', None)

plt.rcParams['font.family'] = 'serif'
plt.rcParams["font.serif"] = ["Latin Modern Roman"]
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['axes.titlepad'] = 10; plt.rcParams['axes.labelpad'] = 10 
plt.rcParams['legend.fancybox'] = False; plt.rcParams['legend.edgecolor'] = "#000000"
plt.rcParams["figure.autolayout"] = True
plt.rc('font', size=SMALL_SIZE); plt.rc('axes', titlesize=MEDIUM_SIZE) 
plt.rc('axes', labelsize=MEDIUM_SIZE); plt.rc('xtick', labelsize=SMALL_SIZE) 
plt.rc('ytick', labelsize=SMALL_SIZE); plt.rc('legend', fontsize=LEGEND_SIZE) 
legend_params = {'frameon': False, 'ncol': 1, 'fontsize': 12, 'handlelength': 2, 'loc': 'best'}

data_ds = np.loadtxt(f"{base_path}/DSenergies.txt", comments="#")
R_vals = data_ds[:, 0]
idx_R0 = np.argmin(np.abs(R_vals - R0_target))
R0_actual = R_vals[idx_R0]

data_gamma = np.loadtxt(f"{base_path}/gamma.txt", comments="#")
E_vals_gamma = data_gamma[:, 0]
Gamma_R0 = data_gamma[:, idx_R0 + 1]

data_delta = np.loadtxt(f"{base_path}/delta.txt", comments="#")
E_vals_delta = data_delta[:, 0]
Delta_R0 = data_delta[:, idx_R0 + 1]

base_path_sparse = cfg.get('base_path_sparse', r"./DATA/DSStateSparce")
try:
    data_gamma_sp = np.loadtxt(f"{base_path_sparse}/gamma.txt", comments="#")
    E_vals_gamma_sp = data_gamma_sp[:, 0]
    Gamma_R0_sp = data_gamma_sp[:, idx_R0 + 1]

    data_delta_sp = np.loadtxt(f"{base_path_sparse}/delta.txt", comments="#")
    E_vals_delta_sp = data_delta_sp[:, 0]
    Delta_R0_sp = data_delta_sp[:, idx_R0 + 1]
    has_sparse = True
except OSError:
    print(f"[WARNING] Sparse data not found at {base_path_sparse}. Continuing without extrapolation check.")
    has_sparse = False

idx_E_pos = np.where(E_vals_gamma > 0.0)[0]
E_pos = E_vals_gamma[idx_E_pos]
Gamma_pos = Gamma_R0[idx_E_pos]
Delta_pos = Delta_R0[idx_E_pos]

A0_fixed = Gamma_pos[0]
c0_fixed = Delta_pos[0] * (2.0 * np.pi / A0_fixed)

def calc_defect(E_eV):
    return -2.17015937 / (1 -0.03881984 * E_eV + 0.00469617 * E_eV**2)

def F_complex(E_eV, A1, A2, A3, A4, A5, A6, A7, c0, c1, c2, c3, c4, c5):
    E_au = E_eV / 27.211386  
    k = np.sqrt(2 * E_au + 0j)
    
    mu = calc_defect(E_eV)
    mu_eff = np.where(E_eV < 0, mu, 0.0)
    
    A_E = np.where(E_eV > 0, (A1 * E_eV + A2 * E_eV ** 2) / (1 + (E_eV / A3)**4) + A4 * np.exp(-((E_eV - A5) / A6)**2), 0.0)
    
    c_E = (c1 * E_eV + c2 * E_eV**2) / (c0 + c3 * E_eV + c4 * E_eV**2 + c5 * E_eV**4)
    
    coulomb_damp = np.where(E_eV > 0, np.exp(-(E_eV / A7)**2), 1.0)
    
    term1 = c_E
    term2 = -0.5 * 1j * k
    term3 = np.log(-1j * k)
    term4 = sp.digamma(1.0 - 1j / k - mu_eff)
    
    F_val = (A0_fixed / (2.0 * np.pi)) * (c0_fixed  + c_E - A_E * 1j + (term2 + term3 + term4) * coulomb_damp)
    return F_val

def fit_gamma(E_eV, A1, A2, A3, A4, A5, A6, A7):
    return -2.0 * np.imag(F_complex(E_eV, A1, A2, A3, A4, A5, A6, A7, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0))

if has_sparse:
    sparse_E_min = cfg.get('fit_sparse_E_min', -50.0)
    sparse_E_max = cfg.get('fit_sparse_E_max', 45.0) 
    
    mask_sp_gamma = (E_vals_gamma_sp >= sparse_E_min) & (E_vals_gamma_sp <= sparse_E_max) & (E_vals_gamma_sp > 0)
    E_vals_gamma_sp_filtered = E_vals_gamma_sp[mask_sp_gamma]
    Gamma_R0_sp_filtered = Gamma_R0_sp[mask_sp_gamma]

    sigma_dense_gamma = np.ones_like(E_pos) * 0.04
    sigma_sparse_gamma = np.ones_like(E_vals_gamma_sp_filtered) * 2.0

    E_all_gamma = np.concatenate([E_pos, E_vals_gamma_sp_filtered])
    Gamma_all_data = np.concatenate([Gamma_pos, Gamma_R0_sp_filtered])
    sigma_all_gamma = np.concatenate([sigma_dense_gamma, sigma_sparse_gamma])
    
    E_all_gamma, unique_idx_g = np.unique(E_all_gamma, return_index=True)
    Gamma_all_data = Gamma_all_data[unique_idx_g]
    sigma_weights_gamma = sigma_all_gamma[unique_idx_g]
else:
    E_all_gamma = E_pos
    Gamma_all_data = Gamma_pos
    sigma_weights_gamma = np.ones_like(E_pos) * 0.1

p0_gamma = cfg.get('p0_guess_gamma', [0.01, 0.005, 5.0, 2.5, 13.0, 4.0, 20.0])

try:
    popt_gamma, pcov_gamma = curve_fit(
        fit_gamma, 
        E_all_gamma, 
        Gamma_all_data, 
        p0=p0_gamma,
        bounds=([-np.inf, -np.inf, 0.1, 0.0, 5.0, 0.1, 1.0], [np.inf, np.inf, np.inf, np.inf, 40.0, np.inf, np.inf]),
        sigma=sigma_weights_gamma,
        absolute_sigma=False
    )
    A1_fit, A2_fit, A3_fit, A4_fit, A5_fit, A6_fit, A7_fit = popt_gamma
    print(f"--> Gamma Fit converged. Parameters [A1...A7]: {popt_gamma}")
except RuntimeError as e:
    print(f"[ERROR] Optimal parameters not found for Gamma fit: {e}")
    sys.exit(1)

def calc_delta(E_eV, A1, A2, A3, A4, A5, A6, A7, c0, c1, c2, c3, c4, c5):
    return np.real(F_complex(E_eV, A1, A2, A3, A4, A5, A6, A7, c0, c1, c2, c3, c4, c5))

def fit_delta_c(E_eV, c0, c1, c2, c3, c4, c5):
    return calc_delta(E_eV, A1_fit, A2_fit, A3_fit, A4_fit, A5_fit, A6_fit, A7_fit, c0, c1, c2, c3, c4, c5) 

if has_sparse:
    sparse_E_min = cfg.get('fit_sparse_E_min', -50.0)
    sparse_E_max = cfg.get('fit_sparse_E_max', 20.0)
    
    mask_sp = (E_vals_delta_sp >= sparse_E_min) & (E_vals_delta_sp <= sparse_E_max)
    E_vals_delta_sp_filtered = E_vals_delta_sp[mask_sp]
    Delta_R0_sp_filtered = Delta_R0_sp[mask_sp]

    sigma_dense = np.where(E_vals_delta > 0, 0.05, 1)
    sigma_sparse = np.where(E_vals_delta_sp_filtered > 0, 2.0, 0.25) 

    E_all_delta = np.concatenate([E_vals_delta, E_vals_delta_sp_filtered])
    Delta_all_data = np.concatenate([Delta_R0, Delta_R0_sp_filtered])
    sigma_all = np.concatenate([sigma_dense, sigma_sparse])
    
    E_all_delta, unique_idx = np.unique(E_all_delta, return_index=True)
    Delta_all_data = Delta_all_data[unique_idx]
    sigma_weights = sigma_all[unique_idx]
else:
    E_all_delta = E_vals_delta
    Delta_all_data = Delta_R0
    sigma_weights = np.where(E_all_delta > 0, 0.1, 1.0)

try:
    popt_delta, pcov_delta = curve_fit(
        fit_delta_c, 
        E_all_delta, 
        Delta_all_data, 
        p0=[1.0, 0.0, 0.0, 0.1, 0.01, 0.001],
        sigma=sigma_weights,
        absolute_sigma=False
    )
    c0_fit, c1_fit, c2_fit, c3_fit, c4_fit, c5_fit = popt_delta
    print(f"--> Delta Fit converged. Parameters [c0...c5]: {popt_delta}")
except RuntimeError:
    print("[WARNING] Optimal c not found, using default zeros")
    c0_fit, c1_fit, c2_fit, c3_fit, c4_fit, c5_fit = 1.0, 0.0, 0.0, 0.1, 0.01, 0.001

E_dense_pos = np.linspace(min(E_pos), max(E_pos), 500)
Gamma_fit = fit_gamma(E_dense_pos, A1_fit, A2_fit, A3_fit, A4_fit, A5_fit, A6_fit, A7_fit)
Delta_fit_pos = calc_delta(E_dense_pos, A1_fit, A2_fit, A3_fit, A4_fit, A5_fit, A6_fit, A7_fit, c0_fit, c1_fit, c2_fit, c3_fit, c4_fit, c5_fit)

E_dense_short = np.linspace(min(E_vals_delta), max(E_vals_delta), 5000)
Delta_fit_short = calc_delta(E_dense_short, A1_fit, A2_fit, A3_fit, A4_fit, A5_fit, A6_fit, A7_fit, c0_fit, c1_fit, c2_fit, c3_fit, c4_fit, c5_fit)

if has_sparse:
    E_min_all = min(min(E_vals_delta), min(E_vals_delta_sp))
    E_max_all = max(max(E_vals_delta), max(E_vals_delta_sp))
else:
    E_min_all = min(E_vals_delta)
    E_max_all = max(E_vals_delta)

E_dense_all = np.linspace(E_min_all, E_max_all, 1500)
Delta_fit_all = calc_delta(E_dense_all, A1_fit, A2_fit, A3_fit, A4_fit, A5_fit, A6_fit, A7_fit, c0_fit, c1_fit, c2_fit, c3_fit, c4_fit, c5_fit)
Gamma_fit_all = fit_gamma(E_dense_all, A1_fit, A2_fit, A3_fit, A4_fit, A5_fit, A6_fit, A7_fit)

if aspect_ratio == 'golden_ratio':
    fig, ax1 = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
else:
    fig, ax1 = plt.subplots(figsize=(SIZE, SIZE))

ax1.plot(E_pos, Gamma_pos, 'o', color=color_data, markersize=4, alpha=0.8, zorder=2)
ax1.plot(E_dense_pos, Gamma_fit, color=color_fit, linestyle='-', linewidth=2, zorder=3)

ax1.set_xlabel(X_axis)
ax1.set_ylabel(Y_axis_gamma)

lines = [
    plt.Line2D([], [], color=color_data, marker='o', linestyle='', markersize=4),
    plt.Line2D([], [], color=color_fit, linestyle='-', linewidth=2)
]
labels = [rf'Data: $\Gamma(R_0, E)$', r'Fit: $\Gamma(E)$' ]
ax1.legend(lines, labels, **legend_params)

ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax1.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

pdf_gamma = f"{output_path}/GammaFit_{args.model}.pdf"
fig.savefig(pdf_gamma, format='pdf')
plt.close(fig)
print(f"--> Gamma graph saved to {pdf_gamma}")

if aspect_ratio == 'golden_ratio':
    fig2, ax2 = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
else:
    fig2, ax2 = plt.subplots(figsize=(SIZE, SIZE))

ax2.plot(E_pos, Delta_pos, 'o', color=color_data, markersize=4, alpha=0.8, zorder=2)
ax2.plot(E_dense_pos, Delta_fit_pos, color=color_fit, linestyle='-', linewidth=2, zorder=3)

ax2.set_xlabel(X_axis)
ax2.set_ylabel(Y_axis_delta)

lines2 = [
    plt.Line2D([], [], color=color_data, marker='o', linestyle='', markersize=4),
    plt.Line2D([], [], color=color_fit, linestyle='-', linewidth=2)
]
labels2 = [rf'Data: $\Delta(R_0, E > 0)$', r'Calc: $\Delta(E)$' ]
ax2.legend(lines2, labels2, **legend_params)

ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax2.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

pdf_delta_pos = f"{output_path}/DeltaCalcPos_{args.model}.pdf"
fig2.savefig(pdf_delta_pos, format='pdf')
plt.close(fig2)
print(f"--> Delta (E > 0) graph saved to {pdf_delta_pos}")

if aspect_ratio == 'golden_ratio':
    fig2b, ax2b = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
else:
    fig2b, ax2b = plt.subplots(figsize=(SIZE, SIZE))

if y_lim_delta is not None:
    ax2b.set_ylim(y_lim_delta)

ax2b.plot(E_vals_delta, Delta_R0, 'o', color=color_data, markersize=4, alpha=0.8, zorder=2)
ax2b.plot(E_dense_short, Delta_fit_short, color=color_fit, linestyle='-', linewidth=2, zorder=3)

ax2b.axvline(x=0.0, color='gray', linestyle='--', linewidth=1.0, zorder=1)
ax2b.axhline(y=0.0, color='gray', linestyle='--', linewidth=1.0, zorder=1)

ax2b.set_xlabel(X_axis)
ax2b.set_ylabel(Y_axis_delta)

lines2b = [
    plt.Line2D([], [], color=color_data, marker='o', linestyle='', markersize=4),
    plt.Line2D([], [], color=color_fit, linestyle='-', linewidth=2)
]
labels2b = [rf'Data: $\Delta(R_0, E)$', r'Calc: $\Delta(E)$' ]
ax2b.legend(lines2b, labels2b, **legend_params)

ax2b.xaxis.set_minor_locator(AutoMinorLocator())
ax2b.yaxis.set_minor_locator(AutoMinorLocator())
ax2b.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax2b.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

pdf_delta_short = f"{output_path}/DeltaCalcShort_{args.model}.pdf"
fig2b.savefig(pdf_delta_short, format='pdf')
plt.close(fig2b)
print(f"--> Delta (Short All) graph saved to {pdf_delta_short}")

if aspect_ratio == 'golden_ratio':
    fig3, ax3 = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
else:
    fig3, ax3 = plt.subplots(figsize=(SIZE, SIZE))

if y_lim_delta is not None:
    ax3.set_ylim(y_lim_delta)

if has_sparse:
    ax3.plot(E_vals_delta_sp, Delta_R0_sp, 'o', markeredgecolor='gray', markerfacecolor='none', markersize=5, alpha=0.7, zorder=1)

ax3.plot(E_vals_delta, Delta_R0, 'o', color=color_data, markersize=4, alpha=0.8, zorder=2)
ax3.plot(E_dense_all, Delta_fit_all, color=color_fit, linestyle='-', linewidth=2, zorder=3)

ax3.axvline(x=0.0, color='gray', linestyle='--', linewidth=1.0, zorder=1)
ax3.axhline(y=0.0, color='gray', linestyle='--', linewidth=1.0, zorder=1)

ax3.set_xlabel(X_axis)
ax3.set_ylabel(Y_axis_delta)

lines3 = [
    plt.Line2D([], [], color=color_data, marker='o', linestyle='', markersize=4),
    plt.Line2D([], [], color=color_fit, linestyle='-', linewidth=2)
]

if has_sparse:
    lines3.insert(0, plt.Line2D([], [], color='gray', marker='o', markerfacecolor='none', linestyle='', markersize=5))
    labels3 = ['Sparse Data', rf'Data: $\Delta(R_0, E)$', r'Calc: $\Delta(E)$']
else:
    labels3 = [rf'Data: $\Delta(R_0, E)$', r'Calc: $\Delta(E)$']

ax3.legend(lines3, labels3, **legend_params)

ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())
ax3.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax3.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

pdf_delta_all = f"{output_path}/DeltaCalcAll_{args.model}.pdf"
fig3.savefig(pdf_delta_all, format='pdf')
plt.close(fig3)
print(f"--> Full Delta graph saved to {pdf_delta_all}")

if aspect_ratio == 'golden_ratio':
    fig4, ax4 = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
else:
    fig4, ax4 = plt.subplots(figsize=(SIZE, SIZE))

if has_sparse:
    idx_E_pos_sp = np.where(E_vals_gamma_sp > 0.0)[0]
    E_pos_sp = E_vals_gamma_sp[idx_E_pos_sp]
    Gamma_pos_sp = Gamma_R0_sp[idx_E_pos_sp]
    ax4.plot(E_pos_sp, Gamma_pos_sp, 'o', markeredgecolor='gray', markerfacecolor='none', markersize=5, alpha=0.7, zorder=1)

ax4.plot(E_pos, Gamma_pos, 'o', color=color_data, markersize=4, alpha=0.8, zorder=2)

idx_dense_pos_all = np.where(E_dense_all > 0.0)[0]
ax4.plot(E_dense_all[idx_dense_pos_all], Gamma_fit_all[idx_dense_pos_all], color=color_fit, linestyle='-', linewidth=2, zorder=3)

ax4.set_xlabel(X_axis)
ax4.set_ylabel(Y_axis_gamma)

lines4 = [
    plt.Line2D([], [], color=color_data, marker='o', linestyle='', markersize=4),
    plt.Line2D([], [], color=color_fit, linestyle='-', linewidth=2)
]

if has_sparse:
    lines4.insert(0, plt.Line2D([], [], color='gray', marker='o', markerfacecolor='none', linestyle='', markersize=5))
    labels4 = ['Sparse Data', rf'Data: $\Gamma(R_0, E)$', r'Fit: $\Gamma(E)$']
else:
    labels4 = [rf'Data: $\Gamma(R_0, E)$', r'Fit: $\Gamma(E)$']

ax4.legend(lines4, labels4, **legend_params)

ax4.xaxis.set_minor_locator(AutoMinorLocator())
ax4.yaxis.set_minor_locator(AutoMinorLocator())
ax4.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax4.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

pdf_gamma_all = f"{output_path}/GammaFitAll_{args.model}.pdf"
fig4.savefig(pdf_gamma_all, format='pdf')
plt.close(fig4)
print(f"--> Full Gamma graph saved to {pdf_gamma_all}")