import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
import scipy.special as sp 
import os
import sys
import argparse

script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.append(script_dir)
import settings_fit

parser = argparse.ArgumentParser(description="Script for isolating resonance from background.")
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
Y_axis_gamma = r'Isolated $\Gamma_{res}(R_0, E)\, (\mathrm{eV})$'
Y_axis_delta = r'Isolated $\Delta_{res}(R_0, E)\, (\mathrm{eV})$'

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

# --- NAČTENÍ HUSTÝCH DAT (PRO ZÍSKÁNÍ A0_fixed, c0_fixed) ---
data_ds = np.loadtxt(f"{base_path}/DSenergies.txt", comments="#")
R_vals = data_ds[:, 0]
idx_R0 = np.argmin(np.abs(R_vals - R0_target))

data_gamma = np.loadtxt(f"{base_path}/gamma.txt", comments="#")
E_vals_gamma = data_gamma[:, 0]
Gamma_R0 = data_gamma[:, idx_R0 + 1]

data_delta = np.loadtxt(f"{base_path}/delta.txt", comments="#")
E_vals_delta = data_delta[:, 0]
Delta_R0 = data_delta[:, idx_R0 + 1]

idx_E_pos = np.where(E_vals_gamma > 0.0)[0]
Gamma_pos = Gamma_R0[idx_E_pos]
Delta_pos = Delta_R0[idx_E_pos]

A0_fixed = Gamma_pos[0]
c0_fixed = Delta_pos[0] * (2.0 * np.pi / A0_fixed)

base_path_sparse = cfg.get('base_path_sparse', r"./DATA/DSStateSparce")
try:
    data_gamma_sp = np.loadtxt(f"{base_path_sparse}/gamma.txt", comments="#")
    E_vals_gamma_sp = data_gamma_sp[:, 0]
    Gamma_R0_sp = data_gamma_sp[:, idx_R0 + 1]

    data_delta_sp = np.loadtxt(f"{base_path_sparse}/delta.txt", comments="#")
    E_vals_delta_sp = data_delta_sp[:, 0]
    Delta_R0_sp = data_delta_sp[:, idx_R0 + 1]
except OSError:
    print(f"[ERROR] Sparse data not found at {base_path_sparse}. Cannot perform background subtraction.")
    sys.exit(1)

def calc_defect(E_eV):
    return -2.17015937 / (1 -0.03881984 * E_eV + 0.00469617 * E_eV**2)

def F_threshold(E_eV):
    E_au = E_eV / 27.211386  
    k = np.sqrt(2 * E_au + 0j)
    
    mu = calc_defect(E_eV)
    mu_eff = np.where(E_eV < 0, mu, 0.0)
    
    term1 = c0_fixed
    term2 = -0.5 * 1j * k
    term3 = np.log(-1j * k)
    term4 = sp.digamma(1.0 - 1j / k - mu_eff)
                                  
    F_val = (A0_fixed / (2.0 * np.pi)) * (term1 + term2 + term3 + term4)
    return F_val

Gamma_thresh = -2.0 * np.imag(F_threshold(E_vals_gamma_sp))
Gamma_isolated = Gamma_R0_sp - Gamma_thresh

Delta_thresh = np.real(F_threshold(E_vals_delta_sp))
Delta_isolated = Delta_R0_sp - Delta_thresh


if aspect_ratio == 'golden_ratio':
    fig1, ax1 = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
    fig2, ax2 = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
else:
    fig1, ax1 = plt.subplots(figsize=(SIZE, SIZE))
    fig2, ax2 = plt.subplots(figsize=(SIZE, SIZE))

ax1.plot(E_vals_gamma_sp, Gamma_isolated, 'o', color='blue', markersize=5, alpha=0.8, label='Isolated Sparse Data')
ax1.axhline(0, color='gray', linestyle='--', linewidth=1.0)
ax1.set_xlabel(X_axis)
ax1.set_ylabel(Y_axis_gamma)
ax1.legend(**legend_params)
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='both', direction='in', which='both', top=True, right=True)

pdf_gamma_iso = f"{output_path}/Gamma_Isolated_Sparse_{args.model}.pdf"
fig1.savefig(pdf_gamma_iso, format='pdf')
plt.close(fig1)
print(f"--> Isolated Gamma graph saved to {pdf_gamma_iso}")

ax2.plot(E_vals_delta_sp, Delta_isolated, 'o', color='green', markersize=5, alpha=0.8, label='Isolated Sparse Data')
ax2.axhline(0, color='gray', linestyle='--', linewidth=1.0)
ax2.axvline(0, color='gray', linestyle='--', linewidth=1.0)
ax2.set_xlabel(X_axis)
ax2.set_ylabel(Y_axis_delta)
if cfg.get('y_lim_delta', None) is not None:
    ax2.set_ylim(cfg.get('y_lim_delta'))
ax2.legend(**legend_params)
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.tick_params(axis='both', direction='in', which='both', top=True, right=True)

pdf_delta_iso = f"{output_path}/Delta_Isolated_Sparse_{args.model}.pdf"
fig2.savefig(pdf_delta_iso, format='pdf')
plt.close(fig2)
print(f"--> Isolated Delta graph saved to {pdf_delta_iso}")