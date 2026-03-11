import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from scipy.optimize import curve_fit
import gc
import os
import sys
import argparse

script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.append(script_dir)
import settings_fit

parser = argparse.ArgumentParser(description="Script for fitting Gamma separability.")
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
X_axis = cfg.get('X_axis_g', r'Internuclear distance $R (a_0)$')
Y_axis = cfg.get('Y_axis_g', r'Scaling factor $g^2(R)$')
color_data = cfg.get('color_data', 'black')
color_fit = cfg.get('color_fit', 'red')

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

try:
    data_ds = np.loadtxt(f"{base_path}/DSenergies.txt", comments="#")
    R_vals = data_ds[:, 0]
    data_gamma = np.loadtxt(f"{base_path}/gamma.txt", comments="#")
    E_vals = data_gamma[:, 0]
    Gamma_matrix = data_gamma[:, 1:]
except OSError:
    print("[ERROR] Could not load data files.")
    sys.exit(1)

idx_R0 = np.argmin(np.abs(R_vals - R0_target))
R0_actual = R_vals[idx_R0]

idx_E = np.where(E_vals > 0.0)[0]
E_min_actual, E_max_actual = E_vals[idx_E[0]], E_vals[idx_E[-1]]

Gamma_ref = Gamma_matrix[idx_E, idx_R0]
S_R = np.zeros(len(R_vals))
S_err = np.zeros(len(R_vals))

for i in range(len(R_vals)):
    ratio = Gamma_matrix[idx_E, i] / Gamma_ref
    S_R[i] = np.mean(ratio)
    S_err[i] = np.std(ratio) 

S_err[S_err == 0] = 1e-8



# General fitting function adapting dynamically to the number of parameters
def custom_shape(R, *params):
    alpha = params[0]
    Rc = params[1]
    base = 1.0 - np.tanh(alpha * (R - Rc))
    
    A_g1 = params[2]
    beta1 = params[3]
    Rw1 = params[4]
    bump1 = A_g1 * np.exp(-beta1 * (R - Rw1)**2)
    
    bump2 = 0.0
    if len(params) >= 8:
        A_g2 = params[5]
        beta2 = params[6]
        Rw2 = params[7]
        bump2 = A_g2 * np.exp(-beta2 * (R - Rw2)**2)
        
    return base + bump1 + bump2

# Normalized fitting function
def fit_func_generic(R, *params):
    val_at_R0 = custom_shape(R0_actual, *params)
    return custom_shape(R, *params) / val_at_R0

p0 = cfg.get('p0_guess_g', [1.0, 2.0, 1.0, 5.0, 0.5])

try:
    popt, pcov = curve_fit(fit_func_generic, R_vals, S_R, p0=p0, sigma=S_err, absolute_sigma=True, maxfev=10000)
    print(f"--> Fit converged with {len(popt)} parameters. Parameters: {popt}")
except RuntimeError:
    print("[ERROR] Optimal parameters not found for fit.")
    sys.exit(1)

R_dense = np.linspace(min(R_vals), max(R_vals), 500)
S_fit = fit_func_generic(R_dense, *popt)

if aspect_ratio == 'golden_ratio':
    fig, ax1 = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
else:
    fig, ax1 = plt.subplots(figsize=(SIZE, SIZE))

ax1.errorbar(R_vals, S_R, yerr=S_err, fmt='o', color=color_data, 
             ecolor='gray', elinewidth=1, capsize=3, markersize=4, alpha=0.8, zorder=2)

ax1.plot(R_dense, S_fit, color=color_fit, linestyle='-', linewidth=2, zorder=3)
ax1.plot(R0_actual, 1.0, marker='D', color='blue', markersize=6, zorder=4)

ax1.set_xlabel(X_axis)
ax1.set_ylabel(Y_axis)

lines = [
    plt.Line2D([], [], color=color_data, marker='o', linestyle='', markersize=4),
    plt.Line2D([], [], color=color_fit, linestyle='-', linewidth=2),
    plt.Line2D([], [], color='blue', marker='D', linestyle='', markersize=6)
]
labels = [
    r'Data: $\langle \Gamma(R) / \Gamma(R_0) \rangle_\epsilon$',
    r'Fit',
    rf'Anchor $(R_0, 1.0)$'
]

ax1.legend(lines, labels, **legend_params)

ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax1.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

pdf_filename = f"{output_path}/SeparabilityFit_{args.model}.pdf"
pgf_filename = f"{output_path}/SeparabilityFit_{args.model}.pgf"
fig.savefig(pdf_filename, format='pdf')
fig.savefig(pgf_filename, format='pgf')
plt.close(fig)
print(f"--> Graph saved to {pdf_filename}")


if len(popt) == 5:
    param_names = ['g_alpha', 'g_Rc', 'g_A_g', 'g_beta', 'g_Rw']
elif len(popt) == 8:
    param_names = ['g_alpha', 'g_Rc', 'g_A_g1', 'g_beta1', 'g_Rw1', 'g_A_g2', 'g_beta2', 'g_Rw2']
else:
    param_names = [f'g_param_{i}' for i in range(len(popt))]

fortran_lines = []
for name, val in zip(param_names, popt):
    fortran_val = f"{val:.15e}".replace('e', 'd')
    fortran_lines.append(f"REAL(KIND = idk), PARAMETER :: {name} = {fortran_val}")

val_at_R0 = custom_shape(R0_actual, *popt)
g_norm = 1.0 / val_at_R0

fortran_g_norm = f"{g_norm:.15e}".replace('e', 'd')
fortran_lines.append(f"REAL(KIND = idk), PARAMETER :: g_norm = {fortran_g_norm}")

print("\n--> Fortran parameters ready to copy:")
for line in fortran_lines:
    print(line)

txt_filename = f"{output_path}/SeparabilityFit_{args.model}.txt"
with open(txt_filename, 'w') as f:
    for line in fortran_lines:
        f.write(line + "\n")

print(f"--> Fortran parameters saved to {txt_filename}")