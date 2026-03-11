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

parser = argparse.ArgumentParser(description="Script for fitting background defect vs Energy.")
parser.add_argument('--model', type=str, default='2DModel2', help='Select the calculation mode defined in settings_fit.py')
args = parser.parse_args()

if args.model not in settings_fit.MODELS:
    print(f"\n[ERROR] Mode '{args.model}' is not defined in settings_fit.py!")
    sys.exit(1)

cfg = settings_fit.MODELS[args.model]
print(f"--> Running {os.path.basename(__file__)} in mode: {args.model}")

R0_target = cfg.get('Fit_R0', 2.0)

base_path = cfg.get('defect_path', r"./DATA/RydbergStates")
output_path = cfg.get('output_path', r"./Graphs/Fitting")
os.makedirs(output_path, exist_ok=True)
aspect_ratio = cfg.get('aspect_ratio', 'square')
SIZE = 6
SMALL_SIZE = 14; MEDIUM_SIZE = 18; BIGGER_SIZE = 24; LEGEND_SIZE = 12
X_axis = r'Electron energy $\epsilon\,(\mathrm{eV})$'
Y_axis = r'Background quantum defect phase $\pi\mu_\mathrm{bg}(R_0,\epsilon)$'
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
    data_ds = np.loadtxt(r"./DATA/DSState/DSenergies.txt", comments="#")
    R_vals = data_ds[:, 0]
except OSError:
    print("[ERROR] Could not load DSenergies.txt to find R values.")
    sys.exit(1)

idx_R0 = np.argmin(np.abs(R_vals - R0_target))
R0_actual = R_vals[idx_R0]

try:
    defect_data = np.loadtxt(f"{base_path}/defects_PHP.txt", comments="#")
    n_vals = defect_data[:, 0]
    mu_matrix = defect_data[:, 1:]
except OSError:
    print("[ERROR] Could not load defects_PHP.txt.")
    sys.exit(1)

try:
    E_data = np.loadtxt(f"{base_path}/eigE_PHP.txt", comments="#")
    E_matrix = E_data[:, 1:]
except OSError:
    print("[ERROR] Could not load eigE_PHP.txt.")
    sys.exit(1)

mu_R0 = mu_matrix[:, idx_R0]
E_R0 = E_matrix[:, idx_R0]

mu_std_other = np.zeros(len(n_vals))
E_std_other = np.zeros(len(n_vals))

for i in range(len(n_vals)):
    mu_row_data = mu_matrix[i, :]
    mu_other_R = np.delete(mu_row_data, idx_R0)
    mu_std_other[i] = np.std(mu_other_R)
    
    E_row_data = E_matrix[i, :]
    E_other_R = np.delete(E_row_data, idx_R0)
    E_std_other[i] = np.std(E_other_R)

mu_std_other[mu_std_other == 0] = 1e-8
E_std_other[E_std_other == 0] = 1e-8

def fit_func_quadratic_E(E, a, b, c):
    return c / (1 + a * E**2 + b * E)

p0_defect = cfg.get('p0_guess_defect', [0.0, 0.0, -2.0])
if len(p0_defect) == 2:
    p0_defect = [0.0] + p0_defect

try:
    popt, pcov = curve_fit(fit_func_quadratic_E, E_R0, mu_R0, p0=p0_defect, sigma=mu_std_other, absolute_sigma=True)
    print(f"--> Defect Fit converged. Parameters [a, b, c]: {popt}")
except RuntimeError:
    print("[ERROR] Optimal parameters not found for defect fit.")
    sys.exit(1)

E_dense = np.linspace(min(E_R0), max(E_R0), 500)
mu_fit = fit_func_quadratic_E(E_dense, *popt)

if aspect_ratio == 'golden_ratio':
    fig, ax1 = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
else:
    fig, ax1 = plt.subplots(figsize=(SIZE, SIZE))

ax1.errorbar(E_R0, mu_R0, xerr=E_std_other, yerr=mu_std_other, fmt='o', color=color_data, 
             ecolor='gray', elinewidth=1, capsize=3, markersize=4, alpha=0.8, zorder=2)

ax1.plot(E_dense, mu_fit, color=color_fit, linestyle='-', linewidth=2, zorder=3)

ax1.set_xlabel(X_axis)
ax1.set_ylabel(Y_axis)

lines = [
    plt.Line2D([], [], color=color_data, marker='o', linestyle='', markersize=4),
    plt.Line2D([], [], color=color_fit, linestyle='-', linewidth=2)
]
labels = [
    r'Data: $\mu_\mathrm{bg}(R_0) \pm \Delta\mu_\mathrm{bg}(R)$',
    r'Fit'
]

ax1.legend(lines, labels, **legend_params)

ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax1.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

pdf_filename = f"{output_path}/DefectFit_{args.model}.pdf"
pgf_filename = f"{output_path}/DefectFit_{args.model}.pgf"
fig.savefig(pdf_filename, format='pdf')
fig.savefig(pgf_filename, format='pgf')
plt.close(fig)
print(f"--> Graph saved to {pdf_filename}")