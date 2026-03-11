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

parser = argparse.ArgumentParser(description="Script for fitting DS energies vs R.")
parser.add_argument('--model', type=str, default='2DModel2', help='Select the calculation mode defined in settings_fit.py')
args = parser.parse_args()

if args.model not in settings_fit.MODELS:
    print(f"\n[ERROR] Mode '{args.model}' is not defined in settings_fit.py!")
    sys.exit(1)

cfg = settings_fit.MODELS[args.model]
print(f"--> Running {os.path.basename(__file__)} in mode: {args.model}")

base_path = cfg.get('base_path', r"./DATA/DSState")
output_path = cfg.get('output_path', r"./Graphs/Fitting")
os.makedirs(output_path, exist_ok=True)
aspect_ratio = cfg.get('aspect_ratio', 'square')
SIZE = 6
SMALL_SIZE = 14; MEDIUM_SIZE = 18; BIGGER_SIZE = 24; LEGEND_SIZE = 12

X_axis = cfg.get('X_axis_R', r'Internuclear distance $R\,(a_0)$')
Y_axis = cfg.get('Y_axis_E', r'Discrete state energy $\epsilon_d\,(\mathrm{eV})$')
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
    E_vals = data_ds[:, 1]
except OSError:
    print("[ERROR] Could not load DSenergies.txt.")
    sys.exit(1)

p0_E = cfg.get('p0_guess_ds_energy', [-15.0, -2.0, 1.0, 2.0, 1.0, 0.5, 2.0])

if len(p0_E) == 4:
    def fit_func_E(R, const, A, alpha, Rc):
        return const + A * np.tanh(alpha * (R - Rc))
else:
    def fit_func_E(R, const, A, alpha, Rc, B, beta, Rw):
        return const + A * np.tanh(alpha * (R - Rc)) + B * np.exp(-beta * (R - Rw)**2)

try:
    popt, pcov = curve_fit(fit_func_E, R_vals, E_vals, p0=p0_E, maxfev=10000)
    if len(p0_E) == 4:
        print(f"--> DS Energy Fit converged. Parameters [const, A, alpha, Rc]: {popt}")
    else:
        print(f"--> DS Energy Fit converged. Parameters [const, A, alpha, Rc, B, beta, Rw]: {popt}")
except RuntimeError:
    print("[ERROR] Optimal parameters not found for DS Energy fit.")
    sys.exit(1)

R_dense = np.linspace(min(R_vals), max(R_vals), 1000)
E_fit = fit_func_E(R_dense, *popt)

if aspect_ratio == 'golden_ratio':
    fig, ax1 = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
else:
    fig, ax1 = plt.subplots(figsize=(SIZE, SIZE))

ax1.plot(R_vals, E_vals, 'o', color=color_data, markersize=4, alpha=0.8, zorder=2)
ax1.plot(R_dense, E_fit, color=color_fit, linestyle='-', linewidth=2, zorder=3)

ax1.set_xlabel(X_axis)
ax1.set_ylabel(Y_axis)

lines = [
    plt.Line2D([], [], color=color_data, marker='o', linestyle='', markersize=4),
    plt.Line2D([], [], color=color_fit, linestyle='-', linewidth=2)
]
labels = [
    r'Data: $\epsilon_d(R)$',
    r'Fit'
]

ax1.legend(lines, labels, **legend_params)

ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax1.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

pdf_filename = f"{output_path}/DSEnergyFit_{args.model}.pdf"
pgf_filename = f"{output_path}/DSEnergyFit_{args.model}.pgf"
fig.savefig(pdf_filename, format='pdf')
fig.savefig(pgf_filename, format='pgf')
plt.close(fig)
print(f"--> Graph saved to {pdf_filename}")

if len(p0_E) == 4:
    param_names = ['DSE_const', 'DSE_A', 'DSE_alpha', 'DSE_Rc']
else:
    param_names = ['DSE_const', 'DSE_A', 'DSE_alpha', 'DSE_Rc', 'DSE_B', 'DSE_beta', 'DSE_Rw']

phys_h0 = 27.2113961

fortran_lines = []
for i, (name, val) in enumerate(zip(param_names, popt)):
    if i in [0, 1, 4]:
        val_converted = val / phys_h0
    else:
        val_converted = val
        
    fortran_val = f"{val_converted:.15e}".replace('e', 'd')
    fortran_lines.append(f"REAL(KIND = idk), PARAMETER :: {name} = {fortran_val}")

print("\n--> Fortran parameters ready to copy (Energies converted to Hartree):")
for line in fortran_lines:
    print(line)

txt_filename = f"{output_path}/DSEnergyFit_{args.model}.txt"
with open(txt_filename, 'w') as f:
    for line in fortran_lines:
        f.write(line + "\n")

print(f"--> Fortran parameters saved to {txt_filename}")