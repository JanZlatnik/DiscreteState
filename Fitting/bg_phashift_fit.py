import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from scipy.optimize import curve_fit
import os
import sys
import argparse

script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.append(script_dir)
import settings_fit

parser = argparse.ArgumentParser(description="Script for fitting 2D background phase shift.")
parser.add_argument('--model', type=str, default='2DModel2')
args = parser.parse_args()

if args.model not in settings_fit.MODELS:
    sys.exit(1)

cfg = settings_fit.MODELS[args.model]

base_path = cfg.get('base_path', r"./DATA/DSState")
output_path = cfg.get('output_path', r"./Graphs/Fitting")
os.makedirs(output_path, exist_ok=True)
aspect_ratio = cfg.get('aspect_ratio', 'square')
SIZE = 6
SMALL_SIZE = 14; MEDIUM_SIZE = 18; BIGGER_SIZE = 24; LEGEND_SIZE = 12

X_axis_E = r'Electron energy $\epsilon\,(\mathrm{eV})$'
X_axis_R = r'Internuclear distance $R\,(a_0)$'
Y_axis = r'Background phase shift $\delta_\mathrm{bg}(R,\epsilon)$'

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
legend_params_R = {'frameon': False, 'ncol': 1, 'fontsize': 12, 'handlelength': 2, 'loc': 'best', 'bbox_to_anchor':(1, 0.90)}

try:
    data_ds = np.loadtxt(f"{base_path}/DSenergies.txt", comments="#")
    R_vals = data_ds[:, 0]
    
    data_tot = np.loadtxt(f"{base_path}/phaseshift.txt", comments="#")
    E_vals = data_tot[:, 0]
    phase_tot = data_tot[:, 1:]
    
    data_res = np.loadtxt(f"{base_path}/DS_phaseshift.txt", comments="#")
    phase_res = data_res[:, 1:]
except OSError:
    sys.exit(1)

phase_bg = phase_tot - phase_res

for i in range(phase_bg.shape[1]):
    phase_bg[:, i] = np.unwrap(phase_bg[:, i] * 2.0) / 2.0

for i in range(1, phase_bg.shape[1]):
    diff = phase_bg[0, i] - phase_bg[0, i-1]
    n_pi = np.round(diff / np.pi)
    phase_bg[:, i] -= n_pi * np.pi

global_shift = np.round(phase_bg[0, 0] / np.pi) * np.pi
phase_bg -= global_shift

E_grid, R_grid = np.meshgrid(E_vals, R_vals, indexing='ij')
E_flat = E_grid.flatten()
R_flat = R_grid.flatten()
Z_flat = phase_bg.flatten()

def fit_func_bg(X, *params):
    E, R = X
    cE0 = params[0]
    cE1 = params[1]
    cE2 = params[2]
    cArctan = params[3]
    cEArctan = params[4]
    cE2Arctan = params[5]
    cW = params[6]
    Rw = params[7]
    return cE0 + cE1*E + cE2*(E**2) + (cArctan + cEArctan*E + cE2Arctan*(E**2))*np.arctan(cW * (R - Rw))

p0_bg = cfg.get('p0_guess_bg_phase', [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.5])

try:
    popt, pcov = curve_fit(fit_func_bg, (E_flat, R_flat), Z_flat, p0=p0_bg, maxfev=20000)
    
    Z_fit = fit_func_bg((E_flat, R_flat), *popt)
    residuals = Z_flat - Z_fit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((Z_flat - np.mean(Z_flat))**2)
    r_squared = 1.0 - (ss_res / ss_tot)
    rmse = np.sqrt(np.mean(residuals**2))
    
    print(f"--> BG Phase Fit converged with {len(popt)} parameters.")
    print(f"    R-squared: {r_squared:.6f}")
    print(f"    RMSE:      {rmse:.6f} rad")
except RuntimeError:
    sys.exit(1)

if aspect_ratio == 'golden_ratio':
    fig1, ax1 = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
    fig2, ax2 = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
else:
    fig1, ax1 = plt.subplots(figsize=(SIZE, SIZE))
    fig2, ax2 = plt.subplots(figsize=(SIZE, SIZE))

target_Rs = [0.7, 2.0, 10.0]
colors_R = ['#1f77b4', '#2ca02c', '#d62728']
custom_lines_R = [plt.Line2D([0], [0], color='black', marker='o', linestyle='', markersize=4, alpha=1.0)]
custom_labels_R = [r'Data: $\delta_\mathrm{bg}(R,\epsilon)$']

for c, tr in zip(colors_R, target_Rs):
    idx = np.argmin(np.abs(R_vals - tr))
    R_current = R_vals[idx]
    ax1.plot(E_vals, phase_bg[:, idx], 'o', color='black', markersize=4, alpha=1.0, markeredgecolor='none', zorder=2)
    fit_curve = fit_func_bg((E_vals, R_current), *popt)
    ax1.plot(E_vals, fit_curve, '-', color=c, linewidth=2.5, zorder=3)
    custom_lines_R.append(plt.Line2D([0], [0], color=c, lw=2.5))
    custom_labels_R.append(f'$R = {R_current:.2f}\ a_0$')

ax1.set_xlabel(X_axis_E)
ax1.set_ylabel(Y_axis)
ax1.legend(custom_lines_R, custom_labels_R, **legend_params)
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax1.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

fig1.savefig(f"{output_path}/BG_PhaseShift_Fit_E_{args.model}.pdf", format='pdf')
fig1.savefig(f"{output_path}/BG_PhaseShift_Fit_E_{args.model}.pgf", format='pgf')
plt.close(fig1)

target_Es = [0.0, 1.0, 2.0]
colors_E = ['#1f77b4', '#2ca02c', '#d62728']
custom_lines_E = [plt.Line2D([0], [0], color='black', marker='o', linestyle='', markersize=4, alpha=1.0)]
custom_labels_E = [r'Data: $\delta_\mathrm{bg}(R,\epsilon)$']

for c, te in zip(colors_E, target_Es):
    idx = np.argmin(np.abs(E_vals - te))
    E_current = E_vals[idx]
    ax2.plot(R_vals, phase_bg[idx, :], 'o', color='black', markersize=4, alpha=1.0, markeredgecolor='none', zorder=2)
    fit_curve = fit_func_bg((E_current, R_vals), *popt)
    ax2.plot(R_vals, fit_curve, '-', color=c, linewidth=2.5, zorder=3)
    custom_lines_E.append(plt.Line2D([0], [0], color=c, lw=2.5))
    custom_labels_E.append(rf'$\epsilon = {E_current:.2f}\ \mathrm{{eV}}$')

ax2.set_xlabel(X_axis_R)
ax2.set_ylabel(Y_axis)
ax2.legend(custom_lines_E, custom_labels_E, **legend_params_R)
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax2.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

fig2.savefig(f"{output_path}/BG_PhaseShift_Fit_R_{args.model}.pdf", format='pdf')
fig2.savefig(f"{output_path}/BG_PhaseShift_Fit_R_{args.model}.pgf", format='pgf')
plt.close(fig2)

param_names = ['bg_cE0', 'bg_cE1', 'bg_cE2', 'bg_cArctan', 'bg_cEArctan', 'bg_cE2Arctan', 'bg_cW', 'bg_Rw']

fortran_lines = []
for name, val in zip(param_names, popt):
    fortran_val = f"{val:.15e}".replace('e', 'd')
    fortran_lines.append(f"REAL(KIND = idk), PARAMETER :: {name} = {fortran_val}")

print("\n--> Fortran parameters ready to copy:")
for line in fortran_lines:
    print(line)

txt_filename = f"{output_path}/BG_PhaseShift_Fit_{args.model}.txt"
with open(txt_filename, 'w') as f:
    for line in fortran_lines:
        f.write(line + "\n")

print(f"--> Fortran parameters saved to {txt_filename}")