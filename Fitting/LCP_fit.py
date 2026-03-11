import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import curve_fit
import os
import sys
import argparse

script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.append(script_dir)
import settings_fit

parser = argparse.ArgumentParser(description="Script for extracting and fitting LCP data.")
parser.add_argument('--model', type=str, default='2DModel2')
args = parser.parse_args()

if args.model not in settings_fit.MODELS:
    sys.exit(1)

cfg = settings_fit.MODELS[args.model]
base_path = cfg.get('LCP_base_path', r"./DATA/DSState")
sparse_path = cfg.get('LCP_sparse_path', r"./DATA/DSStateSparce")
output_path = cfg.get('output_path', r"./Graphs/Fitting")
os.makedirs(output_path, exist_ok=True)

aspect_ratio = cfg.get('aspect_ratio', 'square')
SIZE = 6; SMALL_SIZE = 14; MEDIUM_SIZE = 18; LEGEND_SIZE = 12
color_data = cfg.get('color_data_LCP', 'black')
color_fit = cfg.get('color_fit_LCP', 'red')

plt.rcParams['font.family'] = 'serif'
plt.rcParams["font.serif"] = ["Latin Modern Roman"]
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['axes.titlepad'] = 10; plt.rcParams['axes.labelpad'] = 10 
plt.rcParams['legend.fancybox'] = False; plt.rcParams['legend.edgecolor'] = "#000000"
plt.rcParams["figure.autolayout"] = True
plt.rc('font', size=SMALL_SIZE); plt.rc('axes', titlesize=MEDIUM_SIZE) 
plt.rc('axes', labelsize=MEDIUM_SIZE); plt.rc('xtick', labelsize=SMALL_SIZE) 
plt.rc('ytick', labelsize=SMALL_SIZE); plt.rc('legend', fontsize=LEGEND_SIZE) 

def get_combined_data(b_dir, s_dir, filenames):
    E_list, V_list = [], []
    for directory in [b_dir, s_dir]:
        for fname in filenames:
            fpath = os.path.join(directory, fname)
            if os.path.exists(fpath):
                try:
                    data = np.loadtxt(fpath, comments="#")
                    E_list.append(data[:, 0])
                    V_list.append(data[:, 1:])
                except Exception:
                    pass
    if not E_list:
        return None, None
    expected_cols = V_list[0].shape[1]
    valid_E, valid_V = [], []
    for i, (e, v) in enumerate(zip(E_list, V_list)):
        if v.shape[1] == expected_cols:
            valid_E.append(e)
            valid_V.append(v)
    E_all = np.concatenate(valid_E)
    V_all = np.vstack(valid_V)
    unique_E, unique_idx = np.unique(E_all, return_index=True)
    return unique_E, V_all[unique_idx, :]

E_grid, delta_grid = get_combined_data(base_path, sparse_path, ["deltaA.txt", "delta.txt"])
Eg_grid, gamma_grid = get_combined_data(base_path, sparse_path, ["gamma.txt"])

try:
    ds_data = np.loadtxt(f"{base_path}/DSenergies.txt", comments="#")
    R_vals, Ed_vals = ds_data[:, 0], ds_data[:, 1]
except OSError:
    sys.exit(1)

if E_grid is None:
    sys.exit(1)

R_valid, E_res_vals, Delta_res_vals, Gamma_res_vals = [], [], [], []

for i, R in enumerate(R_vals):
    if i >= delta_grid.shape[1]:
        continue
    Ed = Ed_vals[i]
    delta_col = delta_grid[:, i]
    g_func = E_grid - Ed - delta_col
    crossings = np.where((g_func[:-1] <= 0) & (g_func[1:] >= 0))[0]
    if len(crossings) > 0:
        idx = min(crossings, key=lambda x: abs(E_grid[x] - Ed))
        E1, E2 = E_grid[idx], E_grid[idx+1]
        g1, g2 = g_func[idx], g_func[idx+1]
        if E2 == E1: continue 
        E_res = E1 - g1 * (E2 - E1) / (g2 - g1)
        d1, d2 = delta_col[idx], delta_col[idx+1]
        Delta_res = d1 + (E_res - E1) * (d2 - d1) / (E2 - E1)
        if E_res <= 0.0 or Eg_grid is None or i >= gamma_grid.shape[1]:
            Gamma_res = 0.0
        else:
            gamma_col = gamma_grid[:, i]
            g_idx = max(0, min(np.searchsorted(Eg_grid, E_res) - 1, len(Eg_grid)-2))
            Eg1, Eg2 = Eg_grid[g_idx], Eg_grid[g_idx+1]
            ga1, ga2 = gamma_col[g_idx], gamma_col[g_idx+1]
            if Eg2 == Eg1: continue 
            Gamma_res = ga1 + (E_res - Eg1) * (ga2 - ga1) / (Eg2 - Eg1)
        R_valid.append(R)
        E_res_vals.append(E_res)
        Delta_res_vals.append(Delta_res)
        Gamma_res_vals.append(Gamma_res)

R_valid = np.array(R_valid)
E_res_vals = np.array(E_res_vals)
Delta_res_vals = np.array(Delta_res_vals)
Gamma_res_vals = np.array(Gamma_res_vals)

def fit_E(R, const, A, alpha, Rc, B, beta):
    return const + A * np.tanh(alpha * (R - Rc)) + B * np.exp(-beta * R)

def fit_Delta(R, A, alpha, Rc):
    return A * (np.tanh(alpha * (R - Rc)) - 1)

def fit_Gamma(X, A, B, beta, gamma_exp):
    R, E = X[0], X[1]
    val = np.zeros_like(R)
    mask = E > 0
    val[mask] = (A + B * np.exp(-beta * R[mask])) * (E[mask]**gamma_exp)
    return val

targets = [
    ("E_res", E_res_vals, fit_E, cfg.get('p0_guess_E_res', [0.5, -4.0, 0.5, 3.0, 20.0, 1.0]), r'Resonant energy $\epsilon_\mathrm{res}\,(\mathrm{eV})$', ['const', 'A', 'alpha', 'Rc', 'B', 'beta'], r'Data: $\epsilon_\mathrm{res}$'),
    ("Delta", Delta_res_vals, fit_Delta, cfg.get('p0_guess_Delta', [1.3, 0.8, 1.5]), r'Level shift $\Delta_\mathrm{res}\,(\mathrm{eV})$', ['A', 'alpha', 'Rc'], r'Data: $\Delta_\mathrm{res}$'),
    ("Gamma", Gamma_res_vals, fit_Gamma, cfg.get('p0_guess_Gamma', [0.1, 0.5, 0.5, 1.5]), r'Resonance width $\Gamma_\mathrm{res}\,(\mathrm{eV})$', ['A', 'B', 'beta', 'gamma_exp'], r'Data: $\Gamma_\mathrm{res}$')
]

phys_h0 = 27.2113961
all_fortran_lines = []

for name, y_vals, fit_func, p0, ylabel, param_names, data_label in targets:
    if len(R_valid) == 0: continue
    R_dense = np.linspace(0.7, 10.5, 1000)
    popt = None
    y_fit = None
    try:
        if name == "Gamma":
            x_data = np.vstack((R_valid, E_res_vals))
            popt, _ = curve_fit(fit_func, x_data, y_vals, p0=p0, maxfev=100000)
            E_dense = np.interp(R_dense, R_valid, E_res_vals)
            y_fit = fit_func(np.vstack((R_dense, E_dense)), *popt)
        else:
            popt, _ = curve_fit(fit_func, R_valid, y_vals, p0=p0, maxfev=100000)
            y_fit = fit_func(R_dense, *popt)
    except RuntimeError:
        pass

    if aspect_ratio == 'golden_ratio':
        fig, ax = plt.subplots(figsize=(SIZE * ((1 + 5 ** 0.5) / 2), SIZE))
    else:
        fig, ax = plt.subplots(figsize=(SIZE, SIZE))

    ax.plot(R_valid, y_vals, 'o', color=color_data, markersize=4, alpha=0.8, zorder=2)
    lines = [plt.Line2D([],[], color=color_data, marker='o', ls='')]
    labels = [data_label]
    if y_fit is not None:
        ax.plot(R_dense, y_fit, color=color_fit, linestyle='-', linewidth=2, zorder=3)
        lines.append(plt.Line2D([],[], color=color_fit, ls='-'))
        labels.append('Fit')
    ax.set_xlabel(r'Internuclear distance $R\,(a_0)$')
    ax.set_ylabel(ylabel)
    ax.set_xlim(0.2, 10.5)
    ax.legend(lines, labels, frameon=False)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', direction='in', which='both', top=True, right=True)
    fig.savefig(f"{output_path}/LCP_{name}_{args.model}.pdf", format='pdf')
    plt.close(fig)

    if popt is not None:
        all_fortran_lines.append(f"! --- {name} ---")
        for i, (pname, val) in enumerate(zip(param_names, popt)):
            val_converted = val
            if name in ["E_res", "Delta"] and pname in ["const", "A", "B"]:
                val_converted = val / phys_h0
            elif name == "Gamma" and (pname == "A" or pname == "B"):
                val_converted = val / phys_h0
            fortran_val = f"{val_converted:.15e}".replace('e', 'd')
            all_fortran_lines.append(f"REAL(KIND = idk), PARAMETER :: LCP_{name}_{pname} = {fortran_val}")
        all_fortran_lines.append("")

print("\n".join(all_fortran_lines))
txt_filename = f"{output_path}/LCP_Parameters_{args.model}.txt"
with open(txt_filename, 'w') as f:
    f.write("\n".join(all_fortran_lines))