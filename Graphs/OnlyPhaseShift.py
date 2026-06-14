import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator
import gc
import os

import sys
import argparse
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.append(script_dir)
import settings

parser = argparse.ArgumentParser(description="Script for generating plots.")
parser.add_argument('--model', type=str, default='2DModel', help='Select the calculation mode defined in settings.py')
args = parser.parse_args()

if args.model not in settings.MODELS:
    print(f"\n[ERROR] Mode '{args.model}' is not defined in settings.py!")
    print(f"Available modes: {list(settings.MODELS.keys())}\n")
    sys.exit(1)

cfg = settings.MODELS[args.model]
print(f"--> Running {os.path.basename(__file__)} in mode: {args.model}")

VA_name = cfg.get('legend_variable_name', 'R')
VA_unit = cfg.get('legend_unit', '')

base_path = r"./DATAcut"
output_path = r"./Graphs/PhaseshiftsDefects"
os.makedirs(output_path, exist_ok=True)

potentials = cfg['potentials']

# Shifts
shifts = np.array(cfg['PhaseShiftsDefects_shifts']) * np.pi
shifts_H = np.array(cfg['PhaseShiftsDefects_shifts_H']) * np.pi          
shifts_PHP = np.array(cfg['PhaseShiftsDefects_shifts_PHP']) * np.pi          

aspect_ratio = 'square'

# Definice parametrů legendy
legend_params = {
    'frameon': False,
    'ncol': cfg.get('PhaseShiftsDefects_ncol', 2),              
    'fontsize': cfg.get('PhaseShiftsDefects_font_size',10),  
    'handlelength': 2,
    'title_fontsize': cfg.get('PhaseShiftsDefects_font_size',10),
    'loc': 'best',
    'columnspacing': 1.0,
    'handletextpad': 0.5
}

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 12
SIZE = 6

X_axis = r'Electron energy $\epsilon\,(\mathrm{eV})$'
Y_axis = r'Fixed nuclei phase shift $\delta(\epsilon)$'

limit_x = True 
limit_y = True 
x_rangeAll = [0,15]
y_rangeAll = cfg['PhaseShiftsDefects_y_rangeAll']  

plt.rcParams['font.family'] = 'serif'
plt.rcParams["font.serif"] = ["Latin Modern Roman"]
plt.rcParams['mathtext.fontset'] = 'cm'

plt.rcParams['axes.titlepad'] = 10 
plt.rcParams['axes.labelpad'] = 10 
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.edgecolor'] = "#000000"
plt.rcParams["figure.autolayout"] = True
plt.rcParams["legend.handlelength"] = 3
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.borderpad"] = 0.8

plt.rc('font', size=SMALL_SIZE) 
plt.rc('axes', titlesize=MEDIUM_SIZE) 
plt.rc('axes', labelsize=MEDIUM_SIZE) 
plt.rc('xtick', labelsize=SMALL_SIZE) 
plt.rc('ytick', labelsize=SMALL_SIZE) 
plt.rc('legend', fontsize=LEGEND_SIZE) 
plt.rc('figure', titlesize=BIGGER_SIZE) 

if aspect_ratio == 'golden_ratio':
    golden_ratio = (1 + 5 ** 0.5) / 2
    fig_all, ax_all = plt.subplots(figsize=(SIZE*golden_ratio, SIZE))
else:
    fig_all, ax_all = plt.subplots(figsize=(SIZE, SIZE))

for (potindx, VA) in enumerate(potentials):
    col_idx = potindx + 1 
    
    # 1. Total phase shift line
    try:
        data_ps = np.loadtxt(f"{base_path}/DSState/phaseshift.txt")
        x_ps = data_ps[:, 0]
        y_ps = data_ps[:, col_idx] + shifts[potindx]
        
        lbl_all = r'$\delta(\epsilon)$' if potindx == 0 else None
        ax_all.plot(x_ps, y_ps, label=lbl_all, color='blue', linestyle='-', linewidth=1.5)
    except OSError:
        print(f"Skipping phaseshift.txt (not found)")

    # 2. Total defects H
    try:
        h_eig = np.loadtxt(f"{base_path}/RydbergStates/eigE_H.txt")[:, col_idx]
        h_def = np.loadtxt(f"{base_path}/RydbergStates/defects_H.txt")[:, col_idx]
        h_def = (h_def * np.pi) + shifts_H[potindx] 
        
        lbl_H = r'$\pi\mu^H_n$' if potindx == 0 else None
        ax_all.plot(h_eig, h_def, label=lbl_H, color='purple', marker='x', linestyle='None', markersize=6, markeredgewidth=1)
    except OSError:
        pass

    # 3. Total defects PHP
    try:
        php_eig = np.loadtxt(f"{base_path}/RydbergStates/eigE_PHP.txt")[:, col_idx]
        php_def = np.loadtxt(f"{base_path}/RydbergStates/defects_PHP.txt")[:, col_idx]
        php_def = (php_def * np.pi) + shifts_PHP[potindx] 
        
        lbl_PHP = r'$\pi\mu^{PHP}_n$' if potindx == 0 else None
        ax_all.plot(php_eig, php_def, label=lbl_PHP, color='mediumaquamarine', marker='+', linestyle='None', markersize=6, markeredgewidth=1)
    except OSError:
        pass

ax_all.set_xlabel(X_axis)
ax_all.set_ylabel(Y_axis)

if limit_x: ax_all.set_xlim(x_rangeAll)

ax_all.legend(**legend_params)

ax_all.xaxis.set_minor_locator(AutoMinorLocator())
ax_all.yaxis.set_minor_locator(AutoMinorLocator())
ax_all.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax_all.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)
ax_all.grid(False)

fig_all.savefig(rf'{output_path}/PhaseshiftDefectsAll.pdf', format='pdf')
fig_all.savefig(rf'{output_path}/PhaseshiftDefectsAll.pgf', format='pgf')
plt.close(fig_all)
gc.collect()