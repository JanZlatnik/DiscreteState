import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
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

base_path = cfg['LevelShiftFar_base_path']
output_path = r"./Graphs/LevelShiftFar"
os.makedirs(output_path, exist_ok=True)

potentials = cfg['potentials']
aspect_ratio = 'square'


legend_params = {
    'frameon': False,
    'ncol': 1,              
    'fontsize': 12,  
    'handlelength': 1.5,
    'title_fontsize': 12,
    'loc': 'best',
    'columnspacing': 1.0,
    'handletextpad': 0.5
}

legend_all_params = {
    'frameon': False,
    'ncol': 1,
    'fontsize': 12,
    'handlelength': 2.0,
    'title_fontsize': 12,
    'loc': 'best',
    'columnspacing': 1.0,
    'handletextpad': 0.5
}

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 12 
SIZE = 6

# Definice souborů křivek
files_curves = [
    f"{base_path}/DSState/delta.txt",
    f"{base_path}/DSState/gamma.txt",
    f"{base_path}/DSState/Vde.txt"
]

X_axis = r'Electron energy $\epsilon\,(\mathrm{eV})$'
Y_axis = r'Level shift $F(\epsilon)\,(\mathrm{eV})$' 

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

limit_x = True 
limit_y = True 
x_range = cfg['LevelShiftFar_x_range']
y_range = cfg['LevelShiftFar_y_range'] 


colors_curves = ['blue', 'red', 'green']
styles_curves = ['-', '-', '--'] 
labels_curves = [r'$\Delta(\epsilon)$', r'$\Gamma(\epsilon)$', r'$2\pi|V_{d\epsilon}|^2$']
widths_curves = [1.5, 1.5, 1.2]

for (potindx, VA) in enumerate(potentials):

    name = rf'LevelShiftFar_{VA_name}={VA:.2f}'
    
    col_idx = potindx + 1

    if aspect_ratio == 'golden_ratio':
        golden_ratio = (1 + 5 ** 0.5) / 2
        fig, ax = plt.subplots(figsize=(SIZE*golden_ratio, SIZE ))
    else:
        fig, ax = plt.subplots(figsize=(SIZE, SIZE))

    for i, fname in enumerate(files_curves):
        try:
            data = np.loadtxt(fname)
            x = data[:, 0]
            y_raw = data[:, col_idx]
            
            mask = (x >= x_range[0]) & (x <= x_range[1])
            x_plot = x[mask]
            y_plot = y_raw[mask]

            if i == 2: # Vde je teď na indexu 2
                y_plot = 2 * np.pi * (y_plot**2)
            
            y_plot[np.abs(y_plot) > 50] = np.nan

            ax.plot(x_plot, y_plot, label=labels_curves[i], color=colors_curves[i], linestyle=styles_curves[i], linewidth=widths_curves[i])
            
            if i < 2: # Jen Delta a Gamma do All
                lbl = labels_curves[i] if potindx == 0 else None
                ax_all.plot(x_plot, y_plot, label=lbl, color=colors_curves[i], linestyle=styles_curves[i], linewidth=1.0, alpha=0.7)

        except OSError:
            print(f"Skipping {fname} (not found)")


    ax.set_xlabel(X_axis)
    ax.set_ylabel(Y_axis)
    
    if limit_x: ax.set_xlim(x_range)
    if limit_y: ax.set_ylim(y_range)
    
    legend_title = f"${VA_name} = {VA:.2f}{VA_unit}$"
    ax.legend(title=legend_title, **legend_params)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
    ax.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)
    ax.grid(False)

    fig.savefig(f'{output_path}/{name}.pdf', format='pdf')
    fig.savefig(f'{output_path}/{name}.pgf', format='pgf')
    plt.close(fig)
    gc.collect()


ax_all.set_xlabel(X_axis)
ax_all.set_ylabel(Y_axis)
if limit_x: ax_all.set_xlim(x_range)
if limit_y: ax_all.set_ylim(y_range)

ax_all.legend(**legend_all_params)

ax_all.xaxis.set_minor_locator(AutoMinorLocator())
ax_all.yaxis.set_minor_locator(AutoMinorLocator())
ax_all.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax_all.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)
ax_all.grid(False)

fig_all.savefig(f'{output_path}/LevelShiftFar_All.pdf', format='pdf')
fig_all.savefig(f'{output_path}/LevelShiftFar_All.pgf', format='pgf')
plt.close(fig_all)
gc.collect()