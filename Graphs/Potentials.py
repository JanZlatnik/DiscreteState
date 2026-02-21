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


base_path = r"./DATAcut"
output_path = r"./Graphs/Potentials"
os.makedirs(output_path, exist_ok=True)

potentials = cfg['potentials']
aspect_ratio = 'square'

HARTREE_EV = 27.211386
Z = cfg['Potentials_Z']
l = cfg['Potentials_l']
mask_start = 5  

legend_params = {
    'frameon': False,
    'ncol': 1,              
    'fontsize': 12,  
    'handlelength': 2,
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

X_axis = r'Radial Coordinate$\,(a_0)$'
Y_axis_V = r'Potential$\,(\mathrm{eV})$'
Y_axis_Psi = r'Discrete State Wavefunction$\,(\mathrm{u.a.})$'

limit_x = True 
limit_y_V = True 
x_range = cfg['Potentials_x_range'] 
y_range_V = cfg['Potentials_y_range_V']    
y_range_scale = cfg['Potentials_y_range_scale']

color_V = 'green'
color_Asy = 'purple'
color_Harm = 'red'
color_Psi = 'blue'

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

ds_energies = np.loadtxt(f"{base_path}/DSState/DSenergies.txt", usecols=(1,))


for (potindx, VA) in enumerate(potentials):

    name = rf'Potential_{VA_name}={VA:.2f}'
    col_idx = potindx + 1 
    legend_title = f"${VA_name} = {VA:.2f}$"

    e_ds = ds_energies[potindx]

    try:
        data_V = np.loadtxt(f"{base_path}/DSState/V.txt")
        r_V = data_V[:, 0]
        V_vals = data_V[:, col_idx]
        
        r_plot = r_V[mask_start:]
        V_plot = V_vals[mask_start:]
    except OSError:
        print(f"Skipping V.txt for {VA_name}={VA}")
        continue

    try:
        data_psi = np.loadtxt(f"{base_path}/DSState/dstate.txt")
        r_psi = data_psi[:, 0]
        psi_vals = data_psi[:, col_idx]
    except OSError:
        psi_vals = None

    V_asy = (-Z/r_plot + (l*(l+1))/(2 * r_plot**2)) * HARTREE_EV

    r0, Omega, V0 = cfg['Potentials_harmonic_params'](VA)
    V_harm_au = V0 + 0.5 * (Omega**2) * ((r_plot - r0)**2)
    V_harm = V_harm_au * HARTREE_EV


    if aspect_ratio == 'golden_ratio':
        fig, ax1 = plt.subplots(figsize=(SIZE*golden_ratio + 1.0, SIZE))
    else:
        fig, ax1 = plt.subplots(figsize=(SIZE + 0.7, SIZE))

    ax1.plot(r_plot, V_plot, color=color_V, linestyle='-', linewidth=2)
    ax1.plot(r_plot, V_asy, color=color_Asy, linestyle='--', linewidth=1.2)
    ax1.plot(r_plot, V_harm, color=color_Harm, linestyle=':', linewidth=1.2)

    ax1.set_xlabel(X_axis)
    ax1.set_ylabel(Y_axis_V)
    if limit_x: ax1.set_xlim(x_range)
    if limit_y_V: ax1.set_ylim(y_range_V)
    
    ax2 = ax1.twinx()
    if psi_vals is not None:
        ax2.plot(r_psi, psi_vals, color=color_Psi, linestyle='-.', linewidth=1.5, alpha=1.0)
        ax2.fill_between(r_psi, 0, psi_vals, color=color_Psi, alpha=0.1) 
    
    ax2.set_ylabel(Y_axis_Psi)
    if limit_y_V: 
        V_min, V_max = y_range_V
        frac = (e_ds - V_min) / (V_max - V_min)
        ax2.set_ylim([-frac * y_range_scale, (1-frac) * y_range_scale]) 

    # Sjednocená legenda pomocí dummy lines na ax1
    lines = []
    labels = []
    
    lines.append(plt.Line2D([], [], color=color_V, linestyle='-', linewidth=2))
    labels.append(r'$V(r)$')

    
    if psi_vals is not None:
        lines.append(plt.Line2D([], [], color=color_Psi, linestyle='-.', linewidth=1.5))
        labels.append(r'$\phi_d(r)$')

    ax1.legend(lines, labels, title=legend_title, **legend_params)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.tick_params(axis='both', direction='in', which='major', top=True, length=6)
    ax2.tick_params(axis='y', direction='in', which='major', right=True, length=6)
    ax1.tick_params(axis='both', direction='in', which='minor', top=True, length=3)
    ax2.tick_params(axis='y', direction='in', which='minor', right=True, length=3)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())

    fig.savefig(f'{output_path}/{name}.pdf', format='pdf')
    plt.close(fig)
    gc.collect()

    
    ax_all.plot(r_plot, V_plot, color=color_V, linewidth=1.2, alpha=1.0)


if 'r_plot' in locals():
    V_asy_all = (-Z/r_plot + (l*(l+1))/(2 * r_plot**2)) * HARTREE_EV
    ax_all.plot(r_plot, V_asy_all, color=color_Asy, linestyle='--', linewidth=1.0, zorder=10)

ax_all.set_xlabel(X_axis)
ax_all.set_ylabel(Y_axis_V)
if limit_x: ax_all.set_xlim(x_range)
if limit_y_V: ax_all.set_ylim(y_range_V)

lines_all = [
    plt.Line2D([], [], color=color_V, linestyle='-', linewidth=1.2),
]
labels_all = [r'$V(r)$']
ax_all.legend(lines_all, labels_all, **legend_params)

ax_all.xaxis.set_minor_locator(AutoMinorLocator())
ax_all.yaxis.set_minor_locator(AutoMinorLocator())
ax_all.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax_all.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)
ax_all.xaxis.set_minor_locator(AutoMinorLocator())
ax_all.yaxis.set_minor_locator(AutoMinorLocator())

fig_all.savefig(rf'{output_path}/PotentialsAll.pdf', format='pdf')
plt.close(fig_all)
gc.collect()