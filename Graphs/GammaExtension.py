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

base_path = r"./DATAcut"
output_path = r"./Graphs/GammaExtension"
os.makedirs(output_path, exist_ok=True)

potentials = cfg['potentials']
max_points_to_plot = 35

aspect_ratio = 'square'

# Definice parametrů legendy
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

files_lines = [
    f"{base_path}/DSState/gamma.txt",
    f"{base_path}/DSState/gammaA.txt",
    f"{base_path}/DSState/defect.txt",
    f"{base_path}/DSState/DS_defect.txt",
]

files_points = [
    f"{base_path}/RydbergStates/eigE_PHP.txt",
    f"{base_path}/RydbergStates/Vdn.txt"
]

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 12
SIZE = 6

X_axis = r'Electron energy $\epsilon\,(\mathrm{eV})$'
Y_axis = r'Resonance width $\Gamma(\epsilon)\,(\mathrm{eV})$' 

legends = [
    r'$\Gamma(\epsilon)$',
    r'$\Gamma_\mathrm{anc}(\epsilon)$',
    r'$2\pi|V_{dn}|^2\rho(\epsilon)$'
]

colors = ['blue', 'red', 'green'] 
markers = [None, None, 'x'] 
line_styles = ['-', '-', '']    
line_widths = [1.5, 1.5, 0]
ms_value = 8
mew_value = 1.2
ms_value_all = 6
mew_value_all = 1

ylog = False 
xlog = False

limit_x = True 
limit_y = True 
x_range = cfg['GammaExtension_x_range']
y_range = cfg['GammaExtension_y_range']     
x_rangeAll = x_range
y_rangeAll = y_range


RYDBERG_EV = 13.605693
HARTREE_EV = 27.211386

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

    name = rf'GammaExtension_{VA_name}={VA:.2f}'
    x_column = 0 
    col_idx = potindx + 1 
    
    show_legend = True 
    legend_title = f"${VA_name} = {VA:.2f}{VA_unit}$"
    
    if aspect_ratio == 'golden_ratio':
        golden_ratio = (1 + 5 ** 0.5) / 2
        fig, ax = plt.subplots(figsize=(SIZE*golden_ratio, SIZE ))
    else:
        fig, ax = plt.subplots(figsize=(SIZE, SIZE))

    
    raw_data = [np.loadtxt(f) for f in files_lines]
    
    y_gamma = raw_data[0][:, col_idx]
    x_gamma = raw_data[0][:, 0] 
    
    y_gammaA = raw_data[1][:, col_idx]
    x_gammaA = raw_data[1][:, 0]

    y_defect_cont = raw_data[2][:, col_idx] - raw_data[3][:, col_idx]
    x_defect_cont = raw_data[2][:, 0]
    
    
    energies_discrete = np.loadtxt(files_points[0])[:, col_idx]
    vdn_discrete = np.loadtxt(files_points[1])[:, col_idx]
    
    mask = energies_discrete < 0
    En = energies_discrete[mask]
    Vdn = vdn_discrete[mask] 

    if len(En) > max_points_to_plot:
        En = En[:max_points_to_plot]
        Vdn = Vdn[:max_points_to_plot]
    
    nu_n = np.sqrt(RYDBERG_EV / np.abs(En))
    
    dmu_dE = np.zeros_like(En)
    
    for i, e_val in enumerate(En):
        idx = (np.abs(x_defect_cont - e_val)).argmin()
        
        if idx == 0:
            idx = 1
        elif idx == len(x_defect_cont) - 1:
            idx = len(x_defect_cont) - 2
            
        dmu = (y_defect_cont[idx+1] - y_defect_cont[idx-1])/np.pi
        de  = x_defect_cont[idx+1] - x_defect_cont[idx-1]
        
        if de != 0:
            dmu_dE[i] = dmu / de
        else:
            dmu_dE[i] = 0 

    rho = (nu_n**3 / HARTREE_EV) + dmu_dE
    
    y_gamma_points = 2 * np.pi * (Vdn**2) * rho
    
    x_gamma_points = En


    plot_data_list = [
        (x_gamma, y_gamma),
        (x_gammaA, y_gammaA),
        (x_gamma_points, y_gamma_points)
    ]
    

    for i, (xx, yy) in enumerate(plot_data_list):
        label = legends[i]
        color = colors[i]
        ls = line_styles[i]
        mk = markers[i]
        lw = line_widths[i]
        
        if ls != '':
            if label:
                ax.plot(xx, yy, label=label, color=color, linestyle=ls, linewidth=lw)
            else:
                ax.plot(xx, yy, color=color, linestyle=ls, linewidth=lw)
        else:
            if label:
                ax.plot(xx, yy, label=label, color=color, marker=mk, linestyle='None', markersize=ms_value, markeredgewidth=mew_value, markerfacecolor='none') 
            else:
                ax.plot(xx, yy, color=color, marker=mk, linestyle='None', markersize=ms_value, markeredgewidth=mew_value, markerfacecolor='none')
        

        lbl_all = label if (potindx == 0) else None
        if ls != '':
            color = colors[i]
            ax_all.plot(xx, yy, label=lbl_all, color=color, linestyle=ls, linewidth=lw) 
        else:
            color = colors[i]
            ax_all.plot(xx, yy, label=lbl_all, color=color, marker=mk, linestyle='None', markersize=ms_value_all, markeredgewidth=mew_value_all, markerfacecolor='none')


    ax.set_xlabel(X_axis)
    ax.set_ylabel(Y_axis)
    if limit_x:
        ax.set_xlim(x_range)
    if limit_y:
        ax.set_ylim(y_range)
    ax.set_title('')
    if show_legend:
        if legend_title != 0:
            ax.legend(title=legend_title, **legend_params)
        else:
            ax.legend(**legend_params)
            
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
    ax.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    ax.axvline(0, color='gray', linestyle=':', linewidth=0.8) 

    ax.grid(False)
    if ylog:
        ax.set_yscale('log')
    if xlog:
        ax.set_xscale('log')
        
    fig.savefig(rf'{output_path}/{name}.pdf', format='pdf')
    fig.savefig(rf'{output_path}/{name}.pgf', format='pgf')
    
    plt.close(fig)
    gc.collect()


ax_all.set_xlabel(X_axis)
ax_all.set_ylabel(Y_axis)
if limit_x: ax_all.set_xlim(x_rangeAll)
if limit_y: ax_all.set_ylim(y_rangeAll)
ax_all.axvline(0, color='gray', linestyle=':', linewidth=0.8)

ax_all.legend(**legend_params)

ax_all.xaxis.set_minor_locator(AutoMinorLocator())
ax_all.yaxis.set_minor_locator(AutoMinorLocator())
ax_all.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax_all.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)
ax_all.xaxis.set_minor_locator(AutoMinorLocator())
ax_all.yaxis.set_minor_locator(AutoMinorLocator())
ax_all.grid(False)

fig_all.savefig(rf'{output_path}/GammaExtensionAll.pdf', format='pdf')
fig_all.savefig(rf'{output_path}/GammaExtensionAll.pgf', format='pgf')
plt.close(fig_all)
gc.collect()