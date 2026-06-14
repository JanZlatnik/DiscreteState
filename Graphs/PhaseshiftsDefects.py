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
output_path = r"./Graphs/PhaseshiftsDefects"
os.makedirs(output_path, exist_ok=True)


potentials = cfg['potentials']

# Shifts
shifts = np.array(cfg['PhaseShiftsDefects_shifts']) * np.pi
DS_shifts = np.array(cfg['PhaseShiftsDefects_DS_shifts']) * np.pi          
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

files_lines = [
    f"{base_path}/DSState/phaseshift.txt",
    f"{base_path}/DSState/DS_phaseshift.txt",
    f"{base_path}/DSState/defect.txt",
    f"{base_path}/DSState/DS_defect.txt"
]

files_points = [
    (f"{base_path}/RydbergStates/eigE_H.txt", f"{base_path}/RydbergStates/defects_H.txt"),
    (f"{base_path}/RydbergStates/eigE_PHP.txt", f"{base_path}/RydbergStates/defects_PHP.txt")
]

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 12
SIZE = 6

X_axis = r'Electron energy $\epsilon\,(\mathrm{eV})$'
Y_axis = r'Phase shift $\delta(\epsilon)$'

legends = [
    r'$\delta(\epsilon)$', r'$\delta_\mathrm{res}(\epsilon)$', r'$\delta_\mathrm{bg}(\epsilon)$', 
    r'$\pi\mu(\epsilon)$', r'$\pi\mu_\mathrm{res}(\epsilon)$', r'$\pi\mu_\mathrm{bg}(\epsilon)$',
    r'$\pi\mu_n$', r'$\pi\mu^{\mathcal{P}}_n$'
]
labels = [0] * 8 
labels_position = [(0,0)] * 8 


colors = ['blue', 'red', 'green', 'purple', 'orange', 'mediumaquamarine', 'magenta', 'midnightblue']
colorsAll = ['blue', 'red', 'green', 'purple', 'orange', 'mediumaquamarine', 'purple', 'mediumaquamarine']
markers = [None, None, None, None, None, None, 'x', '+'] 
line_styles = ['-', '-', '-', '-', '-', '-', '', '']    
line_widths = [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 0, 0]
ms_value = 10  
mew_value = 1.5  
ms_value_all = 6
mew_value_all = 1

ylog = False
xlog = False

bg_colour = 'green'

limit_x = True 
limit_y = True 
x_range = cfg['PhaseShiftsDefects_x_range'] 
y_range = cfg['PhaseShiftsDefects_y_range'] 
x_rangeAll = cfg['PhaseShiftsDefects_x_rangeAll'] 
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

    name = rf'PhaseShiftsDefects{VA_name}={VA:.2f}'
    x_column = 1
    col_idx = potindx + 1 
    
    show_legend = True 
    legend_title = f"${VA_name} = {VA:.2f}{VA_unit}$"
    legend_frame = False 

    if aspect_ratio == 'golden_ratio':
        golden_ratio = (1 + 5 ** 0.5) / 2
        fig, ax = plt.subplots(figsize=(SIZE*golden_ratio, SIZE ))
    else:
        fig, ax = plt.subplots(figsize=(SIZE, SIZE))

    
    raw_data = []
    plot_data_list = [] 
    
    loaded_datasets = []
    
    for f in files_lines:
        try:
            loaded_datasets.append(np.loadtxt(f))
        except OSError:
            print(f"Skipping {f} (not found)")
            loaded_datasets.append(None)
            
    if loaded_datasets[0] is not None:
        y_ps = loaded_datasets[0][:, col_idx] + shifts[potindx]
        x_ps = loaded_datasets[0][:, 0]
        plot_data_list.append((x_ps, y_ps))
    else:
        plot_data_list.append(None)

    if loaded_datasets[1] is not None:
        y_ds_ps = loaded_datasets[1][:, col_idx] + DS_shifts[potindx]
        x_ds = loaded_datasets[1][:, 0]
        plot_data_list.append((x_ds, y_ds_ps))
    else:
        plot_data_list.append(None)
        
    if loaded_datasets[0] is not None and loaded_datasets[1] is not None:
        y_ps_bg = y_ps - y_ds_ps
        plot_data_list.append((x_ps, y_ps_bg))
    else:
        plot_data_list.append(None)

    if loaded_datasets[2] is not None:
        y_def = loaded_datasets[2][:, col_idx] + shifts[potindx]
        x_def = loaded_datasets[2][:, 0]
        plot_data_list.append((x_def, y_def))
    else:
        plot_data_list.append(None)

    if loaded_datasets[3] is not None:
        y_ds_def = loaded_datasets[3][:, col_idx] + DS_shifts[potindx]
        x_ds_def = loaded_datasets[3][:, 0]
        plot_data_list.append((x_ds_def, y_ds_def))
    else:
        plot_data_list.append(None)
        
    if loaded_datasets[2] is not None and loaded_datasets[3] is not None:
        y_def_bg = y_def - y_ds_def
        plot_data_list.append((x_def, y_def_bg))
    else:
        plot_data_list.append(None)

    try:
        h_eig = np.loadtxt(files_points[0][0])[:, col_idx]
        h_def = np.loadtxt(files_points[0][1])[:, col_idx]
        h_def = (h_def * np.pi) + shifts_H[potindx] 
        plot_data_list.append((h_eig, h_def))
    except OSError:
        plot_data_list.append(None)
    
    try:
        php_eig = np.loadtxt(files_points[1][0])[:, col_idx]
        php_def = np.loadtxt(files_points[1][1])[:, col_idx]
        php_def = (php_def * np.pi) + shifts_PHP[potindx] 
        plot_data_list.append((php_eig, php_def))
    except OSError:
        plot_data_list.append(None)


    for i, data_pair in enumerate(plot_data_list):
        if data_pair is None:
            continue
            
        xx, yy = data_pair
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
                ax.plot(xx, yy, label=label, color=color, marker=mk, linestyle='None', markersize=ms_value, markeredgewidth=mew_value)
            else:
                ax.plot(xx, yy, color=color, marker=mk, linestyle='None', markersize=ms_value, markeredgewidth=mew_value)
                
        if labels[i] != 0:
            abs_x_pos = labels_position[i][0]
            abs_y_pos = labels_position[i][1]
            ax.text(abs_x_pos, abs_y_pos, labels[i], color=color, verticalalignment='bottom')

        lbl_all = label if (potindx == 0) else None
        if ls != '':
            color = colorsAll[i]
            ax_all.plot(xx, yy, label=lbl_all, color=color, linestyle=ls, linewidth=lw) 
        else:
            color = colorsAll[i]
            ax_all.plot(xx, yy, label=lbl_all, color=color, marker=mk, linestyle='None', markersize=ms_value_all, markeredgewidth=mew_value_all)


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

    ax.tick_params(axis='x', direction='in', which='major', top=True, bottom=True, labelbottom=True, length = 6)
    ax.tick_params(axis='y', direction='in', which='major', left=True, right=True, labelleft=True, length = 6)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='x', direction='in', which='minor', top=True, bottom=True, labelbottom=True, length = 3)
    ax.tick_params(axis='y', direction='in', which='minor', left=True, right=True, labelleft=True, length = 3)
    ax.grid(False)
    if ylog:
        ax.set_yscale('log')
    if xlog:
        ax.set_xscale('log')
    if limit_x:
        ax.set_xlim(x_range)
    if limit_y:
        ax.set_ylim(y_range)

    fig.savefig(rf'{output_path}/{name}.pdf', format='pdf')
    fig.savefig(rf'{output_path}/{name}.pgf', format='pgf')
    
    plt.close(fig)
    gc.collect()


ax_all.set_xlabel(X_axis)
ax_all.set_ylabel(Y_axis)
if limit_x: ax_all.set_xlim(x_rangeAll)
if limit_y: ax_all.set_ylim(y_rangeAll)

ax_all.legend(**legend_params)

ax_all.xaxis.set_minor_locator(AutoMinorLocator())
ax_all.yaxis.set_minor_locator(AutoMinorLocator())
ax_all.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax_all.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)
ax_all.xaxis.set_minor_locator(AutoMinorLocator())
ax_all.yaxis.set_minor_locator(AutoMinorLocator())
ax_all.grid(False)

fig_all.savefig(rf'{output_path}/PhaseshiftDefectsAll.pdf', format='pdf')
fig_all.savefig(rf'{output_path}/PhaseshiftDefectsAll.pgf', format='pgf')
plt.close(fig_all)
gc.collect()