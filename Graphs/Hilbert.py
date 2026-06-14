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
output_path = r"./Graphs/HilbertTransform"
os.makedirs(output_path, exist_ok=True)

plot_gamma = cfg.get('Hilbert_plot_gamma', False)

potentials = cfg['potentials']
aspect_ratio = 'square'

# Definice parametrů legendy
legend_params = {
    'frameon': False,
    'ncol': 1,              
    'fontsize': 12,  
    'handlelength': 3.0,
    'title_fontsize': 12,
    'loc': 'best',
    'columnspacing': 1.0,
    'handletextpad': 0.5
}

filenames = [f"{base_path}/DSState/delta.txt",
             f"{base_path}/Hilbert/DeltaContinuous.txt",
             f"{base_path}/Hilbert/DeltaRydberg.txt",
             f"{base_path}/Hilbert/DeltaFull.txt"]

legends = [[r'$\Delta(\epsilon)$'],[r'$\Delta_\mathrm{cont.}(\epsilon)$'],[r'$\Delta_\mathrm{dis.}(\epsilon)$'],[r'$\Delta_\mathrm{trans.}(\epsilon)$']]
colors = [['blue'],['green'],['orange'],['red']]
line_styles = [['-'],['-.'],[':'],['--']]
line_widths = [[''],[''],[''],['']]

if plot_gamma:
    filenames.extend([f"{base_path}/DSState/gamma.txt", f"{base_path}/DSState/Vde.txt"])
    legends.extend([[r'$\Gamma(\epsilon)$'], [r'$2\pi|V_{d\epsilon}|^2$']])
    colors.extend([['black'], ['gray']])
    line_styles.extend([['-'], ['--']])
    line_widths.extend([[''], ['']])


# Inicializace All grafu
SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 12
SIZE = 6

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

    # Nastavení základních vlastností grafu
    name = rf'Hilbert_{VA_name}={VA:.2f}'
    x_column = 1
    y_columns = [potindx+2]
    show_legend = True 
    legend_title = f"${VA_name} = {VA:.2f}{VA_unit}$"
    
    # Nastavení velikostí
    SMALL_SIZE = 14
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 24
    LEGEND_SIZE = 12
    SIZE = 6

    # Legendy, popisky, barvy a styly čar
    X_axis = r'Electron energy $\epsilon\,(\mathrm{eV})$'
    Y_axis = r'Level shift $\Delta(\epsilon)\,(\mathrm{eV})$'
    if plot_gamma:
        Y_axis = r'Level shift $F(\epsilon)\,(\mathrm{eV})$'
    labels = [[0]] * len(filenames)
    labels_position = [[(0,0)]] * len(filenames)
    markers = [['o', '^']] * len(filenames) 
    
    ylog = False
    xlog = False

    # Nastavení pro omezení os a krokování
    limit_x = True 
    limit_y = True 
    x_range = cfg['Hilbert_x_range']
    y_range = cfg['Hilbert_y_range'] 


    # Nastavení LaTeX formátování a fontů
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

    # Nastavení velikosti a poměru stran figury
    if aspect_ratio == 'golden_ratio':
        golden_ratio = (1 + 5 ** 0.5) / 2
        fig, ax = plt.subplots(figsize=(SIZE*golden_ratio, SIZE ))
    else:
        fig, ax = plt.subplots(figsize=(SIZE, SIZE))

    # Načtení a vykreslení dat z více souborů
    for file_idx, filename in enumerate(filenames):
        try:
            data = np.loadtxt(filename, comments='#', delimiter=None)
            x_full = data[:, x_column-1]
            mask = (x_full >= x_range[0]) & (x_full <= x_range[1])
            x = x_full[mask]

            for col_idx, y_col in enumerate(y_columns):
                y_full = data[:, y_col-1]
                y = y_full[mask]
                
                if plot_gamma and file_idx == 5:
                    y = 2 * np.pi * (y**2)

                y[np.abs(y) > 10] = np.nan
                label = legends[file_idx][0] 
                
                if label:
                    if line_widths[file_idx][0] == '':
                        line, = ax.plot(x, y, label=label, color=colors[file_idx][0], linestyle=line_styles[file_idx][0])
                    else:
                        line, = ax.plot(x, y, label=label, color=colors[file_idx][0], linestyle=line_styles[file_idx][0], linewidth=line_widths[file_idx][0])
                else:
                    if line_widths[file_idx][0] == '':
                        line, = ax.plot(x, y, color=colors[file_idx][0], linestyle=line_styles[file_idx][0])
                    else:
                        line, = ax.plot(x, y, color=colors[file_idx][0], linestyle=line_styles[file_idx][0], linewidth=line_widths[file_idx][0])
                
                if labels[file_idx][0] != 0:
                    abs_x_pos = labels_position[file_idx][0][0] 
                    abs_y_pos = labels_position[file_idx][0][1] 
                    ax.text(abs_x_pos, abs_y_pos, labels[file_idx][0], color=colors[file_idx][0], verticalalignment='bottom')
                
                include_in_all = (file_idx == 0 or file_idx == 3)
                if plot_gamma and (file_idx == 4 or file_idx == 5):
                    include_in_all = True
                
                if include_in_all:
                    lbl_all = label if potindx == 0 else None
                    if line_widths[file_idx][0] == '':
                        ax_all.plot(x, y, label=lbl_all, color=colors[file_idx][0], linestyle=line_styles[file_idx][0], alpha=0.7)
                    else:
                        ax_all.plot(x, y, label=lbl_all, color=colors[file_idx][0], linestyle=line_styles[file_idx][0], linewidth=line_widths[file_idx][0], alpha=0.7)

        except OSError:
            pass


    # Nastavení vzhledu grafu
    ax.set_xlabel(X_axis)
    ax.set_ylabel(Y_axis)
    if limit_x:
        ax.set_xlim(x_range)
    if limit_y:
        ax.set_ylim(y_range)
    ax.set_title('')
    
    if show_legend:
        ax.legend(title=legend_title, **legend_params)
            
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='x', direction='in', which='major', top=True, bottom=True, labelbottom=True, length = 6)
    ax.tick_params(axis='y', direction='in', which='major', left=True, right=True, labelleft=True, length = 6)
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

    fig.savefig(f'{output_path}/{name}.pdf', format='pdf')
    fig.savefig(f'{output_path}/{name}.pgf', format='pgf')
    plt.close(fig)
    gc.collect()


# Uložení All grafu
ax_all.set_xlabel(X_axis)
ax_all.set_ylabel(Y_axis)
if limit_x: ax_all.set_xlim(x_range)
if limit_y: ax_all.set_ylim(y_range)

ax_all.legend(**legend_params)

ax_all.xaxis.set_minor_locator(AutoMinorLocator())
ax_all.yaxis.set_minor_locator(AutoMinorLocator())
ax_all.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax_all.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)
ax_all.grid(False)

fig_all.savefig(f'{output_path}/HilbertAll.pdf', format='pdf')
fig_all.savefig(f'{output_path}/HilbertAll.pgf', format='pgf')
plt.close(fig_all)
gc.collect()