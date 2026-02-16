import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import AutoMinorLocator
import gc
import os

base_path = r"./DATAcut"
output_path = r"./Graphs/HilbertTransform"
os.makedirs(output_path, exist_ok=True)

plot_gamma = False

potentials = np.linspace(-0.15, 0.3, 10)
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

filenames = [f"{base_path}/Hilbert/DeltaRydberg.txt",
             f"{base_path}/DSState/delta.txt",
             f"{base_path}/Hilbert/DeltaFull.txt",
             f"{base_path}/Hilbert/DeltaContinuous.txt"]

legends = [[r'$\Delta_n(E)$'],[r'$\Delta(E)$'],[r'$\Delta_\mathrm{trans}(E)$'],[r'$\Delta_\epsilon(E)$']]
colors = [['orange'],['blue'],['red'],['green']]
line_styles = [[':'],['-'],['--'],['-.']]
line_widths = [[''],[''],[''],['']]

if plot_gamma:
    filenames.extend([f"{base_path}/DSState/gamma.txt", f"{base_path}/DSState/Vde.txt"])
    legends.extend([[r'$\Gamma(E)$'], [r'$2\pi|V_{d\epsilon}|^2$']])
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
    name = rf'Hilbert10VA={VA:.2f}'
    x_column = 1
    y_columns = [potindx+2]
    show_legend = True 
    legend_title = f"$V_A = {VA:.2f}$" 
    
    # Nastavení velikostí
    SMALL_SIZE = 14
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 24
    LEGEND_SIZE = 12
    SIZE = 6

    # Legendy, popisky, barvy a styly čar
    X_axis = r'Energy$\,(\mathrm{eV})$'
    Y_axis = r'Level Shift$\,(\mathrm{eV})$'
    labels = [[0]] * len(filenames)
    labels_position = [[(0,0)]] * len(filenames)
    markers = [['o', '^']] * len(filenames) 
    
    ylog = False
    xlog = False

    # Nastavení pro omezení os a krokování
    limit_x = True 
    limit_y = True 
    x_range = [-10, 10] 
    y_range = [-4.1, 4] 
    minstep_x = False 
    minstep_y = False 
    step_x = 1 
    step_y = 10 


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
                
                include_in_all = (file_idx == 1 or file_idx == 2)
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
            
    ax.tick_params(axis='x', direction='in', which='both', top=True, bottom=True, labelbottom=True, length = 6)
    ax.tick_params(axis='y', direction='in', which='both', left=True, right=True, labelleft=True, length = 6)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='x', direction='in', which='minor', top=True, bottom=True, labelbottom=True, length = 2)
    ax.tick_params(axis='y', direction='in', which='minor', left=True, right=True, labelleft=True, length = 2)
    ax.grid(False)
    if ylog:
        ax.set_yscale('log')
    if xlog:
        ax.set_xscale('log')
    if limit_x:
        ax.set_xlim(x_range)
        if minstep_x:
            ax.xaxis.set_major_locator(plt.MultipleLocator(step_x))
    if limit_y:
        ax.set_ylim(y_range)
        if minstep_y:
            ax.yaxis.set_major_locator(plt.MultipleLocator(step_y))

    fig.savefig(f'{output_path}/{name}.pdf', format='pdf')
    plt.close(fig)
    gc.collect()


# Uložení All grafu
ax_all.set_xlabel(X_axis)
ax_all.set_ylabel(Y_axis)
if limit_x: ax_all.set_xlim(x_range)
if limit_y: ax_all.set_ylim(y_range)

ax_all.legend(**legend_params)

ax_all.tick_params(axis='both', direction='in', which='both', top=True, right=True, length=6)
ax_all.xaxis.set_minor_locator(AutoMinorLocator())
ax_all.yaxis.set_minor_locator(AutoMinorLocator())
ax_all.grid(False)

fig_all.savefig(f'{output_path}/HilbertAll.pdf', format='pdf')
plt.close(fig_all)
gc.collect()