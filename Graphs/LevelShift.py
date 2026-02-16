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
output_path = r"./Graphs/LevelShift"
os.makedirs(output_path, exist_ok=True)

potentials = np.linspace(-0.15, 0.3, 10)
aspect_ratio = 'square'


legend_params = {
    'frameon': False,
    'ncol': 2,              
    'fontsize': 10,  
    'handlelength': 1.5,
    'title_fontsize': 10,
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

# Nastavení velikostí
SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 10 
SIZE = 6

# Definice souborů křivek (indexy odpovídají smyčce níže)
# 0: Delta, 1: Gamma, 2: DeltaA, 3: GammaA, 4: Vde
files_curves = [
    f"{base_path}/DSState/delta.txt",
    f"{base_path}/DSState/gamma.txt",
    f"{base_path}/DSState/deltaA.txt",
    f"{base_path}/DSState/gammaA.txt",
    f"{base_path}/DSState/Vde.txt"
]

files_eig = [
    f"{base_path}/RydbergStates/eigE_PHP.txt",
    f"{base_path}/RydbergStates/eigE_H.txt"
]

ds_energies_data = np.loadtxt(f"{base_path}/DSState/DSenergies.txt")

# Popisky os
X_axis = r'Energy$\,(\mathrm{eV})$'
Y_axis = r'Level Shift Operator$\,(\mathrm{eV})$' 

# Nastavení LaTeX
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

# Inicializace souhrnného grafu (All)
if aspect_ratio == 'golden_ratio':
    golden_ratio = (1 + 5 ** 0.5) / 2
    fig_all, ax_all = plt.subplots(figsize=(SIZE*golden_ratio, SIZE))
else:
    fig_all, ax_all = plt.subplots(figsize=(SIZE, SIZE))

# Limity os
limit_x = True 
limit_y = True 
x_range = [-10, 10]
y_range = [-4.1, 4] 
minstep_x = False 
minstep_y = False 
step_x = 2 
step_y = 1

# Barvy pro křivky
colors_curves = ['blue', 'red', 'teal', 'orange', 'green']
styles_curves = ['-', '-', '-.', '-.', '--'] 
labels_curves = [r'$\Delta(E)$', r'$\Gamma(E)$', r'$\Delta_\mathrm{anc}(E)$', r'$\Gamma_\mathrm{anc}(E)$', r'$2\pi|V_{de}|^2$']
widths_curves = [1.5, 1.5, 1.2, 1.2, 1.2]

for (potindx, VA) in enumerate(potentials):

    # Nastavení grafu pro jeden potenciál
    name = rf'LevelShift_VA={VA:.2f}'
    
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
            
            # Maskování rozsahu X
            mask = (x >= x_range[0]) & (x <= x_range[1])
            x_plot = x[mask]
            y_plot = y_raw[mask]

            # Speciální úprava pro Vde (index 4)
            if i == 4:
                y_plot = 2 * np.pi * (y_plot**2)
            
            # Oříznutí extrémních hodnot pro graf
            y_plot[np.abs(y_plot) > 20] = np.nan

            # Vykreslení do individuálního grafu
            ax.plot(x_plot, y_plot, label=labels_curves[i], color=colors_curves[i], 
                    linestyle=styles_curves[i], linewidth=widths_curves[i])
            
            # Vykreslení do souhrnného grafu (JEN Delta, Gamma a jejich analytika, bez Vde)
            if i < 4: 
                # Pro 'All' graf nebudeme dávat label ke každé čáře, aby se nepřeplnila legenda
                lbl = labels_curves[i] if potindx == 0 else None
                # Použijeme stejné barvy, možná trochu tenčí
                ax_all.plot(x_plot, y_plot, label=lbl, color=colors_curves[i], 
                           linestyle=styles_curves[i], linewidth=1.0, alpha=0.7)

        except OSError:
            print(f"Skipping {fname} (not found)")

    # 2. Vykreslení vertikálních čar (Vlastní energie) - JEN v individuálním grafu
    # PHP stavy (index 0)
    try:
        data_eig = np.loadtxt(files_eig[0])
        eig_vals = data_eig[:, col_idx]
        eig_vals = eig_vals[(eig_vals >= x_range[0]) & (eig_vals <= x_range[1])]
        
        ax.vlines(eig_vals, ymin=y_range[0], ymax=y_range[1], colors='gray', linestyles='solid', linewidth=0.75, alpha=0.6)
        ax.plot([], [], color='gray', linestyle='solid', linewidth=1.0, label=r'$E_n^{PHP}$')
    except: pass

    # H stavy (index 1)
    try:
        data_eig = np.loadtxt(files_eig[1])
        eig_vals = data_eig[:, col_idx]
        eig_vals = eig_vals[(eig_vals >= x_range[0]) & (eig_vals <= x_range[1])]
        
        ax.vlines(eig_vals, ymin=y_range[0], ymax=y_range[1], colors='black', linestyles='dotted', linewidth=0.75, alpha=0.6)
        ax.plot([], [], color='black', linestyle='dotted', linewidth=1.0, label=r'$E_n^{H}$')
    except: pass

    # 3. Vykreslení přímky E + E_DS - JEN v individuálním grafu
    try:
        # DS energie je ve druhém sloupci (index 1)
        E_DS = ds_energies_data[potindx, 1]
        
        # Vytvoření bodů pro přímku
        x_line = np.linspace(x_range[0], x_range[1], 100)
        y_line = x_line - E_DS
        
        ax.plot(x_line, y_line, color='purple', linestyle='--', linewidth=1.0, label=r'$E + E_{d}$')
    except IndexError: pass


    # Nastavení vzhledu individuálního grafu
    ax.set_xlabel(X_axis)
    ax.set_ylabel(Y_axis)
    
    if limit_x: ax.set_xlim(x_range)
    if limit_y: ax.set_ylim(y_range)
    
    # Legenda ve dvou sloupcích
    legend_title = f"$V_A = {VA:.2f}$"
    ax.legend(title=legend_title, **legend_params)
    
    ax.tick_params(axis='both', direction='in', which='both', top=True, right=True, length=6)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(False)

    fig.savefig(f'{output_path}/{name}.pdf', format='pdf')
    plt.close(fig)
    gc.collect()

# Nastavení a uložení souhrnného grafu (All)
ax_all.set_xlabel(X_axis)
ax_all.set_ylabel(Y_axis)
if limit_x: ax_all.set_xlim(x_range)
if limit_y: ax_all.set_ylim(y_range)

# Legenda pro All graf (jen základní veličiny)
# Protože jsme přidávali labely jen při potindx==0, stačí zavolat legend
ax_all.legend(**legend_all_params)

ax_all.tick_params(axis='both', direction='in', which='both', top=True, right=True, length=6)
ax_all.xaxis.set_minor_locator(AutoMinorLocator())
ax_all.yaxis.set_minor_locator(AutoMinorLocator())
ax_all.grid(False)

fig_all.savefig(f'{output_path}/LevelShift_All.pdf', format='pdf')
plt.close(fig_all)
gc.collect()