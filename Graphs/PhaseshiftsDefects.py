import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import AutoMinorLocator
import gc

base_path = r"/work/zlatnikj/DiscreteStateRuns/2D_Test1/DATA"


potentials = np.linspace(0.7, 10.0, 10)
shifts = [-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1]
shifts = np.array(shifts) * np.pi
aspect_ratio = 'square'
specific_value = None  # Specifická hodnota, None pro vykreslení všech dat
filenames = [f"{base_path}/Coulomb/phaseshift.txt",
             f"{base_path}/Coulomb/DS_phaseshift.txt"]

eigefile = f"{base_path}/BoundStates/eigE.txt"
defectsfile = f"{base_path}/BoundStates/defects.txt"
dmarker = 'x'

for (potindx, VA) in enumerate(potentials):

    # Nastavení základních vlastností grafu
    name = rf'PhaseShiftsDefectsVA={VA:.2f}'
    x_column = 1
    y_columns = [potindx+2]
    show_legend = True  # True pro zobrazení, False pro skrytí legendy
    legend_title = f"$V_A = {VA:.2f}$"  # Nadpis nad legendou, 0 pro skrytí
    legend_frame = False # Rámeček kolem legendy

    # Nastavení velikostí
    SMALL_SIZE = 14
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 24
    LEGEND_SIZE = 12
    SIZE = 6

    # Legendy, popisky, barvy a styly čar s escape sekvencemi pro LaTeX
    X_axis = r'Energy$\,(\mathrm{eV})$'
    Y_axis = r'Phaseshift'
    legends = [[r'$\delta(E)$'],[r'$\delta_\mathrm{res}(E)$'],[r'$\delta_\mathrm{bg}(E)$'],[r'$\pi\mu(E)$']]
    labels = [[0],[0],[0],[0]]
    labels_position = [[(0,0)],[(0,0)],[(0,0)],[(0,0)]]  
    colors = [['blue'],['red'],['green'],['purple']]
    markers = [['o', '^']]
    line_styles = [['-'],['-'],['-'],['-']]
    line_widths = [[''],[''],[''],['']]
    ylog = False
    xlog = False

    bg_colour = 'green'

    # Nastavení pro omezení os a krokování
    limit_x = True  # True pro omezení osy X, False pro ponechání plného rozsahu
    limit_y = True  # True pro omezení osy Y, False pro ponechání plného rozsahu
    x_range = [-10, 10]  # Rozsah pro omezení osy X
    y_range = [-2*np.pi-0.1, 0.1]  # Rozsah pro omezení osy Y
    minstep_x = False  # Zapnutí/vypnutí minimálního kroku na ose X
    minstep_y = False  # Zapnutí/vypnutí minimálního kroku na ose Y
    step_x = 2  # Minimální krok na ose X
    step_y = 10  # Minimální krok na ose Y


    # Nastavení LaTeX formátování a fontů
    #rcParams['text.usetex'] = True
    #rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    #rcParams['font.serif'] = ['Latin Modern Roman']
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

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=LEGEND_SIZE)            # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # Nastavení velikosti a poměru stran figury
    if aspect_ratio == 'golden_ratio':
        golden_ratio = (1 + 5 ** 0.5) / 2
        fig, ax = plt.subplots(figsize=(SIZE*golden_ratio, SIZE ))
    else:
        fig, ax = plt.subplots(figsize=(SIZE, SIZE))

    # Načtení a vykreslení dat z více souborů
    for file_idx, filename in enumerate(filenames):
        data = np.loadtxt(filename, comments='#', delimiter=None)
        if file_idx == 0:
            data[:,y_columns[0]-1] = data[:,y_columns[0]-1] + shifts[potindx]
        x = data[:, x_column-1]

        for col_idx, y_col in enumerate(y_columns):
            y = data[:, y_col-1]
            label = legends[file_idx][col_idx]
            if label:
                if line_widths[file_idx][col_idx] == '':
                    line, = ax.plot(x, y, label=label, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx])
                else:
                    line, = ax.plot(x, y, label=label, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], linewidth=line_widths[file_idx][col_idx])
            else:
                if line_widths[file_idx][col_idx] == '':
                    line, = ax.plot(x, y, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx])
                else:
                    line, = ax.plot(x, y, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], linewidth=line_widths[file_idx][col_idx])
            if labels[file_idx][col_idx] != 0:
                # Použijeme absolutní hodnoty pro pozici popisků
                abs_x_pos = labels_position[file_idx][col_idx][0]  # Absolutní pozice na ose X
                abs_y_pos = labels_position[file_idx][col_idx][1]  # Absolutní pozice na ose Y
                ax.text(abs_x_pos, abs_y_pos, labels[file_idx][col_idx], color=colors[file_idx][col_idx], verticalalignment='bottom')


    data1 = np.loadtxt(filenames[0])
    data2 = np.loadtxt(filenames[1])
    x = data1[:,x_column-1]
    for col_idx, y_col in enumerate(y_columns):
        y = data1[:,y_col-1]-data2[:,y_col-1] + shifts[potindx]
        file_idx = 2
        label = legends[file_idx][col_idx]
        if label:
            if line_widths[file_idx][col_idx] == '':
                line, = ax.plot(x, y, label=label, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx])
            else:
                line, = ax.plot(x, y, label=label, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], linewidth=line_widths[file_idx][col_idx])
        else:
            if line_widths[file_idx][col_idx] == '':
                line, = ax.plot(x, y, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx])
            else:
                line, = ax.plot(x, y, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], linewidth=line_widths[file_idx][col_idx])
        if labels[file_idx][col_idx] != 0:
                # Použijeme absolutní hodnoty pro pozici popisků
                abs_x_pos = labels_position[file_idx][col_idx][0]  # Absolutní pozice na ose X
                abs_y_pos = labels_position[file_idx][col_idx][1]  # Absolutní pozice na ose Y
                ax.text(abs_x_pos, abs_y_pos, labels[file_idx][col_idx], color=colors[file_idx][col_idx], verticalalignment='bottom')


    x = np.loadtxt(eigefile)[:,potindx+1]
    y = np.pi * np.loadtxt(defectsfile)[:,potindx+1]
    file_idx = 3
    label = legends[file_idx][col_idx]
    col_idx = 0
    if label:
        if line_widths[file_idx][col_idx] == '':
            line, = ax.plot(x, y, label=label, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], marker=dmarker)
        else:
            line, = ax.plot(x, y, label=label, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], linewidth=line_widths[file_idx][col_idx], marker=dmarker)
    else:
        if line_widths[file_idx][col_idx] == '':
            line, = ax.plot(x, y, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], marker=dmarker)
        else:
            line, = ax.plot(x, y, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], linewidth=line_widths[file_idx][col_idx], marker=dmarker)
    if labels[file_idx][col_idx] != 0:
            # Použijeme absolutní hodnoty pro pozici popisků
            abs_x_pos = labels_position[file_idx][col_idx][0]  # Absolutní pozice na ose X
            abs_y_pos = labels_position[file_idx][col_idx][1]  # Absolutní pozice na ose Y
            ax.text(abs_x_pos, abs_y_pos, labels[file_idx][col_idx], color=colors[file_idx][col_idx], verticalalignment='bottom')



    # Extra funkce fit
    if False:
        x_range = np.linspace(np.min(x), np.max(x), 300)
        f_x = 1/(2*1822.89) * (np.pi/12)**2 * (x_range+11)**2 * 27.2113961
        ax.plot(x_range, f_x, color='black', linestyle='--', linewidth=1, alpha=0.7)
        relative_x_position = 0.9  # 0=začátek, 1=konec
        relative_y_shift = -0.8  # Vertikální posun v procentech rozsahu y
        label_pos_x = x_range[int(relative_x_position * len(x_range))]
        label_pos_y = f_x[int(relative_x_position * len(f_x))] + (max(f_x) - min(f_x)) * relative_y_shift
        ax.text(label_pos_x, label_pos_y, r'$\propto\frac{j^2\pi^2\hbar^2}{2mL^2}$', color='black', verticalalignment='bottom', horizontalalignment='center', alpha=0.7)

    # Nastavení vzhledu grafu
    ax.set_xlabel(X_axis)
    ax.set_ylabel(Y_axis)
    if limit_x:
        ax.set_xlim(x_range)
    if limit_y:
        ax.set_ylim(y_range)
    ax.set_title('')
    if show_legend:
        if legend_title != 0:
            ax.legend(frameon=legend_frame, title=legend_title)
        else:
            ax.legend(frameon=legend_frame)
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

    fig.savefig(rf'2D/Defects/{name}.pdf', format='pdf')
