import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os

# --- Nastavení ---
output_path = r"./Graphs/Drafts"
os.makedirs(output_path, exist_ok=True)

SIZE = 6
SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24

X_axis = r'$x$ axis'
Y_axis_V = r'$y$ axis'
Y_axis_Psi = r'$y$ secondary axis'

# Inicializace rcParams (stejné jako v tvém hlavním skriptu)
plt.rcParams['font.family'] = 'serif'
plt.rcParams["font.serif"] = ["Latin Modern Roman"]
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['axes.titlepad'] = 10 
plt.rcParams['axes.labelpad'] = 10 
plt.rcParams["figure.autolayout"] = True
plt.rc('font', size=SMALL_SIZE) 
plt.rc('axes', titlesize=MEDIUM_SIZE) 
plt.rc('axes', labelsize=MEDIUM_SIZE) 
plt.rc('xtick', labelsize=SMALL_SIZE) 
plt.rc('ytick', labelsize=SMALL_SIZE) 

# ==========================================
# 1. DRAFT GRAF - JEDNA OSA (jako PotentialsAll)
# ==========================================
# Nastavení velikosti stejné jako v fig_all (SIZE, SIZE)
fig1, ax1 = plt.subplots(figsize=(SIZE, SIZE))

ax1.set_xlim(1, 10)
ax1.set_ylim(1, 10)
ax1.set_xlabel(X_axis)
ax1.set_ylabel(Y_axis_V)

ax1.text(5.5, 5.5, 'DRAFT GRAPH', color='red', fontsize=20, 
         fontweight='bold', ha='center', va='center')

# Tick formatting
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.tick_params(axis='both', direction='in', which='major', top=True, right=True, length=6)
ax1.tick_params(axis='both', direction='in', which='minor', top=True, right=True, length=3)

fig1.savefig(rf'{output_path}/Draft_SingleAxis.pgf', format='pgf')
plt.close(fig1)

# ==========================================
# 2. DRAFT GRAF - DVĚ OSY (jako jednotlivé potenciály)
# ==========================================
# Nastavení velikosti stejné jako v jednotlivých grafech (SIZE + 0.7, SIZE)
fig2, ax2 = plt.subplots(figsize=(SIZE + 0.7, SIZE))

ax2.set_xlim(1, 10)
ax2.set_ylim(1, 10)
ax2.set_xlabel(X_axis)
ax2.set_ylabel(Y_axis_V)

ax2.text(5.5, 5.5, 'DRAFT GRAPH', color='red', fontsize=20, 
         fontweight='bold', ha='center', va='center')

ax3 = ax2.twinx()
ax3.set_ylim(1, 10)
ax3.set_ylabel(Y_axis_Psi)

# Tick formatting
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())
ax2.tick_params(axis='both', direction='in', which='major', top=True, length=6)
ax3.tick_params(axis='y', direction='in', which='major', right=True, length=6)
ax2.tick_params(axis='both', direction='in', which='minor', top=True, length=3)
ax3.tick_params(axis='y', direction='in', which='minor', right=True, length=3)

fig2.savefig(rf'{output_path}/Draft_DoubleAxis.pgf', format='pgf')
plt.close(fig2)

print(f"--> Draft graphs saved to {output_path}")