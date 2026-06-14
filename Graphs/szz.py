import numpy as np
import matplotlib
matplotlib.use('Agg')  # Striktní headless mód pro Linux cluster
import matplotlib.pyplot as plt
import gc
import os

# Bezpečné nastavení pro PGF export bez volání lokálního LaTeXu
plt.rcParams.update({
    'pgf.rcfonts': False,
    'text.usetex': False,
    'font.family': 'serif'
})

# Převodní konstanta
a0 = 0.529177

# Osa R se vždy počítá od 0.1 do kus dál než je limit osy x
# --- 1. Valenční stav: H2 ($X\,^1\Sigma_g^+$) ---
R_data_val = np.linspace(0.1, 10.0, 1000) 
De_H2 = 4.75 # eV
Re_H2 = 0.74 / a0 # cca 1.4 a0
a_H2 = 1.93 * a0 # škálováno do Bohr
E_val = De_H2 * ( (1 - np.exp(-a_H2 * (R_data_val - Re_H2)))**2 - 1 )

# --- 2. Rydbergovy stavy a kation: H2 ($\rightarrow\ H_2^+$) ---
R_data_ryd = np.linspace(0.1, 10.0, 1000)
De_cat = 2.77 # eV
Re_cat = 1.06 / a0 # cca 2.0 a0
a_cat = 1.3 * a0
E_asymp_cat = 13.6 # eV
E_cat = De_cat * ( (1 - np.exp(-a_cat * (R_data_ryd - Re_cat)))**2 - 1 ) + E_asymp_cat

# --- 3. Iontová: LiF ($X\,^1\Sigma^+$) ---
R_data_ion = np.linspace(0.1, 60.0, 2000)
R_a0_ion = R_data_ion 
E_inf_ion = 1.99
E_ion = 1630 * np.exp(-4.58 * (R_data_ion * a0)) - 27.21 / R_a0_ion + E_inf_ion

# --- 4. Van der Waals: Ar2 ($X\,^1\Sigma_g^+$) ---
R_data_vdw = np.linspace(0.1, 20.0, 1000)
eps_Ar = 12.0 # meV
Rm_Ar = 3.8 / a0 # cca 7.18 a0
E_vdw = eps_Ar * ( (Rm_Ar / R_data_vdw)**12 - 2 * (Rm_Ar / R_data_vdw)**6 )

# --- Vykreslení ---
fig, axs = plt.subplots(2, 2, figsize=(11, 8.5), facecolor='white')

xlim_val_ryd = [0.6, 8.0]
xlim_ion = [0.1, 40.0]
xlim_vdw = [6.0, 15.0]

# --- Panel 1: Valenční (H2) ---
axs[0, 0].plot(R_data_val, E_val, color='#1f77b4', lw=2)
axs[0, 0].axhline(0, color='black', lw=0.8, ls='--')
axs[0, 0].set_ylim(-5.5, 2)
axs[0, 0].set_xlim(xlim_val_ryd)
axs[0, 0].set_title(r'Valenční vazba ($\mathrm{H}_2\ X\,^1\Sigma_g^+$)', fontsize=13)
axs[0, 0].set_ylabel('$E$ [eV]')
axs[0, 0].grid(False)

# --- Panel 2: Rydbergovy (H2 -> H2+) ---
ryd_color = 'green'
axs[0, 1].plot(R_data_ryd, E_cat, color='black', lw=2, label=r'$\mathrm{H}_2^+\ X\,^2\Sigma_g^+$ (Ionizační limita)')

for n in range(1, 11):
    if n == 1:
        E_Rn = E_cat - 13.6 / (1.1**2) 
    else:
        E_Rn = E_cat - 13.6 / (n**2)
    
    axs[0, 1].plot(R_data_ryd, E_Rn, color=ryd_color, ls='-.', lw=1.2, alpha=0.8)
    
    if n in [1, 2, 3]:
        idx = np.abs(R_data_ryd - 7.9).argmin()
        y_pos = E_Rn[idx]
        axs[0, 1].text(7.9, y_pos - 0.4, f'$n={n}$', color=ryd_color, ha='right', va='top', fontsize=11)

axs[0, 1].set_ylim(-6, 16)
axs[0, 1].set_xlim(xlim_val_ryd)
axs[0, 1].set_title(r'Rydbergova série ($\mathrm{H}_2\ \rightarrow\ \mathrm{H}_2^+$)', fontsize=13)
axs[0, 1].set_ylabel('$E$ [eV]')
axs[0, 1].legend(loc='lower left', fontsize=9)
axs[0, 1].grid(False)

# --- Panel 3: Iontová (LiF) ---
axs[1, 0].plot(R_data_ion, E_ion, color='purple', lw=2, label=r'Iontový stav $X\,^1\Sigma^+$')
axs[1, 0].axhline(0, color='black', lw=0.8, ls='--', label=r'Asymptota $\mathrm{Li} + \mathrm{F}$')
axs[1, 0].axhline(E_inf_ion, color='gray', lw=0.8, ls=':', label=r'Asymptota $\mathrm{Li}^+ + \mathrm{F}^-$')
axs[1, 0].set_ylim(-7, 4)
axs[1, 0].set_xlim(xlim_ion) 
axs[1, 0].set_title(r'Iontová vazba ($\mathrm{LiF}\ X\,^1\Sigma^+$)', fontsize=13)
axs[1, 0].set_xlabel(r'$R$ [Bohr]')
axs[1, 0].set_ylabel('$E$ [eV]')
axs[1, 0].legend(loc='lower right', fontsize=9)
axs[1, 0].grid(False)

# --- Panel 4: Van der Waals (Ar2) ---
axs[1, 1].plot(R_data_vdw, E_vdw, color='teal', lw=2)
axs[1, 1].axhline(0, color='black', lw=0.8, ls='--')
axs[1, 1].set_ylim(-15, 5) 
axs[1, 1].set_xlim(xlim_vdw)
axs[1, 1].set_title(r'Van der Waalsova interakce ($\mathrm{Ar}_2\ X\,^1\Sigma_g^+$)', fontsize=13)
axs[1, 1].set_xlabel(r'$R$ [Bohr]')
axs[1, 1].set_ylabel('$E$ [meV]')
axs[1, 1].grid(False)

plt.tight_layout()

# --- Ukládání a cleanup ---
output_path = "."
name = "four_molecular_states"

fig.savefig(f'{output_path}/{name}.pdf', format='pdf', bbox_inches='tight')
fig.savefig(f'{output_path}/{name}.pgf', format='pgf', bbox_inches='tight')

plt.close(fig)
gc.collect()