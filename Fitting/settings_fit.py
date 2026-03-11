import numpy as np

MODELS = {
    '2DModel2': {
        'Fit_R0': 2.0,                              # Referencial for fitting, i.e. S(R) = Gamma(E,R)/Gamma(E,R0)
        'base_path': r"./DATA/DSState",
        'output_path': r"./Graphs/Fitting",
        'color_data': 'black',
        'color_fit': 'red',
        'aspect_ratio': 'square',


        # R-scyling factor fitting
        'p0_guess_g': [1.0, 2.0, 1.0, 5.0, 0.5],      
        'X_axis_g': r'Internuclear distance $R\,(a_0)$',
        'Y_axis_g': r'Scaling factor $g^2(R)$',
 

        # Defect fitting
        'defect_path': r"./DATA/RydbergStates",
        'X_axis_defect': r'Energy $(\mathrm{eV})$',
        'Y_axis_defect': r'Background quantum defect $\delta_\mathrm{bg}$',
        'p0_guess_defect': [0.0, 0.0, -2.0],

        # Gamma, Delta energy fitting
        'X_axis_E':  r'Energy $E\,(\mathrm{eV})$',
        'Y_axis_gamma': r'$\Gamma(R_0, E)\, (\mathrm{eV})$',
        'Y_axis_delta': r'$\Delta(R_0, E)\, (\mathrm{eV})$',
        'p0_guess_gamma': [0.01, 0.005, 5.0, 2.5, 13.0, 4.0, 20.0],
        'y_lim_delta': [-1.0, 1.0],


        # DS energy fitting
        'X_axis_R': r'Internuclear distance $R\,(a_0)$',
        'Y_axis_E': r'Discrete state energy $\epsilon_d\,(\mathrm{eV})$',
        'p0_guess_ds_energy': [-15.0, -2.0, 1.0, 2.0, 1.0, 0.5, 2.0],


        # bg phaseshift fitting
        'p0_guess_bg_phase': [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.5],
        'X_axis_bgE': r'Electron energy $\epsilon\,(\mathrm{eV})$',
        'Y_axis_bgPhase': r'Background phase shift $\delta_\mathrm{bg}$',
        'X_axis_bgR': r'Internuclear distance $R\,(a_0)$',


        # LCP Fitting
        'LCP_base_path': r"./DATA/DSState",
        
        
        'color_data_LCP': 'black',
        'color_fit_LCP': 'red',
        'LCP_base_path': r"./DATA/DSState",
        'LCP_sparse_path': r"./DATA/DSStateSparce",

        'p0_guess_E_res': [0.5, -4.0, 0.5, 3.0, 20.0, 1.0],
        'p0_guess_Delta': [1.3, 0.8, 1.5],
        'p0_guess_Gamma': [10.0, 1.0, 1.5, 0.01],
    },


    '2DModel1': {
        'Fit_R0': 2.0,                              # Referencial for fitting, i.e. S(R) = Gamma(E,R)/Gamma(E,R0)
        'base_path': r"./DATA/DSState",
        'output_path': r"./Graphs/Fitting",
        'color_data': 'black',
        'color_fit': 'red',
        'aspect_ratio': 'square',


        # R-scyling factor fitting
        'p0_guess_g': [1.0, 2.0, 1.0, 5.0, 0.5, 1.0, 5.0, 0.5],      
        'X_axis_g': r'Internuclear distance $R\,(a_0)$',
        'Y_axis_g': r'Scaling factor $g^2(R)$',
 

        # Defect fitting
        'defect_path': r"./DATA/RydbergStates",
        'X_axis_defect': r'Energy $(\mathrm{eV})$',
        'Y_axis_defect': r'Background quantum defect $\delta_\mathrm{bg}$',
        'p0_guess_defect': [0.0, 0.0, -2.0],

        # Gamma, Delta energy fitting
        'X_axis_E':  r'Energy $E (\mathrm{eV})$',
        'Y_axis_gamma': r'$\Gamma(R_0, E)\, (\mathrm{eV})$',
        'Y_axis_delta': r'$\Delta(R_0, E)\, (\mathrm{eV})$',
        'p0_guess_gamma': [0.01, 0.005, 5.0, 2.5, 13.0, 4.0, 20.0],
        'y_lim_delta': [-1.0, 1.0],


        # DS energy fitting
        'X_axis_R': r'Internuclear distance $R\,(a_0)$',
        'Y_axis_E': r'Discrete state energy $\epsilon_d\,(\mathrm{eV})$',
        'p0_guess_ds_energy': [-15.0, -2.0, 1.0, 2.0],

        # bg phaseshift fitting
        'p0_guess_bg_phase': [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.5],
        'X_axis_bgE': r'Electron energy $\epsilon\,(\mathrm{eV})$',
        'Y_axis_bgPhase': r'Background phase shift $\delta_\mathrm{bg}$',
        'X_axis_bgR': r'Internuclear distance $R\,(a_0)$',


        # LCP Fitting
        'LCP_base_path': r"./DATA/DSState",
        
        
        'color_data_LCP': 'black',
        'color_fit_LCP': 'red',
        'LCP_base_path': r"./DATA/DSState",
        'LCP_sparse_path': r"./DATA/DSStateSparce",

        'p0_guess_E_res': [0.5, -4.0, 0.5, 3.0, 20.0, 1.0],
        'p0_guess_Delta': [1.3, 0.8, 1.5],
        'p0_guess_Gamma': [10.0, 1.0, 1.5, 0.01],
    }
}