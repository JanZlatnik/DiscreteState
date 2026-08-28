import numpy as np


def harmonic_params_2D(VA):
    a = 0.4
    b = 2.0
    c = 1.5
    d = 0.72
    rc = 5.0
    Rc = 1.5
    r0 = 2.0000399463254763 - 0.23313506121835248*np.tanh(1.6687920618773417*(-1.1726784346053065 + VA)) - 0.47621440751149846*np.tanh(0.699057007565932*(-1.0900928248785635 + VA))
    Omega = 0.4911722132480025 + 0.4788236657617127/np.exp(0.5254434396192675*VA**2) + 0.5824851552152334*np.tanh(0.6948355643578058*(-2.104488212788001 + VA)) + 0.23611621449803974*np.tanh(1.3301337541821736*(-1.7627644298851943 + VA))
    Omega = np.sqrt(Omega)
    V0 = -1/r0 + 1/r0**2 + a*np.exp(-(r0 - rc)**2/b**2) - d*np.exp(-r0**2/4)*np.tanh((VA - Rc)/c)
    return r0, Omega, V0


def harmonic_params_simple(VA):
    r0 = 0.0
    Omega = 1.0
    V0 = -1.5 + VA
    return r0, Omega, V0

def harmonic_turn_off(VA):
    return np.nan, np.nan, np.nan


MODELS = {
    # -------------------------------------------------------------------------------
    '2DModel1' : {
        'potentials': np.linspace(0.7, 10.0, 10),
        'legend_variable_name': 'R',
        'legend_unit': r'\,a_0',
        
        # Potentials
        'Potentials_x_range': [0, 15],
        'Potentials_y_range_V': [-21, 21],
        'Potentials_y_range_scale': 4,
        'Potentials_harmonic_params': harmonic_turn_off,
        'Potentials_Z': 1.0,
        'Potentials_l': 1,

        # PhaseShiftsDefects
        'PhaseShiftsDefects_shifts': [-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2],
        'PhaseShiftsDefects_DS_shifts': [0,0,0,0,0,0,0,0,0,0,0],
        'PhaseShiftsDefects_shifts_H': [0,0,0,0,0,0,0,0,0,0,0],
        'PhaseShiftsDefects_shifts_PHP': [0,0,0,0,0,0,0,0,0,0,0],
        'PhaseShiftsDefects_x_range': [-7.5, 7.5],
        'PhaseShiftsDefects_y_range': [-9, 4],
        'PhaseShiftsDefects_x_rangeAll': [-7.5, 7.5],
        'PhaseShiftsDefects_y_rangeAll': [-9, 4],
        'PhaseShiftsDefects_font_size' : 10,
        'PhaseShiftsDefects_ncol' : 2,

        # LevelShift
        'LevelShift_x_range': [-5, 5],
        'LevelShift_y_range': [-5, 5],

        # LevelShiftsFar
        'LevelShiftFar_base_path': r"./DATA",
        'LevelShiftFar_x_range': [-300, 300],
        'LevelShiftFar_y_range': [-5, 5],

        # GammaExtension
        'GammaExtension_x_range': [-5, 2.5],
        'GammaExtension_y_range': [-0.1, 5.0],

        # Hilbert
        'Hilbert_plot_gamma': False,
        'Hilbert_x_range': [-5, 5],
        'Hilbert_y_range': [-4, 4]

    },

    # -------------------------------------------------------------------------------
    '2DModel2' : {
        'potentials': np.linspace(0.7, 10.0, 10),
        'legend_variable_name': 'R',
        'legend_unit': r'\,a_0',
        
        # Potentials
        'Potentials_x_range': [0, 15],
        'Potentials_y_range_V': [-21, 21],
        'Potentials_y_range_scale': 4,
        'Potentials_harmonic_params': harmonic_turn_off,
        'Potentials_Z': 1.0,
        'Potentials_l': 1,

        # PhaseShiftsDefects
        'PhaseShiftsDefects_shifts': [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
        'PhaseShiftsDefects_DS_shifts': [0,0,0,0,0,0,0,0,0,0,0],
        'PhaseShiftsDefects_shifts_H': [1,1,1,1,1,1,1,1,1,1,1],
        'PhaseShiftsDefects_shifts_PHP': [1,1,1,1,1,1,1,1,1,1,1],
        'PhaseShiftsDefects_x_range': [-7.5, 7.5],
        'PhaseShiftsDefects_y_range': [-6, 4],
        'PhaseShiftsDefects_x_rangeAll': [-7.5, 7.5],
        'PhaseShiftsDefects_y_rangeAll': [-6, 4],
        'PhaseShiftsDefects_font_size' : 10,
        'PhaseShiftsDefects_ncol' : 2,

        # LevelShift
        'LevelShift_x_range': [-2.5, 2.5],
        'LevelShift_y_range': [-3, 3],

        # LevelShiftsFar
        'LevelShiftFar_base_path': r"./DATA",
        'LevelShiftFar_x_range': [-300, 300],
        'LevelShiftFar_y_range': [-5, 5],

        # GammaExtension
        'GammaExtension_x_range': [-5, 2.5],
        'GammaExtension_y_range': [-0.05, 2.0],

        # Hilbert
        'Hilbert_plot_gamma': False,
        'Hilbert_x_range': [-5, 5],
        'Hilbert_y_range': [-4, 4]

    },

    # -------------------------------------------------------------------------------
    'CoulombModel' : {
        'potentials': np.linspace(-0.15, 0.30, 10),
        'legend_variable_name': 'V_A',
        
        # Potentials
        'Potentials_x_range': [0, 7],
        'Potentials_y_range_V': [-50, 50],
        'Potentials_y_range_scale': 4,
        'Potentials_harmonic_params': harmonic_params_simple,
        'Potentials_Z': 1.0,
        'Potentials_l': 0,

        # PhaseShiftsDefects
        'PhaseShiftsDefects_shifts': [0,0,0,0,0,0,0,0,0,0,0],
        'PhaseShiftsDefects_DS_shifts': [1,1,1,1,1,1,1,1,1,1,1],
        'PhaseShiftsDefects_shifts_H': [0,0,0,0,0,0,0,0,0,0,0],
        'PhaseShiftsDefects_shifts_PHP': [0,0,0,0,0,0,0,0,0,0,0],
        'PhaseShiftsDefects_x_range': [-12, 12],
        'PhaseShiftsDefects_y_range': [-4.5, 3.5],
        'PhaseShiftsDefects_x_rangeAll': [-12,12],
        'PhaseShiftsDefects_y_rangeAll': [-6.0, 3.5],
        'PhaseShiftsDefects_font_size' : 10,
        'PhaseShiftsDefects_ncol' : 2,

        # LevelShift
        'LevelShiftFar_base_path': r"./DATAsparse",
        'LevelShift_x_range': [-10, 10],
        'LevelShift_y_range': [-4.1, 4.1],

        # LevelShiftsFar
        'LevelShiftFar_x_range': [-500, 500],
        'LevelShiftFar_y_range': [-5, 5],

        # GammaExtension
        'GammaExtension_x_range': [-5, 2.5],
        'GammaExtension_y_range': [-0.2, 5],

        # Hilbert
        'Hilbert_plot_gamma': False,
        'Hilbert_x_range': [-10, 10],
        'Hilbert_y_range': [-4.1, 4.1]

    },

    # -------------------------------------------------------------------------------
    'NoCoulombModel' : {
        'potentials': np.linspace(-0.15, 0.30, 10),
        'legend_variable_name': 'V_A',
        
        # Potentials
        'Potentials_x_range': [0, 7],
        'Potentials_y_range_V': [-50, 50],
        'Potentials_y_range_scale': 4,
        'Potentials_harmonic_params': harmonic_params_simple,
        'Potentials_Z': 0.0,
        'Potentials_l': 0,

        # PhaseShiftsDefects
        'PhaseShiftsDefects_shifts': [1,1,1,1,1,0,0,0,0,0,0],
        'PhaseShiftsDefects_DS_shifts': [1,1,1,1,1,1,1,1,1,1,1],
        'PhaseShiftsDefects_shifts_H': [0,0,0,0,0,0,0,0,0,0,0],
        'PhaseShiftsDefects_shifts_PHP': [0,0,0,0,0,0,0,0,0,0,0],
        'PhaseShiftsDefects_x_range': [-0.1, 12],
        'PhaseShiftsDefects_y_range': [-3, 3.2],
        'PhaseShiftsDefects_x_rangeAll': [-0.1,12],
        'PhaseShiftsDefects_y_rangeAll': [-3, 3.2],
        'PhaseShiftsDefects_font_size' : 12,
        'PhaseShiftsDefects_ncol' : 1,

        # LevelShift
        'LevelShift_x_range': [-10, 10],
        'LevelShift_y_range': [-4.1, 4.1],

        # LevelShiftsFar
        'LevelShiftFar_base_path': r"./DATAsparse",
        'LevelShiftFar_x_range': [-500, 500],
        'LevelShiftFar_y_range': [-5, 5],

        # GammaExtension
        'GammaExtension_x_range': [-5, 2.5],
        'GammaExtension_y_range': [-0.2, 5],

        # Hilbert
        'Hilbert_plot_gamma': True,
        'Hilbert_x_range': [-10, 10],
        'Hilbert_y_range': [-4.1, 4.1]

    }
}