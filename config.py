import numpy as np

# ────────── デフォルト設定 ──────────
DEFAULT_TMIN    = 700       # 最低温度 [℃]
DEFAULT_TMAX    = 1100      # 最高温度 [℃]
DEFAULT_T_SCALE = 1         # 温度刻み [℃]
DEFAULT_P_MAX   = 0.8       # 総ガス圧 [atm] (pH2 + pN2)
DEFAULT_P_SCALE = 0.001     # 分圧刻み [atm]
DEFAULT_R_MIN   = 100
DEFAULT_R_MAX   = 10000
DEFAULT_R_SCALE = 1
# ────────────────────────────────────

# 全エネルギー定義（Ry単位）
ENERGIES_RY = {
    'target' : {
        'ideal'       : -490.5342757675,
        'Ga_ad'       : -495.0099139979,
        # 'Ga-H'        : ,
        # '2Ga-H'       : ,
        '3Ga-H'       : -494.2674037458,
        # '4Ga-H'       : ,
        # 'Ga-NH'       : ,
        # '2Ga-NH'      : ,
        # 'Ga-NH2'      : ,
        # '2Ga-NH2'     : ,
        # '3Ga-NH2'     : ,
        # '4Ga-NH2'     : ,
        # 'N_ad'        : ,
        # 'N_ad-H'      : ,
        # 'Ga-N+Ga_H'   : ,
        # 'Ga-N+2Ga_H'  : ,
        # 'N_ad-H+Ga-H' : ,
        # 'N_ad-H+Ga-N' : ,
        # 'N_ad+Ga-H'   : ,
        # 'N_br+H2'     : ,
        'monolayer'   : -508.2698749156,
        'bilayer'     : -525.9969245889
    },
    'atoms' : {
        'H2':     -2.33270092,
        'N2':    -39.82813901,
        'III':     -4.21494534
    }
}

# ENERGIES_RY = {
#     'target' : {
#         'VN' : {
#             'ideal':  -685.81114677,
#             '3Ga-H':  -689.55662434,
#             'Ga_ad':  -690.30287423,
#             'bilayer':-721.19582686
#         },
#         'VIII' : {
#             'ideal':  -685.81114677,
#             '3Ga-H':  -689.55662434,
#             'Ga_ad':  -690.30287423,
#             'bilayer':-721.19582686
#         }
#     },
#     'atoms' : {
#         'H2':     -2.33270092,
#         'N2':    -39.82813901,
#         'III':     -4.21494534
#     }
# }

Additions = {
    'ideal'         : {'H':0, 'N':0, 'III':0},
    'Ga_ad'         : {'H':0, 'N':0, 'III':1},
    'Ga-H'          : {'H':1, 'N':0, 'III':0},
    '2Ga-H'         : {'H':2, 'N':0, 'III':0},
    '3Ga-H'         : {'H':3, 'N':0, 'III':0},
    '4Ga-H'         : {'H':3, 'N':0, 'III':0},
    'Ga-NH'         : {'H':0, 'N':1, 'III':0}, # 構造を要確認
    '2Ga-NH'        : {'H':2, 'N':2, 'III':0},
    'Ga-NH2'        : {'H':2, 'N':1, 'III':0},
    '2Ga-NH2'       : {'H':4, 'N':2, 'III':0},
    '3Ga-NH2'       : {'H':6, 'N':3, 'III':0},
    '4Ga-NH2'       : {'H':8, 'N':4, 'III':0},
    'N_ad'          : {'H':0, 'N':1, 'III':0},
    'N_ad-H'        : {'H':1, 'N':1, 'III':0},
    'Ga-N+Ga_H'     : {'H':1, 'N':1, 'III':0},
    'Ga-N+2Ga_H'    : {'H':2, 'N':1, 'III':0},
    'N_ad-H+Ga-H'   : {'H':2, 'N':1, 'III':0},
    'N_ad-H+Ga-N'   : {'H':1, 'N':2, 'III':0},
    'N_ad+Ga-H'     : {'H':1, 'N':1, 'III':0},
    'N_br+H2'       : {'H':2, 'N':1, 'III':0},
    'monolayer'     : {'H':0, 'N':0, 'III':4},
    'bilayer'       : {'H':0, 'N':0, 'III':8}
}

class Constants:
    """
    物理定数および変換係数をまとめたクラス
    """

    # プランク定数
    h = 6.62607015e-34         # [J·s]
    # ボルツマン定数
    kB = 1.38064852e-23        # [J/K]
    # 円周率
    PI = np.pi
    # 統一原子質量単位
    AMU = 1.66053904e-27       # [kg]
    # ボーア半径
    bohr = 0.529177210903e-10  # [m]
    # Joule→eV 変換係数
    J2eV = 6.241509074e18      # [eV/J]
    # atm→Pa 変換係数
    ATM2PA = 101325            # [Pa/atm]
    # Ry→eV 変換係数
    Ry2eV = 13.605662285137    # [eV/Ry]
    # 波数→振動数変換
    cm2Hz = 29979245800        # [Hz/(cm⁻¹)]

    # 分子の縮退度
    g_H = 1
    g_N = 1
    g_Ga = 2

    # 質量 (amu)
    M_H2_AMU = 2.0158
    M_N2_AMU = 28.0134
    M_Ga_AMU = 69.72

    # 質量 (kg)
    M_H2_KG = M_H2_AMU * AMU
    M_N2_KG = M_N2_AMU * AMU
    M_Ga_KG = M_Ga_AMU * AMU

    # 慣性モーメント (amu·Å²)
    I_H2_AMU_A2 = 1.0099
    I_N2_AMU_A2 = 30.72982

    # 慣性モーメント (kg·m²)
    I_H2_KG_M2 = I_H2_AMU_A2 * AMU * bohr**2
    I_N2_KG_M2 = I_N2_AMU_A2 * AMU * bohr**2

    # 振動数 (cm⁻¹ → Hz)
    v_H2_cm = 4381.0
    v_N2_cm = 2378.1
    v_H2_s  = v_H2_cm * cm2Hz
    v_N2_s  = v_N2_cm * cm2Hz

    # 分圧
    p_Ga_atm = 1e-4           # [atm]
    p_Ga_Pa  = p_Ga_atm * ATM2PA
    p_NH3_atm = 0.2         # [atm]
    p_NH3_Pa  = p_NH3_atm * ATM2PA
