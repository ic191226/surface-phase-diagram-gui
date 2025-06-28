import numpy as np
import pandas as pd
from typing import Union, List, Any, Tuple, Dict

def make_T_array(T_min: float, T_max: float, T_scale: float
                 ) -> Tuple[np.ndarray, np.ndarray]:
    T_C = np.arange(T_min, T_max + T_scale, T_scale, dtype=float)
    T_K = T_C + 273.15
    return T_C, T_K

def make_partial_pressures(p_max: float, p_scale: float
                           ) -> Tuple[np.ndarray, np.ndarray]:
    pH = np.arange(p_max, -p_scale, -p_scale, dtype=float)
    pN = p_max - pH

    ε = 1e-4
    pH[0],   pN[0]   = p_max - ε, ε
    pH[-1],  pN[-1]  = ε, p_max - ε
    return pH, pN

def compute_mu(
    T_K:         np.ndarray,
    pH:          Union[float, np.ndarray],
    pN:          Union[float, np.ndarray],
    const:       Any,
    graph_type:  str,
    r_params:    List[float]
) -> Dict[float, Dict[str, Dict[float, float]]]:
    """
    H2, N2, III(Ga) の化学ポテンシャル μ をまとめて計算する。
    graph_type=="F" なら pH ごとに、"V/III" なら比率ごとに Ga 分圧を変える。
    """

    # --- 1. キーと配列化 ---
    T_keys  = [round(float(T), 6) for T in T_K]
    pH_arr  = np.atleast_1d(pH).astype(float)
    pN_arr  = np.atleast_1d(pN).astype(float)
    keys_HN = [round(x, 4) for x in pH_arr]

    # print(f"[DEBUG compute_mu] T_keys ({len(T_keys)}): {T_keys[:3]} ... {T_keys[-3:]}")
    # print(f"[DEBUG compute_mu] pH_arr ({pH_arr.shape}): {pH_arr[:3]} ... {pH_arr[-3:]}")
    # print(f"[DEBUG compute_mu] pN_arr ({pN_arr.shape}): {pN_arr[:3]} ... {pN_arr[-3:]}")
    # print(f"[DEBUG compute_mu] keys_HN: {keys_HN[:3]} ... {keys_HN[-3:]}")

    # --- 2. H/N の ζ を温度依存部のみ計算 ---
    def _precalc_zeta(M, I, v):
        ζt = ((2*np.pi*M*const.kB*T_K)/const.h**2)**1.5
        ζr = (8*np.pi**2*I*const.kB*T_K)/(2*const.h**2)
        ζv = 1.0/(1.0 - np.exp(-const.h*v/(const.kB*T_K)))
        return ζt, ζr, ζv

    ζt_H, ζr_H, ζv_H = _precalc_zeta(const.M_H2_KG, const.I_H2_KG_M2, const.v_H2_s)
    ζt_N, ζr_N, ζv_N = _precalc_zeta(const.M_N2_KG, const.I_N2_KG_M2, const.v_N2_s)

    # print(f"[DEBUG compute_mu] ζt_H shape: {ζt_H.shape}, ζr_H shape: {ζr_H.shape}, ζv_H shape: {ζv_H.shape}")
    # print(f"[DEBUG compute_mu] ζt_N shape: {ζt_N.shape}, ζr_N shape: {ζr_N.shape}, ζv_N shape: {ζv_N.shape}")

    # --- 3. H/N の μ をベクトル演算で一括計算 ---
    def _calc_mu(z_t, z_r, z_v, g, P_arr):
        arg   = (g*const.kB*T_K[:,None]/P_arr) * (z_t[:,None]*z_r[:,None]*z_v[:,None])
        muJ   = -const.kB*T_K[:,None]*np.log(arg)
        mu_eV = muJ * const.J2eV
        return muJ, mu_eV

    P_H           = (pH_arr * const.ATM2PA)[None, :]
    P_N           = (pN_arr * const.ATM2PA)[None, :]
    muJ_H, mu_eV_H = _calc_mu(ζt_H, ζr_H, ζv_H, const.g_H, P_H)
    muJ_N, mu_eV_N = _calc_mu(ζt_N, ζr_N, ζv_N, const.g_N, P_N)

    # print(f"[DEBUG compute_mu] P_H shape: {P_H.shape}, P_N shape: {P_N.shape}")
    # print(f"[DEBUG compute_mu] muJ_H shape: {muJ_H.shape}, mu_eV_H shape: {mu_eV_H.shape}")
    # print(f"[DEBUG compute_mu] muJ_N shape: {muJ_N.shape}, mu_eV_N shape: {mu_eV_N.shape}")

    # --- 4. III の Ga 分圧配列とキーをモード別に用意 ---
    if graph_type == "F":
        keys_III = keys_HN
        P_Ga_arr = np.full_like(pH_arr, const.p_Ga_Pa)
        # print(f"[DEBUG compute_mu] graph_type='F', keys_III = keys_HN")
    else:
        start, stop, step = r_params
        ratio_arr         = np.arange(start, stop + step, step, dtype=float)
        keys_III          = [round(r, 4) for r in ratio_arr]
        P_Ga_arr          = (const.p_NH3_atm / ratio_arr) * const.ATM2PA
    #     print(f"[DEBUG compute_mu] graph_type='V/III', ratio_arr ({len(ratio_arr)}): {ratio_arr[:3]} ... {ratio_arr[-3:]}")

    # print(f"[DEBUG compute_mu] keys_III ({len(keys_III)}): {keys_III[:3]} ... {keys_III[-3:]}")

    P_Ga = P_Ga_arr[None, :]
    # print(f"[DEBUG compute_mu] P_Ga shape: {P_Ga.shape}, P_Ga_arr: {P_Ga_arr[:3]} ... {P_Ga_arr[-3:]}")

    # --- 5. III の μ を一括計算 ---
    prefac      = (2*np.pi*const.M_Ga_KG*const.kB*T_K)**1.5
    num_T       = const.kB*T_K * const.g_Ga * prefac
    den         = const.h**3 * P_Ga
    muJ_III     = -const.kB*T_K[:,None] * np.log(num_T[:,None] / den)
    mu_eV_III   = muJ_III * const.J2eV

    # print(f"[DEBUG compute_mu] muJ_III shape: {muJ_III.shape}, mu_eV_III shape: {mu_eV_III.shape}")

    # --- 6. ネスト辞書を組み立てる ---
    # 6-1. H/N の dict を作成
    mu_H_dict = {
        keys_HN[j]: {
            T_keys[i]: {
                'ζtrans': float(ζt_H[i]),
                'ζrot':   float(ζr_H[i]),
                'ζvibr':  float(ζv_H[i]),
                'mu_J':   float(muJ_H[i,j]),
                'mu_eV':  float(mu_eV_H[i,j]),
            }
            for i in range(len(T_K))
        }
        for j in range(len(pH_arr))
    }
    mu_N_dict = {
        keys_HN[j]: {
            T_keys[i]: {
                'ζtrans': float(ζt_N[i]),
                'ζrot':   float(ζr_N[i]),
                'ζvibr':  float(ζv_N[i]),
                'mu_J':   float(muJ_N[i,j]),
                'mu_eV':  float(mu_eV_N[i,j]),
            }
            for i in range(len(T_K))
        }
        for j in range(len(pN_arr))
    }

    # 6-2. III の dict
    mu_III_dict = {
        keys_III[j]: {
            T_keys[i]: {
                'mu_J':  float(muJ_III[i,j]),
                'mu_eV': float(mu_eV_III[i,j]),
            }
            for i in range(len(T_K))
        }
        for j in range(len(keys_III))
    }

    # 6-3. モード別に最終 result を組み立て
    result: Dict[float, Dict[str, Dict[float, float]]] = {}
    if graph_type == "F":
        for key in keys_HN:
            result[key] = {
                'H':   mu_H_dict[key],
                'N':   mu_N_dict[key],
                'III': mu_III_dict[key],
            }
    else:
        for key in keys_III:
            result[key] = {
                'H':   mu_H_dict,
                'N':   mu_N_dict,
                'III': mu_III_dict[key],
            }

    # print(f"[DEBUG compute_mu] Returning result keys: {list(result.keys())[:3]} ... {list(result.keys())[-3:]}")
    return result

def compute_adsorption_energy(energies_ry: Dict[str, Dict[str, float]],
                              additions: Dict[str, Dict[str, int]],
                              const: Any) -> Dict[str, float]:
    """
    全エネルギー（Ry単位）から吸着エネルギー E_ad (eV) を計算して返す。
    """
    E_target = energies_ry['target']
    E_atoms  = energies_ry['atoms']

    # 各原子種の eV 変換済みエネルギー
    H2_ev  = E_atoms['H2']  * const.Ry2eV * 0.5
    N2_ev  = E_atoms['N2']  * const.Ry2eV * 0.5
    III_ev = E_atoms['III'] * const.Ry2eV

    # 理想面の基準エネルギー（eV）
    E0_ideal_eV = E_target['ideal'] * const.Ry2eV

    E_ad: Dict[str, float] = {}
    for struct, E0_ry in E_target.items():
        E0_eV = E0_ry * const.Ry2eV
        add   = additions[struct]
        E_atoms_eV = (
            add['H']   * H2_ev +
            add['N']   * N2_ev +
            add['III'] * III_ev
        )
        E_ad_val = E0_eV - (E0_ideal_eV + E_atoms_eV)
        E_ad[struct] = float(np.round(E_ad_val, 12))

    # print(f"[DEBUG compute_adsorption_energy] structures: {list(E_target.keys())}")
    # print(f"[DEBUG compute_adsorption_energy] E_ad: {E_ad}")
    return E_ad

def compute_gamma(T_K: np.ndarray,
                  pH:   np.ndarray,
                  pN:   np.ndarray,
                  const: Any,
                  energies_ry: Dict[str, Dict[str, float]],
                  additions: Dict[str, Dict[str, int]],
                  mu_dict: Dict[float, Dict[str, Any]],
                  graph_type: str,
                  r_params: List[int]
                 ) -> Dict[float, Dict[str, Dict[float, float]]]:
    """
    表面形成エネルギー γ を pH（Fモード）または V/III 比（V/IIIモード）ごとにまとめて返す。
    """

    # 1) 吸着エネルギーを取得
    E_ad = compute_adsorption_energy(energies_ry, additions, const)

    # デバッグ: 入力キー・配列の情報
    NT = T_K.size
    # print(f"[DEBUG compute_gamma] T_K length: {NT}, range: [{T_K[0]}, {T_K[-1]}]")
    # print(f"[DEBUG compute_gamma] pH shape: {pH.shape}, range: [{pH[0]}, {pH[-1]}]")
    # print(f"[DEBUG compute_gamma] pN shape: {pN.shape}, range: [{pN[0]}, {pN[-1]}]")

    structures = list(energies_ry['target'].keys())
    gamma_dict: Dict[float, Dict[str, Dict[float, float]]] = {}

    # 共通して使う E_ad 配列化
    Ead_vec = np.array([E_ad[s] for s in structures])  # shape (S,)

    if graph_type == "F":
        P = pH.size
        # print(f"[DEBUG compute_gamma] graph_type='F', pH count: {P}")

        # μ 行列をまとめて取り出し → shape (P, NT)
        mu_H_mat   = np.empty((P, NT))
        mu_N_mat   = np.empty((P, NT))
        mu_III_mat = np.empty((P, NT))
        for i, pH_val in enumerate(pH):
            key = round(float(pH_val), 4)
            mu_H_mat[i]   = [mu_dict[key]['H'][round(float(T),6)]['mu_eV']   for T in T_K]
            mu_N_mat[i]   = [mu_dict[key]['N'][round(float(T),6)]['mu_eV']   for T in T_K]
            mu_III_mat[i] = [mu_dict[key]['III'][round(float(T),6)]['mu_eV'] for T in T_K]

        # additions を配列化 (S,3) の形で [ N/2, H/2, III ]
        adds = np.array([[ additions[s]['N']/2,
                           additions[s]['H']/2,
                           additions[s]['III'] ]
                         for s in structures])  # (S,3)

        # ブロードキャストで γ[k,t,s] を計算: shape (P,NT,S)
        # Ead_vec: (S,) → (1,1,S)
        # mu_N_mat etc: (P,NT) → (P,NT,1)
        # adds: (S,3) →添え字で掛ける
        # まとめて引く
        gamma_mat = (
            Ead_vec[None,None,:]
            - mu_N_mat[:,:,None] * adds[None,:,0]
            - mu_H_mat[:,:,None] * adds[None,:,1]
            - mu_III_mat[:,:,None] * adds[None,:,2]
        )
        # 丸め
        gamma_mat = np.round(gamma_mat, 12)

        # 辞書化（Python 側ではこのループは残りますが，数値演算はすべて上で終わり）
        for i, pH_val in enumerate(pH):
            pH_key = round(float(pH_val), 4)
            gamma_dict[pH_key] = {}
            for si, s in enumerate(structures):
                # T_K ごとの小辞書
                gamma_dict[pH_key][s] = {
                    round(float(T_K[t]),6): float(gamma_mat[i,t,si])
                    for t in range(NT)
                }

        # print(f"[DEBUG compute_gamma] Completed F-mode, keys: {list(gamma_dict.keys())[:5]} ...")

    else:
        # V/III モード
        start, stop, step = r_params
        ratio_arr = np.arange(start, stop + 1, step, dtype=float)
        K = ratio_arr.size
        keys_III = [round(r,4) for r in ratio_arr]
        # print(f"[DEBUG compute_gamma] graph_type='V/III', ratio count: {K}, "
        #       f"range: [{keys_III[0]}, {keys_III[-1]}]")

        # H/N は最初の pH 値で固定
        pH_key = round(float(pH[0]), 4)

        # μ 行列をまとめて取り出し → shape (K, NT)
        mu_H_mat   = np.empty((K, NT))
        mu_N_mat   = np.empty((K, NT))
        mu_III_mat = np.empty((K, NT))
        for i, rk in enumerate(keys_III):
            mu_H_mat[i]   = [mu_dict[rk]['H'][pH_key][round(float(T),6)]['mu_eV']   for T in T_K]
            mu_N_mat[i]   = [mu_dict[rk]['N'][pH_key][round(float(T),6)]['mu_eV']   for T in T_K]
            mu_III_mat[i] = [mu_dict[rk]['III'][round(float(T),6)]['mu_eV']         for T in T_K]

        # additions 配列化 (S,3)
        adds = np.array([[ additions[s]['N']/2,
                           additions[s]['H']/2,
                           additions[s]['III'] ]
                         for s in structures])  # (S,3)

        # ブロードキャストで γ[k,t,s] 計算 → shape (K,NT,S)
        gamma_mat = (
            Ead_vec[None,None,:]
            - mu_N_mat[:,:,None] * adds[None,:,0]
            - mu_H_mat[:,:,None] * adds[None,:,1]
            - mu_III_mat[:,:,None] * adds[None,:,2]
        )
        gamma_mat = np.round(gamma_mat, 12)

        # 辞書化
        for i, rk in enumerate(keys_III):
            gamma_dict[rk] = {}
            for si, s in enumerate(structures):
                gamma_dict[rk][s] = {
                    round(float(T_K[t]),6): float(gamma_mat[i,t,si])
                    for t in range(NT)
                }

        # print(f"[DEBUG compute_gamma] Completed V/III-mode, keys: "
        #       f"{keys_III[:3]} ... {keys_III[-3:]}")

    return gamma_dict

def analyze_stability(
    gamma_dict: Dict[float, Dict[str, Dict[float, float]]],
    pH:         np.ndarray,
    pN:         np.ndarray,
    T_C:        np.ndarray,
    graph_type: str,
    r_params:  List[int]
) -> pd.DataFrame:

    # 1) まず Kelvin array を再現
    T_K = T_C + 273.15
    T_keys = [round(float(tk), 6) for tk in T_K]  # これが gamma_dict のキー

    # 2) pH or Ratio keys
    if graph_type == "F":
        p_keys = [round(float(p), 4) for p in pH]
    else:
        start, stop, step = r_params
        p_keys = list(range(start, stop+step, step))

    # 3) 構造リスト
    structures = list(next(iter(gamma_dict.values())).keys())
    S = len(structures)
    nT = len(T_keys)
    nP = len(p_keys)

    # 4) γ を (nP, S, nT) の NumPy 配列にロード
    gamma_arr = np.zeros((nP, S, nT), dtype=float)
    for pi, p_key in enumerate(p_keys):
        for si, struct in enumerate(structures):
            # Kelvinキーで引く
            row = gamma_dict[p_key][struct]
            gamma_arr[pi, si, :] = [row[tk] for tk in T_keys]

    # 5) 各 p_key での安定相＋遷移温度検出
    rows = []
    for pi, p_key in enumerate(p_keys):
        g_arr = gamma_arr[pi]  # shape (S,nT)

        prev = None
        stable_idxs = []
        trans_C = []  # record in Celsius

        for ti, Tc in enumerate(T_C):
            col = g_arr[:, ti]
            mn = col.min()
            cands = np.where(col == mn)[0]

            if ti == 0:
                curr = cands[0]
            else:
                if prev in cands:
                    curr = prev
                elif cands.size == 1:
                    curr = cands[0]
                else:
                    curr = cands[0]

            if ti == 0 or curr != prev:
                stable_idxs.append(curr)
                if ti > 0:
                    trans_C.append(round(float(Tc), 6))
            prev = curr

        # 6) 1行分の辞書組立
        row: Dict[str, Union[float,str]] = {}
        if graph_type == "F":
            ph = float(p_keys[pi])
            pn = float(round(pN[pi], 4))
            row.update({
                'pH': round(ph, 3),
                'pN': round(pn, 3),
                'F' : round(ph/(ph+pn), 3)
            })
        else:
            row['Ratio'] = p_key

        for j, si in enumerate(stable_idxs, start=1):
            row[f'Stable{j}'] = structures[si]
        for k, t in enumerate(trans_C, start=1):
            row[f'Transition{k}'] = t

        rows.append(row)

    # 7) DataFrame 化＆列順
    df = pd.DataFrame(rows)
    stable_cols = sorted([c for c in df if c.startswith("Stable")],
                         key=lambda x: int(x.replace("Stable","")))
    trans_cols  = sorted([c for c in df if c.startswith("Transition")],
                         key=lambda x: int(x.replace("Transition","")))
    if graph_type == "F":
        df = df[['pH','pN','F'] + stable_cols + trans_cols]
    else:
        df = df[['Ratio'] + stable_cols + trans_cols]

    return df
