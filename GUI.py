import tkinter as tk
from tkinter import ttk, messagebox, filedialog, colorchooser
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from tksheet import Sheet

import pandas as pd
import threading
import numpy as np
import re
import time

import config
import utils

class PhaseDiagramApp(tk.Tk):
    """
    汎用表面相図 GUI アプリケーション
    """
    def __init__(self):
        super().__init__()
        self.title("汎用表面相図 GUI")
        self.minsize(900, 600)

        # -- データ・パラメータ初期化 --
        self.stability_df = pd.DataFrame()
        self.gamma_dict = {}
        self.param_entries = {}

        # 温度・圧力スケール
        self.T_min   = config.DEFAULT_TMIN
        self.T_max   = config.DEFAULT_TMAX
        self.T_scale = config.DEFAULT_T_SCALE
        self.p_max   = config.DEFAULT_P_MAX
        self.p_scale = config.DEFAULT_P_SCALE

        # V/III 比パラメータ
        self.R_params = [
            config.DEFAULT_R_MIN,
            config.DEFAULT_R_MAX,
            config.DEFAULT_R_SCALE
        ]

        # グラフ目盛
        self.X_tick = 100
        self.Y_tick = 0.1

        # エネルギー・定数
        self.energys   = config.ENERGIES_RY
        self.additions = config.Additions
        self.const     = config.Constants()

        # Always-on-top フラグ
        self.var_topmost = tk.BooleanVar(value=False)

        # UI 初期化
        self._init_ui()

        # 起動時に一度だけ計算開始
        self.after(100, self._run_with_waitbar)

    def _init_ui(self):
        """UI要素をすべて初期化"""
        paned = tk.PanedWindow(self, orient='horizontal')
        paned.pack(fill='both', expand=True)

        # 左：matplotlib キャンバス
        left = ttk.Frame(paned)
        left.columnconfigure(0, weight=1)
        left.rowconfigure(0, weight=1)
        paned.add(left, minsize=600)

        self.fig = Figure(figsize=(5,5), dpi=100)
        self.ax  = self.fig.add_subplot(111)
        self.ax.set_facecolor('white')
        self.ax.set_box_aspect(1)
        self.canvas = FigureCanvasTkAgg(self.fig, master=left)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew')

        # 右：コントロール＆テーブル
        right = ttk.Frame(paned)
        right.columnconfigure((0,1), weight=1)
        right.rowconfigure(2, weight=1)
        paned.add(right, minsize=300)

        # ─ コントロール部 ─
        ctrl = ttk.Frame(right, padding=5)
        ctrl.grid(row=0, column=0, columnspan=2, sticky='we')
        ctrl.columnconfigure(0, weight=1)
        ctrl.columnconfigure(1, weight=1)

        # パラメータ枠
        prm = ttk.Labelframe(ctrl, text='パラメータ', padding=5)
        prm.grid(row=0, column=0, sticky='nswe', padx=5, pady=2)
        prm.columnconfigure(1, weight=1)

        # Graph Type
        ttk.Label(prm, text='Graph Type:').grid(row=0, column=0, sticky='e', padx=4, pady=4)
        self.var_graph_type = tk.StringVar(value='F')
        self.comb = ttk.Combobox(
            prm,
            textvariable=self.var_graph_type,
            values=['F', 'V/III Nrich', 'V/III Hrich'],
            state='readonly'
        )
        self.comb.grid(row=0, column=1, sticky='we', padx=4, pady=4)
        self.comb.bind('<<ComboboxSelected>>', self._on_graph_type_change)

        # その他パラメータ入力
        self.param_vars = {}
        defaults = {
            'Title': '',
            'Tmin': self.T_min,
            'Tmax': self.T_max,
            'X Tick': self.X_tick,
            'Y Tick': self.Y_tick,
        }
        for i, key in enumerate(['Title','Tmin','Tmax','X Tick','Y Tick'], start=1):
            var = tk.DoubleVar(value=defaults[key]) if isinstance(defaults[key], (int,float)) else tk.StringVar(value=defaults[key])
            self.param_vars[key] = var
            ttk.Label(prm, text=f'{key}:').grid(row=i, column=0, sticky='e', padx=4, pady=2)
            entry = ttk.Entry(prm, textvariable=var)
            entry.grid(row=i, column=1, sticky='we', padx=4, pady=2)
            self.param_entries[key] = entry

        # 凡例色選択枠
        leg = ttk.Labelframe(ctrl, text='Legend Colors', padding=5)
        leg.grid(row=0, column=1, sticky='nswe', padx=5, pady=2)
        self.leg = leg

        self.leg_canvas = tk.Canvas(leg, borderwidth=0, highlightthickness=0, height=200, width=100)
        self.leg_scrollbar = ttk.Scrollbar(leg, orient="vertical", command=self.leg_canvas.yview)
        self.leg_scrollable_frame = ttk.Frame(self.leg_canvas)
        self.leg_scrollable_frame.bind(
            "<Configure>",
            lambda e: self.leg_canvas.configure(scrollregion=self.leg_canvas.bbox("all"))
        )
        self.leg_canvas.create_window((0,0), window=self.leg_scrollable_frame, anchor="nw")
        self.leg_canvas.configure(yscrollcommand=self.leg_scrollbar.set)
        self.leg_canvas.pack(side="left", fill="both", expand=True)
        self.leg_scrollbar.pack(side="right", fill="y")

        self.phase_names = []
        self.phase_colors = {}

        # ボタン群
        btn_fr = ttk.Frame(right, padding=5)
        btn_fr.grid(row=1, column=0, columnspan=2, sticky='we')
        btn_fr.columnconfigure((0,1,2,3), weight=1)
        ttk.Button(btn_fr, text='更新', command=self._on_redraw).grid(row=0, column=0, sticky='we', padx=3)
        ttk.Button(btn_fr, text='図保存', command=self._save_figure).grid(row=0, column=1, sticky='we', padx=3)
        ttk.Button(btn_fr, text='表保存', command=self._save_table).grid(row=0, column=2, sticky='we', padx=3)
        ttk.Checkbutton(btn_fr, text='常に最前面', variable=self.var_topmost, command=self._update_topmost)\
            .grid(row=0, column=3, sticky='we', padx=3)

        # テーブル表示
        self.sheet = Sheet(right)
        self.sheet.grid(row=2, column=0, columnspan=2, sticky='nsew', padx=5, pady=5)

    def _update_topmost(self):
        """常に最前面モードを切り替え"""
        top = self.var_topmost.get()
        self.attributes('-topmost', top)
        if getattr(self, 'waitbar', None) and self.waitbar.winfo_exists():
            self.waitbar.attributes('-topmost', top)
            if top:
                self.waitbar.lift()

    def _on_graph_type_change(self, event=None):
        """Graph Type 切替時の処理"""
        gt = self.var_graph_type.get()
        state = 'disabled' if gt.startswith('V/III') else 'normal'
        self.param_entries['Y Tick'].config(state=state)
        self._on_redraw()

    def _center_window(self, win, width):
        """ウィンドウwinをメインの中央に、横幅widthで配置"""
        win.update_idletasks()
        mh = win.winfo_height()
        mx, my = self.winfo_rootx(), self.winfo_rooty()
        mw, mh_main = self.winfo_width(), self.winfo_height()
        x = mx + (mw - width)//2
        y = my + (mh_main - mh)//2
        win.geometry(f"{width}x{mh}+{x}+{y}")

    def _show_quick_waitbar(self, title: str, message: str, callback):
        """
        計算を伴わない UI 更新時に、簡易的な待ちバーを表示して
        callback() を呼び出します。
        """
        wb = tk.Toplevel(self)
        wb.title(title)
        wb.transient(self)
        wb.grab_set()
        lbl = ttk.Label(wb, text=message)
        lbl.pack(padx=10, pady=(10,0))
        pb = ttk.Progressbar(wb, mode='indeterminate', length=200)
        pb.pack(padx=10, pady=10)
        pb.start(50)
        self._center_window(wb, 300)
        
        wb.update_idletasks()
        wb.update()

        def finish():
            try:
                callback()
            finally:
                wb.destroy()

        # メインスレッドのアイドル時に実際の更新を行い、待ちバーを閉じる
        self.after(20, finish)

    def _on_redraw(self):
        """パラメータ変更時の再計算判定・再描画"""
        new_Tmin = self.param_vars['Tmin'].get()
        new_Tmax = self.param_vars['Tmax'].get()
        new_gt   = self.var_graph_type.get()
        need_recalc = (
            not hasattr(self, '_last_params') or
            new_Tmin != self._last_params['Tmin'] or
            new_Tmax != self._last_params['Tmax'] or
            new_gt   != self._last_params['Graph Type']
        )
        self._last_params = {'Tmin': new_Tmin, 'Tmax': new_Tmax, 'Graph Type': new_gt}
        if need_recalc:
            self._run_with_waitbar()
        else:
            def update_ui():
                self.generate_phase_diagram(self.stability_df)
                self.show_table()
            
            self._show_quick_waitbar(
                title="更新中...",
                message="UIを更新しています...",
                callback=update_ui
            )

    def _run_with_waitbar(self):
        """バックグラウンド計算＋モーダル待ちバー"""
        wb = tk.Toplevel(self)
        wb.title("計算中…")
        wb.transient(self)
        wb.grab_set()

        lbl = ttk.Label(wb, text="")
        lbl.pack(padx=10, pady=(10,0))
        pb = ttk.Progressbar(wb, mode='indeterminate')
        pb.pack(fill='x', expand=True, padx=10, pady=10)
        pb.start(10)

        self._update_topmost()
        self._center_window(wb, 400)

        def set_text(msg):
            lbl.config(text=msg)
            lbl.update_idletasks()

        def task():
            try:
                print(f"[INFO] {self.comb.get()} 計算開始")
                set_text("1/5 温度・分圧配列生成中…")
                t_start = time.perf_counter()
                T_C, T_K = utils.make_T_array(self.T_min, self.T_max, self.T_scale)
                if self.comb.get() == "F":
                    pH_, pN_ = utils.make_partial_pressures(self.p_max, self.p_scale)
                    comb_str = "F"
                else:
                    if self.comb.get() == "V/III Hrich":
                        pH_ = np.array([0.7999])
                        pN_ = np.array([0.8 - 0.7999])
                    else:
                        pH_ = np.array([0.8 - 0.7999])
                        pN_ = np.array([0.7999])
                    comb_str = "V/III"
                t_end = time.perf_counter()
                print(f"[TIMING] 1/5 温度・分圧配列生成: {t_end - t_start:.3f} 秒")

                set_text("2/5 化学ポテンシャル計算中…")
                t_start = time.perf_counter()
                mu_dict = utils.compute_mu(
                    T_K=T_K, pH=pH_, pN=pN_,
                    const=self.const,
                    graph_type=comb_str,
                    r_params=self.R_params
                )
                self.mu_dict = mu_dict
                t_end = time.perf_counter()
                print(f"[TIMING] 2/5 compute_mu: {t_end - t_start:.3f} 秒")

                set_text("3/5 表面形成エネルギー計算中…")
                t_start = time.perf_counter()
                gamma_dict = utils.compute_gamma(
                    T_K=T_K, pH=pH_, pN=pN_,
                    const=self.const,
                    energies_ry=self.energys,
                    additions=self.additions,
                    mu_dict=mu_dict,
                    graph_type=comb_str,
                    r_params=self.R_params
                )
                self.gamma_dict = gamma_dict
                t_end = time.perf_counter()
                print(f"[TIMING] 3/5 compute_gamma: {t_end - t_start:.3f} 秒")

                set_text("4/5 安定構造解析中…")
                t_start = time.perf_counter()
                self.stability_df = utils.analyze_stability(
                    gamma_dict=gamma_dict,
                    pH=pH_, pN=pN_,
                    T_C=T_C,
                    graph_type=comb_str,
                    r_params=self.R_params
                )
                t_end = time.perf_counter()
                print(f"[TIMING] 4/5 analyze_stability: {t_end - t_start:.3f} 秒")

                def finish():
                    lbl.config(text="5/5 図と表の更新中…")
                    wb.update_idletasks()
                    t1 = time.perf_counter()
                    self.generate_phase_diagram(self.stability_df)
                    self.show_table()
                    t2 = time.perf_counter()
                    print(f"[TIMING] 5/5 GUI 更新: {t2-t1:.3f} 秒")
                    
                    self._last_params = {
                        'Tmin': self.param_vars['Tmin'].get(),
                        'Tmax': self.param_vars['Tmax'].get(),
                        'Graph Type': self.var_graph_type.get()
                    }
                    wb.destroy()

                self.after(0, finish)
            except Exception as e:
                self.after(0, lambda err=e: messagebox.showerror('Error', str(err)))
                self.after(0, wb.destroy)

        threading.Thread(target=task, daemon=True).start()

    def generate_phase_diagram(self, df):
        """DataFrame から相図を生成表示"""
        if df is None or df.empty:
            self.ax.clear()
            self.canvas.draw()
            return

        if not any(c.startswith('Stable') for c in df.columns):
            messagebox.showerror("Error", "Stable列がありません")
            self.ax.clear()
            self.canvas.draw()
            return

        struct_cols = sorted(
            [c for c in df.columns if c.startswith('Stable')],
            key=lambda x: int(re.search(r"\d+", x).group())
        )
        temp_cols = sorted(
            [c for c in df.columns if c.startswith('Transition')],
            key=lambda x: int(re.search(r"\d+", x).group())
        )

        phase_names = []
        for col in struct_cols:
            for ph in df[col].dropna().unique():
                s = str(ph).strip()
                if s and s not in phase_names:
                    phase_names.append(s)
        self.phase_names = phase_names

        cmap = matplotlib.colormaps['Set2'](np.linspace(0,1,len(phase_names)))
        if not self.phase_colors:
            self.phase_colors = {}
        for i, name in enumerate(phase_names):
            if name not in self.phase_colors:
                self.phase_colors[name] = matplotlib.colors.to_hex(cmap[i])

        graph_type = self.var_graph_type.get()
        graph_lb = "F" if graph_type=="F" else "Ratio"
        if graph_lb not in df.columns:
            messagebox.showerror("Error", f"{graph_lb}列がありません")
            self.ax.clear()
            self.canvas.draw()
            return

        F_vals = df[graph_lb].values
        N = len(F_vals)

        dF = np.diff(F_vals)
        half = np.empty_like(F_vals)
        half[1:-1] = 0.5*(dF[:-1] + dF[1:])
        half[0], half[-1] = dF[0], dF[-1]
        F_edges = np.concatenate([
            [F_vals[0] - half[0]/2],
            (F_vals[:-1] + F_vals[1:]) / 2,
            [F_vals[-1] + half[-1]/2]
        ])

        Tmin = self.param_vars['Tmin'].get()
        Tmax = self.param_vars['Tmax'].get()
        temps = [Tmin] + [float(x) for col in temp_cols for x in df[col].dropna()] + [Tmax]
        T_edges = np.unique(temps)

        region = np.zeros((N, len(T_edges)-1), dtype=int)
        for i, row in df.iterrows():
            phases = row[struct_cols].dropna().tolist()
            bounds = [Tmin] + [float(row[c]) for c in temp_cols if not np.isnan(row[c])] + [Tmax]
            for k, ph in enumerate(phases):
                mask = (T_edges[:-1] >= bounds[k]) & (T_edges[1:] <= bounds[k+1])
                region[i, mask] = phase_names.index(str(ph).strip())

        self.ax.clear()
        cmap2 = ListedColormap([self.phase_colors[n] for n in phase_names])
        X, Y = np.meshgrid(T_edges, F_edges)
        self.ax.pcolormesh(X, Y, region, cmap=cmap2, shading='flat')

        for i in range(N):
            for j in range(region.shape[1]):
                if j < region.shape[1]-1 and region[i,j] != region[i,j+1]:
                    x  = T_edges[j+1]
                    y0 = F_edges[i]
                    y1 = F_edges[i+1]
                    self.ax.plot([x,x],[y0,y1], color='k', linewidth=1)
                if i < N-1 and region[i,j] != region[i+1,j]:
                    y  = F_edges[i+1]
                    x0 = T_edges[j]
                    x1 = T_edges[j+1]
                    self.ax.plot([x0,x1],[y,y], color='k', linewidth=1)

        self.ax.set_xlim(Tmin, Tmax)
        self.ax.xaxis.set_major_locator(MultipleLocator(self.param_vars['X Tick'].get()))
        self.ax.set_xlabel('Temperature [℃]')
        if graph_type=="F":
            self.ax.set_ylim(0,1)
            self.ax.yaxis.set_major_locator(MultipleLocator(self.param_vars['Y Tick'].get()))
            self.ax.set_ylabel('F = pH₂/(pH₂ + pN₂)')
        else:
            self.ax.set_yscale('log')
            self.ax.set_ylim(100,10000)
            self.ax.set_yticks([100,1000,10000])
            self.ax.yaxis.set_major_formatter(ScalarFormatter())
            self.ax.set_ylabel('Ⅴ / Ⅲ Ratio')

        handles = [Patch(facecolor=self.phase_colors[n], edgecolor='k', label=n)
                   for n in phase_names]
        self.ax.legend(handles=handles, loc='upper right', frameon=True)

        title = self.param_vars['Title'].get().strip()
        if title:
            self.ax.set_title(title)

        self.canvas.draw()
        self.update_legend_colors()

    def show_table(self):
        """tksheet で再安定構造テーブルを表示"""
        df = self.stability_df
        cols = df.columns.tolist()
        data = df.values.tolist()
        self.sheet.headers(cols)
        self.sheet.set_sheet_data(data)

    def update_legend_colors(self):
        """凡例色ボタン再生成"""
        for w in self.leg_scrollable_frame.winfo_children():
            w.destroy()
        for i, name in enumerate(self.phase_names):
            color = self.phase_colors.get(name, "#cccccc")
            btn = tk.Button(
                self.leg_scrollable_frame,
                text=name, bg=color, width=12,
                command=lambda n=name: self._choose_color(n)
            )
            btn.grid(row=i, column=0, sticky='we', padx=2, pady=2)

    def _choose_color(self, phase_name):
        """色選択ダイアログ"""
        c = colorchooser.askcolor(title=f"{phase_name}の色を選択")[1]
        if c:
            self.phase_colors[phase_name] = c
            def update_legend():
                self.generate_phase_diagram(self.stability_df)
                self.update_legend_colors()
            self._show_quick_waitbar(
                title="凡例色変更中…",
                message=f"{phase_name} の色を更新しています…",
                callback=update_legend
            )

    def _save_figure(self):
        """図を保存するときにも waitbar を表示"""
        path = filedialog.asksaveasfilename(defaultextension='.png', filetypes=[('PNG','*.png')])
        if not path:
            return
        wb = tk.Toplevel(self)
        wb.title("図を保存中…")
        wb.transient(self)
        wb.grab_set()
        lbl = ttk.Label(wb, text="保存中…")
        lbl.pack(padx=10, pady=(10,0))
        pb = ttk.Progressbar(wb, mode='indeterminate')
        pb.pack(fill='x', expand=True, padx=10, pady=10)
        pb.start(10)
        self._update_topmost()
        self._center_window(wb, 400)

        def do_save():
            self.fig.savefig(path)
            wb.destroy()
            messagebox.showinfo('保存', f'図を保存しました: {path}')

        threading.Thread(target=do_save, daemon=True).start()

    def _save_table(self):
        """
        再安定構造・表面形成エネルギー・化学ポテンシャルのいずれかを
        1, 2, 3 のボタンで選択して、waitbar付きで保存します。
        """
        # --- 1) カスタム選択ダイアログ ---
        choice_win = tk.Toplevel(self)
        choice_win.title("保存する表を選択")
        choice_win.transient(self)
        choice_win.grab_set()
        ttk.Label(choice_win, text="保存する表を選んでください：").pack(padx=20, pady=(10,5))
        sel = {'value': None}
        def select(n):
            sel['value'] = n
            choice_win.destroy()

        btn1 = ttk.Button(choice_win, text="1: 再安定構造", command=lambda: select(1))
        btn2 = ttk.Button(choice_win, text="2: 表面形成エネルギー (Excel)", command=lambda: select(2))
        btn3 = ttk.Button(choice_win, text="3: 化学ポテンシャル (Excel)", command=lambda: select(3))
        for b in (btn1, btn2, btn3):
            b.pack(fill='x', padx=20, pady=5)

        self._center_window(choice_win, 300)
        self.wait_window(choice_win)
        choice = sel['value']
        if choice is None:
            return  # キャンセル

        # --- 2) ファイルダイアログ ---
        if choice == 1:
            filetypes = [('CSV','*.csv'), ('Excel','*.xlsx')]
        else:
            filetypes = [('Excel','*.xlsx')]
        path = filedialog.asksaveasfilename(defaultextension=filetypes[0][1], filetypes=filetypes)
        if not path:
            return

        # --- 3) モーダルな待ちバーウィンドウ ---
        wb = tk.Toplevel(self)
        wb.title("表を保存中…")
        wb.transient(self)
        wb.grab_set()
        lbl = ttk.Label(wb, text="保存中…")
        lbl.pack(padx=10, pady=(10,0))
        pb = ttk.Progressbar(wb, mode='indeterminate')
        pb.pack(fill='x', expand=True, padx=10, pady=10)
        pb.start(10)
        self._update_topmost()
        self._center_window(wb, 400)

        def do_save():
            try:
                # 1 → stability_df
                if choice == 1:
                    if path.lower().endswith('.csv'):
                        self.stability_df.to_csv(path, index=False)
                    else:
                        self.stability_df.to_excel(path, index=False)
                # 2 → gamma_dict
                elif choice == 2:
                    keys = list(self.gamma_dict.keys())
                    total = len(keys)
                    pb.stop()
                    pb.config(mode='determinate', maximum=total, value=0)
                    lbl.config(text=f"保存中… 0/{total}")
                    wb.update_idletasks()
                    with pd.ExcelWriter(path) as writer:
                        for idx, key in enumerate(keys, start=1):
                            pb['value'] = idx
                            lbl.config(text=f"保存中… {idx}/{total}")
                            wb.update_idletasks()
                            gdict = self.gamma_dict[key]
                            first_struct = next(iter(gdict))
                            T_list = list(gdict[first_struct].keys())
                            data = {'T[K]': T_list}
                            for struct, vals in gdict.items():
                                data[struct] = [vals[T] for T in T_list]
                            pd.DataFrame(data).to_excel(writer, sheet_name=str(key), index=False)
                # 3 → mu_dict
                else:
                    keys = list(self.mu_dict.keys())
                    total = len(keys)
                    pb.stop()
                    pb.config(mode='determinate', maximum=total, value=0)
                    lbl.config(text=f"保存中… 0/{total}")
                    wb.update_idletasks()
                    with pd.ExcelWriter(path) as writer:
                        for idx, key in enumerate(keys, start=1):
                            pb['value'] = idx
                            lbl.config(text=f"保存中… {idx}/{total}")
                            wb.update_idletasks()
                            mdict = self.mu_dict[key]
                            T_list = list(mdict['H'].keys())
                            data = {'T[K]': T_list}
                            # H, N は ζtrans, ζrot, ζvibr, mu_J, mu_eV
                            for sp in ['H', 'N']:
                                for metric in ['ζtrans', 'ζrot', 'ζvibr', 'mu_J', 'mu_eV']:
                                    col = f"{sp}_{metric}"
                                    data[col] = [mdict[sp][T][metric] for T in T_list]

                            # III は mu_J, mu_eV のみ
                            for metric in ['mu_J', 'mu_eV']:
                                col = f"III_{metric}"
                                data[col] = [mdict['III'][T][metric] for T in T_list]
                            df_mu = pd.DataFrame(data)
                            df_mu.to_excel(writer, sheet_name=str(key), index=False)
            except Exception as e:
                self.after(0, lambda: messagebox.showerror('保存エラー', str(e)))
            finally:
                wb.destroy()
                messagebox.showinfo('保存', f'表を保存しました: {path}')

        threading.Thread(target=do_save, daemon=True).start()

if __name__ == '__main__':
    app = PhaseDiagramApp()
    app.mainloop()
