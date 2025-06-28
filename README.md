# 汎用表面相図 GUI

汎用的な表面相図（Surface Phase Diagram）を簡単に描画・解析できるデスクトップGUIアプリケーションです。

## 🎯 概要

* **温度（℃）** と **分圧（H₂/N₂比またはV/III比）** に基づいて、

  * 化学ポテンシャルの計算
  * 表面形成エネルギーの算出
  * 安定構造の解析
  * 相図の生成とテーブル表示
* 凡例カラーのカスタマイズや、グラフ・表の保存機能を備えます。

## 🚀 特徴

* バックグラウンドスレッドで重い計算を実行しつつ、進捗状況をモーダルウィンドウで表示
* `tksheet` を用いたインタラクティブなテーブル表示
* `matplotlib` による高品質な相図描画
* カラーピッカーでフェーズごとの色を変更
* 図（PNG）・表（CSV/Excel）を簡単に保存

## 🛠️ 前提条件

* Python 3.7 以上
* 必要なパッケージ:

  * 標準ライブラリ: `tkinter`, `threading`, `re`, `time`
  * 外部ライブラリ:

    * `numpy`
    * `pandas`
    * `matplotlib`
    * `tksheet`

## 📦 インストール

```bash
pip install numpy pandas matplotlib tksheet
```

## ⚙️ 設定

プロジェクトルートの `config.py` で各種パラメータを調整できます。

* **温度設定**

  * `DEFAULT_TMIN` : 最低温度 \[℃]
  * `DEFAULT_TMAX` : 最高温度 \[℃]
  * `DEFAULT_T_SCALE` : 温度刻み \[℃]
* **分圧設定**

  * `DEFAULT_P_MAX` : 総ガス圧 \[atm]
  * `DEFAULT_P_SCALE` : 分圧刻み \[atm]
* **V/III 比**

  * `DEFAULT_R_MIN`, `DEFAULT_R_MAX`, `DEFAULT_R_SCALE`
* **全エネルギー（Ry単位）**

  * `ENERGIES_RY['target']`, `ENERGIES_RY['atoms']`
* **吸着構造ごとの原子数**

  * `Additions`
* **物理定数**

  * `Constants` クラスにまとめて定義

## 🚀 使い方

```bash
python GUI.py
```

1. アプリ起動後、ウィンドウ右上で **Graph Type** を選択
2. **Tmin/Tmax** や **X Tick/Y Tick** を必要に応じて変更
3. **更新** ボタンをクリックして再描画
4. **図保存** / **表保存** ボタンから出力形式（PNG, CSV, XLSX）を選択

## 📁 ファイル構成

```
├── config.py       # 設定値定義
├── utils.py        # 計算ロジック実装
├── GUI.py          # TkinterベースのGUIアプリ本体
└── README.md       # 本ドキュメント
```

---

ご不明点やバグ報告などありましたら、Issue をお立てください。
