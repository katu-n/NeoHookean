import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# 1. データの読み込み確認
csv_file = "muscle_comparison_v2.csv"
if not os.path.exists(csv_file):
    print(f"エラー: '{csv_file}' が見つかりません。先にC++プログラムを実行してCSVを生成してください。")
    exit()

df = pd.read_csv(csv_file)

# 2. 理論モデル（解析解）のパラメータ設定 (論文 Table 1 準拠)
p_pass = {'mu': 755.5, 'mu_p': 86.13, 'lambda_0': 0.95, 'c1': 470.0, 'c2': 0.95, 'r': 1.05, 'm': 2500.0, 'a': 2500.0, 'b': 2500.0}
p_act  = {'mu': 755.5, 'mu_p': 86.13, 'lambda_0': 0.95, 'c1': 470.0, 'c2': 0.95, 'r': 1.55, 'm': 30000.0, 'a': 30000.0, 'b': 15000.0}

# Eq. 82: 活性化パラメータαによる補間
def get_params(alpha):
    r = (1.0 - alpha) * p_pass['r'] + alpha * p_act['r']
    m = (1.0 - alpha) * p_pass['m'] + alpha * p_act['m']
    a = (1.0 - alpha) * p_pass['a'] + alpha * p_act['a']
    b = (1.0 - alpha) * p_pass['b'] + alpha * p_act['b']
    return r, m, a, b

# Eq. 72, Eq. 75: 基準エネルギー W0 と 一次負荷応力 sigma_0 の計算
def get_W0_sigma0(lam, alpha):
    mu, mu_p, lambda_0, c1, c2 = p_pass['mu'], p_pass['mu_p'], p_pass['lambda_0'], p_pass['c1'], p_pass['c2']
    
    # Isotropic Matrix 
    I1 = lam**2 + 2.0/lam
    W_iso = (mu/2.0)*(I1 - 3.0)
    sigma_iso = mu * (lam**2 - 1.0/lam)
    
    # Passive Fibers
    W_pass, sigma_pass = 0.0, 0.0
    if lam > 1.0:
        I4 = lam**2
        W_pass = (mu*mu_p/2.0) * (1.0 - alpha) * (I4 - 1.0)**2
        sigma_pass = 2.0 * mu * mu_p * (1.0 - alpha) * (lam**2) * (lam**2 - 1.0)
            
    # Active Fibers
    I4_hat = (lam / lambda_0)**2
    I4_0 = lambda_0**2
    mu_a_val = c1 * np.exp(-(I4_hat - I4_0)/c2)
    
    W_act = (mu * mu_a_val / 2.0) * alpha * (I4_hat - I4_0)**2
    bracket = 1.0 - (I4_hat - I4_0)/(2.0 * c2)
    sigma_act = 2.0 * mu * mu_a_val * alpha * I4_hat * (I4_hat - I4_0) * bracket
    
    W0 = W_iso + W_pass + W_act
    sigma0 = sigma_iso + sigma_pass + sigma_act
    
    return W0, sigma0

# 3. 理論解（解析解）のステップ計算
analytical_stress = []

# 状態遷移用履歴変数
state = "PRIMARY_LOADING"
lambda_prev = 1.0
lambda_max = 1.0
W0_max = 0.0
W0_reversal_reload = 0.0
eta_reversal_reload = 1.0
W0_reversal_unload = 0.0
eta_reversal_unload = 1.0

for i, row in df.iterrows():
    lam = row['Lambda']
    alpha = row['Alpha']
    
    W0_curr, sigma0 = get_W0_sigma0(lam, alpha)
    d_lambda = lam - lambda_prev
    
    # 状態の判定
    if abs(d_lambda) > 1e-9:
        if d_lambda < 0:
            if state == "PRIMARY_LOADING":
                state = "UNLOADING"
            elif state == "RELOADING":
                state = "SECONDARY_UNLOADING"
                W0_reversal_unload, _ = get_W0_sigma0(lambda_prev, alpha)
        else:
            if state in ["UNLOADING", "SECONDARY_UNLOADING"]:
                state = "RELOADING"
                W0_reversal_reload, _ = get_W0_sigma0(lambda_prev, alpha)
                
    lambda_prev = lam
    r, m, a, b = get_params(alpha)
    
    # Active化によりエネルギーが過去の最大を超えた場合のPrimary Loading復帰判定
    if W0_curr > W0_max:
        state = "PRIMARY_LOADING"
        W0_max = W0_curr
        lambda_max = lam
        
    # 各状態に応じた eta の計算
    if state == "PRIMARY_LOADING":
        if lam > lambda_max:
            lambda_max = lam
            W0_max = W0_curr
        eta = 1.0
    elif state == "UNLOADING":
        arg = (W0_max - W0_curr) / m
        eta = 1.0 - (1.0 / r) * np.tanh(arg)
    elif state == "RELOADING":
        arg = (W0_curr - W0_reversal_reload) / a
        term = ((1.0 - eta_reversal_reload) / eta_reversal_reload) * np.tanh(arg)
        eta = eta_reversal_reload * (1.0 + term)
        if eta >= 1.0 or lam >= lambda_max or W0_curr >= W0_max:
            state = "PRIMARY_LOADING"
            eta = 1.0
            W0_max = W0_curr
            lambda_max = lam
    elif state == "SECONDARY_UNLOADING":
        arg = (W0_reversal_unload - W0_curr) / b
        eta = eta_reversal_unload * (1.0 - np.tanh(arg))
    
    # 反転時のeta保存（次のフェーズ移行時に使用）
    if state in ["UNLOADING", "RELOADING", "SECONDARY_UNLOADING"]:
        eta_reversal_unload = eta
        eta_reversal_reload = eta
        
    # Eq. 77: 最終的なCauchy応力
    sigma = eta * sigma0
    analytical_stress.append(sigma / 1e6) # MPa単位

df['Analytical_Stress(MPa)'] = analytical_stress

# 4. グラフのプロット
plt.figure(figsize=(10, 7))

# C++出力結果のプロット（太い青線）
plt.plot(df['Lambda'], df['Stress_Model(MPa)'], '-', label='Paetsch theoretical', linewidth=2, color='blue')

# 理論解のプロット（赤い点：10個飛ばし）
plt.plot(df['Lambda'], df['Analytical_Stress(MPa)'], marker='o', linestyle='none', 
         label='IGA simulation', markersize=5, color='red', markevery=10)

# Neo-Hookeanモデルのプロット（黒色の破線）
plt.plot(df['Lambda'], df['Stress_NH(MPa)'], '--', label='Neo-Hookean Model', linewidth=2, color='black')

plt.xlabel('Stretch $\lambda$', fontsize=16)
plt.ylabel('Cauchy Stress (MPa)', fontsize=16)
plt.legend(fontsize=14)
plt.grid(True, linestyle=':', alpha=0.7)

# 画像として保存＆表示
plt.savefig('comparison_plot.png', dpi=300, transparent=True)
print("グラフを 'comparison_plot.png' として保存しました。")
plt.show()