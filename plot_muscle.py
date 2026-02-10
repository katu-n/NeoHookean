import pandas as pd
import matplotlib.pyplot as plt

# CSVファイルの読み込み
filename = "muscle_transition_test.csv"
try:
    data = pd.read_csv(filename)
except FileNotFoundError:
    print(f"Error: {filename} not found. Please run the simulation first.")
    exit()

# データの取得
lambda_val = data["Lambda"]
sigma_pass = data["Sigma_Passive(Pa)"] / 1000.0 # kPaに変換
sigma_act = data["Sigma_Active(Pa)"] / 1000.0   # kPaに変換

# プロットの設定
plt.figure(figsize=(10, 6))

# Passive状態のプロット (破線)
plt.plot(lambda_val, sigma_pass, label="Passive State", linestyle="--", color="blue", linewidth=2)

# Active状態のプロット (実線)
plt.plot(lambda_val, sigma_act, label="Active State", linestyle="-", color="red", linewidth=2)

# グラフの装飾
plt.title("Muscle Stress-Stretch Response (Uniaxial)", fontsize=16)
plt.xlabel("Stretch Ratio $\lambda$", fontsize=14)
plt.ylabel("Cauchy Stress $\sigma$ (kPa)", fontsize=14)
plt.legend(fontsize=12)
plt.grid(True, linestyle=':', alpha=0.6)

# 保存と表示
plt.savefig("muscle_transition.png", dpi=300)
print("Graph saved as muscle_response_plot.png")
# GUI環境がある場合のみ表示
# plt.show()

