import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# CSVファイルの読み込み
filename = "comparison_experiment.csv"
try:
    data = pd.read_csv(filename)
except FileNotFoundError:
    print(f"Error: {filename} not found. Please run the C++ simulation first.")
    exit()

# グラフのセットアップ
fig, ax = plt.subplots(figsize=(10, 7))

# --- プロット 1: Zhao Model (Passive Only) ---
# Zhaoモデルは刺激の影響を受けないため、全データをプロットして基準線とする
# 見やすくするため、少し間引いてプロット (markevery)
ax.plot(data["Lambda"], data["Stress_Zhao(kPa)"], 
        label="Zhao Model (Passive Ref.)", color="blue", linestyle="--", linewidth=2, alpha=0.7)

# --- プロット 2: Integrated Model (Phaseごとに色分け) ---

# Phase 1: Passive Loading (Gray)
phase_pass_load = data[data["Phase_Label"] == "Passive_Load"]
ax.plot(phase_pass_load["Lambda"], phase_pass_load["Stress_Integrated(kPa)"],
        label="Integrated (Passive Load)", color="gray", linewidth=3, alpha=0.8)

# Phase 2: Passive Unloading (Gray Dashed - option, usually overlaps)
# (省略: グラフが混雑するため)

# Phase 3: Activation (Red Vertical Line)
# lambda=1.0 での垂直立ち上がり
phase_act = data[data["Phase_Label"] == "Activation"]
ax.plot(phase_act["Lambda"], phase_act["Stress_Integrated(kPa)"],
        label="Integrated (Activation)", color="red", linewidth=3)

# Phase 4: Active Loading (Dark Red)
phase_act_load = data[data["Phase_Label"] == "Active_Load"]
ax.plot(phase_act_load["Lambda"], phase_act_load["Stress_Integrated(kPa)"],
        label="Integrated (Active Load)", color="darkred", linewidth=3)

# Phase 5: Active Unloading (Pink/Light Red)
phase_act_unload = data[data["Phase_Label"] == "Active_Unload"]
ax.plot(phase_act_unload["Lambda"], phase_act_unload["Stress_Integrated(kPa)"],
        label="Integrated (Active Unload)", color="salmon", linestyle="-.", linewidth=2)

# --- 論文での主張点を強調する矢印とテキスト ---

# 1. 残留応力 (Residual Stress)
# Phase 3の立ち上がり部分を指す
ax.annotate('1. Residual Stress\n(Activation at Rest)', 
            xy=(1.0, 30), xytext=(1.02, 30),
            arrowprops=dict(facecolor='black', shrink=0.05),
            fontsize=10, color='red', fontweight='bold')

# 2. 剛性の増加 (Stiffening)
# Active Loadの傾きが急であることを示す
# 比較点としてLambda=1.1付近を使用
stress_zhao = data.loc[(data["Lambda"] > 1.09) & (data["Lambda"] < 1.11), "Stress_Zhao(kPa)"].mean()
stress_integ = data.loc[(data["Phase_Label"]=="Active_Load") & (data["Lambda"] > 1.09) & (data["Lambda"] < 1.11), "Stress_Integrated(kPa)"].mean()

ax.annotate('2. Increased Stiffness\n(Steeper Slope)', 
            xy=(1.1, stress_integ), xytext=(1.05, stress_integ + 20),
            arrowprops=dict(facecolor='black', shrink=0.05),
            fontsize=10, color='darkred', fontweight='bold')

# --- グラフの装飾 ---
ax.set_title("Comparison of Constitutive Models: Passive vs Active Response", fontsize=14)
ax.set_xlabel("Stretch Ratio $\lambda$", fontsize=12)
ax.set_ylabel("Cauchy Stress $\sigma$ (kPa)", fontsize=12)
ax.set_xlim(0.98, 1.17) # ここで範囲指定 (エラーの原因だった部分をPythonで制御)
ax.set_ylim(0, None)
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc="upper left")

# 保存と表示
plt.tight_layout()
plt.savefig("comparison_final_python.png", dpi=300)
print("Graph saved as 'comparison_final_python.png'")