# ============================================================
# Replicate Paetsch et al. (2012) Figure 5
# ============================================================
set terminal pngcairo size 1000,500 enhanced font 'Arial,12'
set output 'Figure5_Replication.png'
set datafile separator ","

set multiplot layout 1,2 title "Replication of Figure 5: Passive-Active Transition"

# --- 1. Stress vs Stretch ---
set title "Cauchy Stress vs Stretch"
set xlabel "Stretch {/Symbol l}"
set ylabel "Stress (MPa)"  # 論文に合わせてMPa表示にする
set grid
set key top left

# 列4(Stress Pa)を 1e6 で割って MPa に変換
plot 'passive_active_transition.csv' using 2:($4/1e6) every ::1 w l lw 3 lc rgb "blue" title "Model Response"

# --- 2. Energy vs Stretch ---
set title "Pseudo-Energy vs Stretch"
set xlabel "Stretch {/Symbol l}"
set ylabel "Energy (kPa)"
set grid

# 列5(Energy Pa)を 1000 で割って kPa に変換
plot 'passive_active_transition.csv' using 2:($5/1000) every ::1 w l lw 3 lc rgb "red" title "Strain Energy"

unset multiplot