set title "Comparison of Zhao Model vs. Integrated Model\n(Demonstrating Active Stiffening & Residual Stress)" font ",14"
set xlabel "Stretch Ratio {/Symbol l}" font ",14"
set ylabel "Cauchy Stress {/Symbol s} (kPa)" font ",14"

set grid lc rgb "#dddddd" lt 1
set key top left box opaque font ",11" spacing 1.5

# 軸の範囲
set xrange [0.98:1.16]
set yrange [0:*]

# --- 重要な特徴を強調する矢印とラベル ---

# 1. 残留応力 (Residual Stress) の発生
set arrow from 1.0, 5 to 1.0, 45 head filled lc rgb "red" lw 2
set label "1. Activation at rest\n(Residual Stress)" at 1.01, 25 textcolor rgb "red" font ",10"

# 2. 剛性 (Stiffness) の増加
# Active Loadの傾き付近
set arrow from 1.08, 60 to 1.05, 90 heads size 0.5,90 front lc rgb "black" lw 1 dt 2
set label "2. Increased Stiffness\n(Steeper Slope)" at 1.02, 85 textcolor rgb "dark-red" font ",10"

set datafile separator ","

# プロット内容:
# 1. Zhao Model (Passive): 青い破線。刺激が入っても変わらない基準線。
# 2. Integrated (Passive Phase): 灰色。初期のPassive挙動。
# 3. Integrated (Active Phase): 赤い実線。Activation後の挙動。

plot \
    "comparison_experiment.csv" using 1:4 every 5 with linespoints pt 7 ps 0.5 lc rgb "blue" dt 2 title "Zhao Model (Passive Only)", \
    "comparison_experiment.csv" using (stringcolumn(6) eq "Passive_Load" ? $2 : 1/0):5 with lines lw 4 lc rgb "gray" title "Integrated (Passive State)", \
    "comparison_experiment.csv" using (stringcolumn(6) eq "Activation" ? $2 : 1/0):5 with lines lw 4 lc rgb "red" title "Integrated (Activation Phase)", \
    "comparison_experiment.csv" using (stringcolumn(6) eq "Active_Load" ? $2 : 1/0):5 with lines lw 4 lc rgb "dark-red" title "Integrated (Active Loading)", \
    "comparison_experiment.csv" using (stringcolumn(6) eq "Active_Unload" ? $2 : 1/0):5 with lines lw 3 dt 1 lc rgb "pink" title "Integrated (Active Unloading)"