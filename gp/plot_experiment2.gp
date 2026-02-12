set terminal pngcairo size 1000,600 enhanced font 'Arial,12'
set output 'experiment2_revised.png'

set title "Active Force-Length Relationship (Improved Visualization)" font ",14"
set xlabel "Stretch Ratio {/Symbol l}" font ",14"
set ylabel "Cauchy Stress {/Symbol s} (kPa)" font ",14"

set grid
set key top left box opaque font ",10"

# 重要: Y軸の範囲を固定して、Active成分が見えるようにする
# Passiveが急上昇しても、関心領域（0-200kPa）にフォーカスする
set yrange [-50:250] 
set xrange [0.7:1.4]

# ゼロラインを強調
set xzeroaxis lt 1 lc rgb "black" lw 1

set datafile separator ","

plot "experiment2_revised.csv" every ::1 using 1:2 w l lw 2 dt 2 lc rgb "blue" ti "Passive Tension (Matrix + Fiber)", \
     "experiment2_revised.csv" every ::1 using 1:3 w l lw 2 lc rgb "black" ti "Total Tension (Active + Passive)", \
     "experiment2_revised.csv" every ::1 using 1:4 w l lw 4 lc rgb "red" ti "Active Tension (Contractile)"
