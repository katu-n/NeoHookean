# ==========================================
# Experiment 1: Passive Biaxial Extension Comparison
# ==========================================

set terminal pngcairo size 1000,600 enhanced font 'Arial,12'
set output 'experiment1_passive.png'

set title "Experiment 1: Passive Biaxial Extension (Equibiaxial)\nComparison of Constitutive Models" font ",14"
set xlabel "Stretch Ratio {/Symbol l} (= {/Symbol l}_{circ} = {/Symbol l}_{axial})" font ",14"
set ylabel "Cauchy Stress {/Symbol s} (kPa)" font ",14"

set grid
set key top left box opaque font ",10"

set datafile separator ","

# 線種の設定
# Neo-Hookean: Black
# Paetsch: Blue (Circ), Cyan (Axial)
# Zhao: Red (Circ), Orange (Axial)

plot "experiment1_passive.csv" every ::1 using 1:2 w l lw 2 dt 2 lc rgb "gray" ti "Neo-Hookean (Isotropic)", \
     "experiment1_passive.csv" every ::1 using 1:4 w l lw 2 lc rgb "blue" ti "Paetsch (Poly) - Circ", \
     "experiment1_passive.csv" every ::1 using 1:5 w l lw 2 dt 3 lc rgb "blue" ti "Paetsch (Poly) - Axial", \
     "experiment1_passive.csv" every ::1 using 1:6 w l lw 3 lc rgb "red" ti "Zhao (Exp) - Circ", \
     "experiment1_passive.csv" every ::1 using 1:7 w l lw 3 dt 3 lc rgb "red" ti "Zhao (Exp) - Axial"