# 出力設定 (プレゼン資料用に解像度を高めに設定)
set terminal pngcairo size 800,600 enhanced font "Arial,14" transparent
set output "stress_comparison.png"

# CSVファイルの読み込み設定 (1行目はヘッダーなのでスキップ)
set datafile separator ","
set key left top box
set grid

# 軸のラベル設定
set xlabel "Stretch {/Symbol l} [-]"
set ylabel "Cauchy Stress [kPa]"

# 描画範囲の調整 (必要に応じて変更してください)
set xrange [1.0:1.10]
set yrange [0:350]

# 線のスタイル定義
# Neo-Hookean (黒の破線)
set style line 1 lc rgb "black" dt 2 lw 2
# 分析解_円周方向 (青の実線)
set style line 2 lc rgb "blue" lw 2
# 分析解_軸方向 (赤の実線)
set style line 3 lc rgb "red" lw 2
# シミュレーション_円周方向 (青の点)
set style line 4 lc rgb "blue" pt 7 ps 1.2
# シミュレーション_軸方向 (赤の点)
set style line 5 lc rgb "red" pt 7 ps 1.2

# プロット実行
plot "stress_results.csv" every ::1 using 1:6 with lines ls 1 title "Neo-Hookean (Isotropic)", \
     "stress_results.csv" every ::1 using 1:3 with lines ls 2 title "Zhao theoretical ({/Symbol q})", \
     "stress_results.csv" every ::1 using 1:5 with lines ls 3 title "Zhao theoretical (z)", \
     "stress_results.csv" every ::1 using 1:2 with points ls 4 title "IGA simulation ({/Symbol q})", \
     "stress_results.csv" every ::1 using 1:4 with points ls 5 title "IGA simulation (z)"