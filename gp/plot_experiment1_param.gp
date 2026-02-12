# --- 設定開始 ---
set datafile separator ","

# 出力形式の設定 (PNG画像として保存)
set terminal pngcairo size 900,600 enhanced font "Arial,12"
set output "passive_biaxial_comparison.png"

# グラフのタイトルと軸ラベル
set title "Passive Biaxial Extension Test Comparison" font ",14"
set xlabel "Stretch Ratio ({/Symbol l})" font ",12"
set ylabel "Cauchy Stress (kPa)" font ",12"

# グリッドと凡例の設定
set grid
set key top left box opaque width 2

# --- 線種スタイルの定義 ---

# 1. Neo-Hookean (灰色) - 比較対象
set style line 11 lc rgb "gray40" lw 2 lt 1  # 実線 (Circ)
set style line 12 lc rgb "gray40" lw 2 dt 2  # 点線 (Axial)

# 2. Paetsch et al. (青色) - 多項式型
set style line 21 lc rgb "blue"   lw 2 lt 1  # 実線 (Circ)
set style line 22 lc rgb "blue"   lw 2 dt 2  # 点線 (Axial)

# 3. Zhao (赤色) - 提案モデル (指数型) ※線を少し太く強調
set style line 31 lc rgb "red"    lw 3 lt 1  # 実線 (Circ)
set style line 32 lc rgb "red"    lw 3 dt 2  # 点線 (Axial)

# --- プロット実行 ---
# CSV列: 1:Lam, 2:Neo_C, 3:Neo_A, 4:Pae_C, 5:Pae_A, 6:Zha_C, 7:Zha_A

plot "experiment1_passive_adjusted.csv" using 1:2 with lines ls 11 title "Neo-Hookean (Circ)", \
     "" using 1:3 with lines ls 12 title "Neo-Hookean (Axial)", \
     "" using 1:4 with lines ls 21 title "Paetsch (Circ)", \
     "" using 1:5 with lines ls 22 title "Paetsch (Axial)", \
     "" using 1:6 with lines ls 31 title "Zhao (Circ)", \
     "" using 1:7 with lines ls 32 title "Zhao (Axial)"

# --- 設定終了 ---