#!/usr/bin/gnuplot

set terminal pdfcairo size 6cm, 6cm \
    transparent enhanced color \
    font "Latin Modern Math,10"

load "third_party/gnuplot-palettes/moreland.pal"

# set xrange [0.8:1]
# set xtics 0.8,0.05,1.0
# set yrange [100:5500]
# set logscale y
# set border lt 1 lw 1.5 lc "black" back
# set key bottom right
# set grid
# set lmargin at screen 0.15
# set rmargin at screen 0.95
# set tmargin at screen 1
# set bmargin at screen 0.15
# set xrange [:10]

# set xlabel "τ, × system size"
# set ylabel "C(τ)"
# set datafile separator ','
# unset border
# unset xtics
# unset ytics
set size ratio 1.0
set xtics ("-4π" -4 * pi, "-2π" -2 * pi, "0" 0, "2π" 2 * pi, "4π" 4 * pi)
set ytics ("-4π" -4 * pi, "-2π" -2 * pi, "0" 0, "2π" 2 * pi, "4π" 4 * pi)
set cbtics scale 0.5

# set cbrange [0:1]
# set title "λ=2.5"
n = 25
λ = ARG1 # "4.0" # "1.7320508" # 1.5
β = ARG2 # 0.2
k = 201 # 51
period = 8 * pi
basename = sprintf("data/mean/new/structure_fourier_n=%d_λ=%s_β=%s", n, λ, β, n, λ, β)
input_filename = basename . ".dat"
output_filename = basename . ".pdf"

set output output_filename

f(q) = -0.5 * period + (q / (k - 1.0)) * period
plot input_filename matrix u (f($1)):(f($2)):3 with image notitle

set output
command = sprintf("convert -density 600 %s -quality 00 %s.png", output_filename, basename)
system command
