#!/usr/bin/gnuplot

set terminal pdfcairo size 8cm, 6cm \
    transparent enhanced color \
    font "Latin Modern Math,10"

load "third_party/gnuplot-palettes/moreland.pal"

# set xrange [0.8:1]
# set xtics 0.8,0.05,1.0
# set yrange [100:5500]
# set logscale y
# set border lt 1 lw 1.5 lc "black" back
set key bottom right
# set grid
# set lmargin at screen 0.15
# set rmargin at screen 0.95
# set tmargin at screen 1
# set bmargin at screen 0.15
# set xrange [:10]

set xlabel "τ, × system size"
set ylabel "C(τ)"
set datafile separator ','

# set cbrange [0:1]
# set title "λ=2.5"
n = 25
λ = ARG1 # "2.5" # "1.7320508" # 1.5
input_filename(β) = sprintf("data/mean/autocorr_n=%d_λ=%s_β=%.1f.csv", n, λ, β)
basename = sprintf("data/mean/autocorr_n=%d_λ=%s", n, λ)
output_filename = basename . ".pdf"

set output output_filename

# set yrange [1e-3:1]
# set logscale y
# plot [:30000] \
#     'data/autocorr_n=25_λ=7.5_β=0.29_seed=47588.csv' u 0:1 w l title "λ=7.5, β=0.29", \
#     'data/autocorr_n=25_λ=7.5_β=0.3_seed=47588.csv' u 0:1 w l title "λ=7.5, β=0.3", \
#     'data/autocorr_n=25_λ=7.5_β=0.35_seed=47588.csv' u 0:1 w l title "λ=7.5, β=0.35", \
#     'data/autocorr_n=25_λ=7.5_β=0.7_seed=47588.csv' u 0:1 w l title "λ=7.5, β=0.7"

# plot \
#     'data/autocorr_n=25_λ=2.5_β=0.2_seed=47589.csv' u 0:1 w l title "λ=2.5, β=0.2", \
#     'data/autocorr_n=25_λ=2.5_β=0.25_seed=47589.csv' u 0:1 w l title "λ=2.5, β=0.25", \
#     'data/autocorr_n=25_λ=2.5_β=0.35_seed=47589.csv' u 0:1 w l title "λ=2.5, β=0.35", \
#     'data/autocorr_n=25_λ=2.5_β=0.45_seed=47589.csv' u 0:1 w l title "λ=2.5, β=0.45", \
#     'data/autocorr_n=25_λ=2.5_β=0.55_seed=47589.csv' u 0:1 w l title "λ=2.5, β=0.55"

# plot \
#     'data/autocorr_n=25_λ=2.5_β=0.6_seed=47590.csv' u 0:1 w l title "λ=2.5, β=0.6", \
#     'data/autocorr_n=25_λ=2.5_β=0.8_seed=47590.csv' u 0:1 w l title "λ=2.5, β=0.8", \
#     'data/autocorr_n=25_λ=2.5_β=1.0_seed=47590.csv' u 0:1 w l title "λ=2.5, β=1.0", \
#     'data/autocorr_n=25_λ=2.5_β=1.2_seed=47590.csv' u 0:1 w l title "λ=2.5, β=1.2"

# plot [:200000] \
#     'remote/autocorr_n=25_λ=7.5_β=0.27_t=32_seed=47607.csv' u 0:1 w l title "t_w=32", \
#     'remote/autocorr_n=25_λ=7.5_β=0.27_t=128_seed=47607.csv' u 0:1 w l title "t_w=128", \
#     'remote/autocorr_n=25_λ=7.5_β=0.27_t=512_seed=47607.csv' u 0:1 w l title "t_w=512", \
#     'remote/autocorr_n=25_λ=7.5_β=0.27_t=2048_seed=47607.csv' u 0:1 w l title "t_w=2048", \
#     'remote/autocorr_n=25_λ=7.5_β=0.27_t=8192_seed=47607.csv' u 0:1 w l title "t_w=8192", \
#     'remote/autocorr_n=25_λ=7.5_β=0.27_t=32768_seed=47607.csv' u 0:1 w l title "t_w=32768"

# plot \
#     'data/autocorr_n=400_λ=0.0_β=0.15_t=32_seed=1234.csv' u 0:1 w l title "t_w=32", \
#     'data/autocorr_n=400_λ=0.0_β=0.15_t=128_seed=1234.csv' u 0:1 w l title "t_w=128", \
#     'data/autocorr_n=400_λ=0.0_β=0.15_t=512_seed=1234.csv' u 0:1 w l title "t_w=512", \
#     'data/autocorr_n=400_λ=0.0_β=0.15_t=2048_seed=1234.csv' u 0:1 w l title "t_w=2048", \
#     'data/autocorr_n=400_λ=0.0_β=0.15_t=8192_seed=1234.csv' u 0:1 w l title "t_w=8192"

# δt = 100
# 
# plot \
#     'data/mean_autocorr_n=400_β=1.0_t=32.csv' u (δt * $0):1 every δt w l title "t_w=32", \
#     'data/mean_autocorr_n=400_β=1.0_t=128.csv' u (δt * $0):1 every δt w l title "t_w=128", \
#     'data/mean_autocorr_n=400_β=1.0_t=512.csv' u (δt * $0):1 every δt w l title "t_w=512", \
#     'data/mean_autocorr_n=400_β=1.0_t=2048.csv' u (δt * $0):1 every δt w l title "t_w=2048", \
#     'data/mean_autocorr_n=400_β=1.0_t=8192.csv' u (δt * $0):1 every δt w l title "t_w=8192", \
#     'data/mean_autocorr_n=400_β=1.0_t=32768.csv' u (δt * $0):1 every δt w l title "t_w=32768"

# set logscale y
# plot [0:1000][5e-5:1e0] \
#     'data/mean_autocorr_n=400_β=0.2.csv' u 0:1 w l ls 1 title "β=0.2", \
#     'data/mean_autocorr_n=400_β=0.6.csv' u 0:1 w l ls 2 title "β=0.6", \
#     'data/mean_autocorr_n=400_β=1.0.csv' u 0:1 w l ls 3 title "β=1.0", \
#     'data/mean_autocorr_n=400_β=1.4.csv' u 0:1 w l ls 4 title "β=1.4", \
#     'data/mean_autocorr_n=400_β=1.8.csv' u 0:1 w l ls 5 title "β=1.8", \
#     'data/mean_autocorr_n=400_β=2.2.csv' u 0:1 w l ls 6 title "β=2.2", \
#     'data/mean_autocorr_n=400_β=2.6.csv' u 0:1 w l ls 7 title "β=2.6", \
#     'data/mean_autocorr_n=400_β=3.0.csv' u 0:1 w l ls 8 title "β=3.0"

# set title "λ=2.5"
set logscale y
# [0:10][5e-4:1e0]
plot [0:1000][1e-4:1e0] \
    input_filename(0.2) u 0:1 w l ls 1 lw 3 title "β=0.2", \
    input_filename(0.6) u 0:1 w l ls 2 lw 3 title "β=0.6", \
    input_filename(1.0) u 0:1 w l ls 3 lw 3 title "β=1.0", \
    input_filename(1.4) u 0:1 w l ls 4 lw 3 title "β=1.4", \
    input_filename(2.2) u 0:1 w l ls 6 lw 3 title "β=2.2", \
    input_filename(3.0) u 0:1 w l ls 8 lw 3 title "β=3.0"

set output
command = sprintf("convert -density 600 %s -quality 00 %s.png", output_filename, basename)
system command

# 'data/autocorr_n=400_β=0.15_t=128_seed=1234.csv' u 0:1 w l title "t_w=128",
# 'data/autocorr_n=400_β=0.15_t=512_seed=1234.csv' u 0:1 w l title "t_w=512",
# 'data/autocorr_n=400_β=0.15_t=2048_seed=1234.csv' u 0:1 w l title "t_w=2048",
# 'data/autocorr_n=400_β=0.15_t=8192_seed=1234.csv' u 0:1 w l title "t_w=8192"


# plot [:30000] \
#     'data/autocorr_n=10_λ=20.0_β=0.1.csv' u 0:1 w l, \
#     'data/autocorr_n=10_λ=20.0_β=0.11.csv' u 0:1 w l, \
#     'data/autocorr_n=10_λ=20.0_β=0.12.csv' u 0:1 w l, \
#     'data/autocorr_n=10_λ=20.0_β=0.13.csv' u 0:1 w l, \
#     'data/autocorr_n=10_λ=20.0_β=0.14.csv' u 0:1 w l, \
#     'data/autocorr_n=10_λ=20.0_β=0.15.csv' u 0:1 w l, \
#     'data/autocorr_n=10_λ=20.0_β=0.16.csv' u 0:1 w l


# plot "data/observables_n=10_λ=20.0_β=0.1" u 0:2

# set output
# system("for m in pyrochlore kagome sk; do \
#           convert -density 600 scatter_${m}.pdf -quality 00 scatter_${m}.png; \
#         done")


