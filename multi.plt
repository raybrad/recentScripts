set terminal postscript eps size 8.0,6.0 enhanced color font 'Helvetica,20' lw 2
set output 'temp.eps'
load 'matlab.pal'

set multiplot layout 1,2

set lmargin 5
set rmargin 5
set tmargin 5
set bmargin 5
set size square
set border
set xtics format "%2.1f"  
set ytics format "%2.1f"  
set xlabel 'eV'
set ylabel '|E|/|E_{inc}|'
set xlabel 'nm'
set xr [300:800]
set yr [0.8:1.6]
set key off
p '../monitorScan.dat' u (1239.84187/$1):($4)  every 2 w lp pt 6 ps 2 t 'QM'

set lmargin 5
set rmargin 5
set tmargin 5
set bmargin 5
set size square
set pm3d map
#set contour
#set cntrparam levels 20
#set cntrparam bspline
#set cbr[0:4]
unset border
unset xtics
unset ytics
unset xlabel
unset ylabel
set xr[-20:20]
set yr[-20:20]
set format x ''
set format y ''
sp 'xyField.dat10000' u 1:2:3  every 2 w  image



unset multiplot
