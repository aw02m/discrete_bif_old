set encoding utf8

set terminal pdf font "Nimbus Roman, 16"
set output "output.pdf"
set minussign

set size square
set colorsequence classic

set xlabel "a"
set ylabel "d"

set xrange[-5:-1]
set yrange[0.4:1]

set xtics 1.0
set ytics 0.1

set grid
unset key

plot "./build/out" using 1:2 w l lc rgb "black"

set terminal x11
replot
pause -1
