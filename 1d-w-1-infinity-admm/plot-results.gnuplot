#!/bin/bash gnuplot 

set terminal png size 1024, 580

set grid
set xlabel 'x'
set ylabel "u, u', f, g, q"
set output "out/results.png"
plot 'out/point_data.dat' using 1:2 with lines linewidth 3 title 'Solution u',\
     'out/point_data.dat' using 1:3 with lines linewidth 2 title 'Right hand side f',\
     'out/point_data.dat' using 1:4 with lines linewidth 2 title 'Right hand side g',\
     'out/element_data.dat' using 1:2 with lines linewidth 2 title "Gradient u'",\
     'out/element_data.dat' using 1:3 with points pointsize 2 pointtype 2 title 'Slack variable q'

set output 'out/objective.png'
set title 'Normalized objective'
unset key
set xlabel '# Steps'
set ylabel "J_k / |J_0|"
plot 'out/obj.dat' using 1:2 with lines linewidth 3


set output "out/residual.png"
set log y
set title 'Normalized residual'
set ylabel 'R_k / |R_0|'
plot 'out/res.dat' using 1:2 with lines linewidth 3
