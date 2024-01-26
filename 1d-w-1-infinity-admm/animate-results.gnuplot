#!/bin/bash gnuplot 

set terminal png size 1024, 580 # HD

system("rm -rf out/animation/*")
system("mkdir -p out/animation")

set yrange [-4:6]
set xrange [0:2*pi]
set ylabel "Solutions"
set xlabel "x"
set grid 

n=system("ls -1 out/steps/element_data* | wc -l")

do for [j=1:n-1] {
     set output sprintf("out/animation/%04d.png", j)
     set title sprintf("Step: %d", j)
     print(sprintf("Creating frame for step %d", j))
     plot sprintf("out/steps/point_data_%04d.dat", j) using 1:2 with lines linewidth 3 title 'Solution u',\
          sprintf("out/steps/element_data_%04d.dat", j) using 1:2 with lines linewidth 2 title 'Gradient du', \
          sprintf("out/steps/element_data_%04d.dat", j) using 1:3 with points pointsize 1 title 'Slack q', \
          sprintf("out/steps/element_data_%04d.dat", j) using 1:4 with points pointsize 1 title 'Multiplier lambda'
}