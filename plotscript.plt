#!/usr/bin/gnuplot -persist
#
set key center top title " "
set ylabel "concentration"
set y2label "energy"
set xlabel "z"
set grid x y
set ytics nomirror
set tics out
set y2tics
set autoscale  y
set yrange [0.0:1.0]
#set autoscale y2
plot 'data' using 1:2 with linespoints t 'fi(A)' axis x1y1, 'data' using 1:3 with linespoints t 'fi(B)' lt -1 axis x1y1, 'data' using 1:4 with linespoints t 'u(A)' axis x1y2, 'data' using 1:5 with linespoints t 'u(B)' axis x1y2  
#pause 5
#reread
pause 100 
#pause -1 "Hit RETURN to continue"
