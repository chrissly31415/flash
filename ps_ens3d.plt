#!/usr/bin/gnuplot 
#eps versions
#set terminal postscript eps color enhanced font "helvetica, 16" 
#set output "ads1.eps"
#set rmargin 0.5
#set lmargin 0.5
#unset label
set style data lines 
set style data dots 
set nokey
set hidden3d
#set parametric
#set isosamples 2,33
set surface
#set key center top title " "
set zlabel "volume\nfraction" offset -8.0 
#set zrange [0:0.5]
set xlabel "z" offset 2.0 
set ytics scale 0.8 offset 0 
set ylabel "degree of polymerization" offset 3.0
set xtics 
set ytics
set ztics
set ticslevel 0.5
#set autoscale y 
#set autoscale x
#view angle and distance and height 
set view 50,80,0.8,1.5
# x,y and weight factor
#map
#set dgrid3d 23,20,100 
set dgrid3d 50,50,20 
#set dgrid3d 20,320,80 
#set contour base
#show contour
#introduces colout
#set palette rgbformulae 7,5,0 
#set palette defined (0 "dark-blue", 1 "blue", 3 "green", 4 "yellow", 5 "red", 6 "dark-red")
#set palette defined (1 "blue", 3 "green", 4 "yellow", 8 "red")
#set xrange [0:20]
#set noxtics
#set noytics
#set autoscale x 
#set cblabel "{/Symbol j}, volume fraction" -1.0 font "helvetica, 12" 
#set cbtics font "helvetica, 10"
#unset colorbox
#set format cb '%.3f'
#set sample 20 
#set isosamples 20 
set xrange [0:9]
#set yrange [0:]
#set cbrange [:0.5]
#set pm3d  
#set pm3d flush begin
#set pm3d flush begin noftriangles scansforward
set pm3d at s 
splot 'data_fi3D' 
#set output "ads1.png"
#set terminal png  
#pause 10000 
pause -1 "Hit RETURN to continue"
