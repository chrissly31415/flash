#!/usr/bin/gnuplot 
#eps versions
set terminal postscript eps color enhanced font "helvetica, 16" 
set output "ads1.eps"
#set rmargin 0.5
#set lmargin 0.5
#unset label
set nokey
set hidden3d
#set parametric
#set isosamples 2,33
#set surface
#set key center top title " "
set zlabel "{/Symbol j}\nvolume fraction" offset -5.0 
#set zrange [0:0.1]
set xlabel "z, distance from surface" 
set ytics scale 0.2 offset 0 
set ylabel "X_{n}, degree of polymerization"
set xtics 
set ytics
set ztics
set ticslevel 0.5
#set autoscale y 
#set autoscale x
#view angle and distance and height 
#set view 45,60,0.8,1.5
# x,y and weight factor
#map
#set dgrid3d 23,20,100 
#set dgrid3d 23,19,90 
#set dgrid3d 20,320,80 
#set contour base
#show contour
#introduces colout
#set palette rgbformulae 7,5,0 
#set palette defined (0 "dark-blue", 1 "blue", 3 "green", 4 "yellow", 5 "red", 6 "dark-red")
set palette defined (1 "blue", 3 "green", 4 "yellow", 8 "red")
#set xrange [0:20]
#set noxtics
#set noytics
#set autoscale x 
set cblabel "{/Symbol j}, volume fraction" offset -1.0 font "helvetica, 12" 
set cbtics font "helvetica, 10"
#set colorbox vertical user origin .1,.02 size .8,.04 
#set format cb '%.3f'
set sample 6
set isosamples 6 
set xrange [0:9]
set yrange [0:]
set cbrange [:0.5]
set pm3d map  
set pm3d flush begin
set pm3d flush begin noftriangles scansforward
#set pm3d at s 
splot 'data_fi3D' 
#set output "ads1.png"
#set terminal png  
#pause 10000 
#pause -1 "Hit RETURN to continue"
