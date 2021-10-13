set terminal pngcairo  transparent enhanced font "arial,14" fontscale 1.7 size 1200, 600 
set output 'graphics/radar.png'
set clip points
unset border
set raxis
set angles degrees
set polar
set style data lines
set multiplot layout 1, 2 title "Positionnement des scenarii possibles" font ",18"
#set xzeroaxis
#set yzeroaxis
#set zzeroaxis
unset xtics
unset ytics 
unset rtics
set size ratio -1
set style line 1 linetype 1 linecolor rgb "#265c75"
set style line 2 linetype 1 linecolor rgb "#4bbdf2"
set style line 3 linetype 3 linecolor rgb "grey" linewidth 0.5
set rrange [0:6]
set label 1 center rotate by -30 at 1.8,3.12 "Performance"
set label 2 center at -3.6,0 rotate by 90  "Peu coûteux"
set label 3 center at 1.8,-3.12 rotate by 30  "Maitrise des couches basses"

set arrow from 0,0 to 6,0 nohead  filled front lw 3 lc rgb "black"
set arrow from 0,0 to -3,5.196152422706632  nohead filled front lw 3 lc rgb "black"
set arrow from 0,0 to -3,-5.196152422706632  nohead filled front lw 3 lc rgb "black"
set title "Scenario 1 :\n base open source"
plot 'radar.dat' using 1:2 with filledcurves below ls 2 notitle , 'radar.dat' using 1:4 ls 3 notitle
set title "Scenario 2 :\n base CEA"
plot 'radar.dat' using 1:3 with filledcurves below ls 1 notitle , 'radar.dat' using 1:4 ls 3 notitle
#cos 0 : 1
# sin 0 : 0
#cos 120 : -0.5
# sin 120 0.8660254037844387   4.330127018922194
#cos 240 -0.5
# sin 240

