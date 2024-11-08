set terminal gif animate delay 10 optimize size 800,600
set output 'ACE-equation.gif'

# Nastavení titulku a popisků os
set title "ACE equation"
set xlabel "X"
set ylabel "Y"
set zlabel "f(x,y)"
set size ratio -1

# Povolení mřížky
set pm3d map
set palette rgbformulae 33,13,10
set cbrange [0:1]

set grid

# Vykreslení povrchu
splot 'Results/ACE-equation-00005.txt' using 1:2:3 with pm3d title 'f(x,y)'
do for [i=0:1000:10]{
    filename = sprintf('Results/ACE-equation-%05d.txt', i)
    set title sprintf("ACE-equation - Frame %d", i)
    splot filename using 1:2:3 with pm3d notitle,\
          filename using 1:2:(($3 > 0.49 && $3 < 0.51) ? 0.5 : 1/0) with points pointtype 7 pointsize 0.5 linecolor "black" notitle
}