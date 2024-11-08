set terminal png size 800,600
set output 'ACE-equation.png'

# Nastavení titulku a popisků os
set title "ACE equation"
set xlabel "X"
set ylabel "Y"
set zlabel "f(x,y)"
set size ratio -1

# Povolení mřížky
set pm3d map
set palette rgbformulae 33,13,10

# Vykreslení povrchu
splot 'Results/ACE-equation-01000.txt' using 1:2:3 with pm3d title 'f(x,y)'