if (!exists("frame")) frame = 0
if (!exists("fixed_color_range")) fixed_color_range = 1
if (!exists("interface")) interface = 1
if (!exists("folder")) folder = "Results"
output_file = sprintf('ACE-equation-%05d.png', frame)

set terminal png size 800,600
set output output_file

# Nastavení titulku a popisků os
set title "ACE equation"
set xlabel "X"
set ylabel "Y"
set zlabel "f(x,y)"
set size ratio -1

# Povolení mřížky
set pm3d map
set palette rgbformulae 33,13,10
if(fixed_color_range) set cbrange[0:1] 

# Vykreslení povrchu
filename = sprintf('%s/ACE-equation-%05d.txt', folder, frame)
set title sprintf("ACE-equation - Frame %d", frame)
if(interface) {
    splot filename using 1:2:3 with pm3d notitle,\
          filename using 1:2:(($3 > 0.499 && $3 < 0.501) ? 0.5 : 1/0) with points pointtype 7 pointsize 0.5 linecolor "black" notitle

} else {
    splot filename using 1:2:3 with pm3d notitle
}