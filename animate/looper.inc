j = i + 0
t = 0.01 * j
sp \
   'bond.dat' u 2:3:(stringcolumn(1) eq 'H' ? $4 : 1/0) w p pt 7 ps 0.25 lc rgb 'gray' notit, \
   'bond.dat' u 2:3:(stringcolumn(1) eq 'C' ? $4 : 1/0) w p pt 7 ps 0.25 lc rgb 'gray' notit, \
   'bond.dat' u 2:3:(stringcolumn(1) eq 'S' ? $4 : 1/0) w p pt 7 ps 0.25 lc rgb 'gray' notit, \
    "dat/c".j.".dat" u 1:2:3:4 w p palette ps 1 pt 7 t sprintf("Time: %3.2f fs", t)
# step
i = i + 10
# end point
if (i < 10000) reread

#   'atom.dat' u 2:3:(stringcolumn(1) eq 'H' ? $4 : 1/0) w p pt 7 ps 0.7 lc rgb 'black'  notit, \
#   'atom.dat' u 2:3:(stringcolumn(1) eq 'C' ? $4 : 1/0) w p pt 7 ps 1.0 lc rgb 'orange' notit, \
#   'atom.dat' u 2:3:(stringcolumn(1) eq 'S' ? $4 : 1/0) w p pt 7 ps 1.5 lc rgb 'yellow' notit, \

