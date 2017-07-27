set xlabel 'time(fs)
#set ylabel 'Voltage (arb)'
#set ylabel 'dAt (arb)'
set ylabel 'div dAt (arb)'
unset ytics
#p '../../../../nometal/nometal_epdf1e-15/dump/Vzline.dat' u ($1*1e15):($4*1e20) w l  t 'no metal(above)', 'Vzline.dat' u ($1*1e15):($4*1e20) w l  t 'metal(above)
#rep '../../../../nometal/nometal_epdf1e-15/dump/Vzline.dat' u ($1*1e15):($6*1e20) w l t 'no metal(in)', 'Vzline.dat' u ($1*1e15):($6*1e20) w l t 'metal(in)'
#rep '../../../../nometal/nometal_epdf1e-15/dump/Vzline.dat' u ($1*1e15):($8*1e20) w l t 'no metal(below)', 'Vzline.dat' u ($1*1e15):($8*1e20) w l t 'metal(below)'


#p '../../../../nometal/nometal_epdf1e-15/dump/EzdAtline.dat' u ($1*1e15):($4*1e20) w l  t 'no metal(above)', 'EzdAtline.dat' u ($1*1e15):($4*1e20) w l  t 'metal(above)
#rep '../../../../nometal/nometal_epdf1e-15/dump/EzdAtline.dat' u ($1*1e15):($6*1e20) w l t 'no metal(in)', 'EzdAtline.dat' u ($1*1e15):($6*1e20) w l t 'metal(in)'
#rep '../../../../nometal/nometal_epdf1e-15/dump/EzdAtline.dat' u ($1*1e15):($8*1e20) w l t 'no metal(below)', 'EzdAtline.dat' u ($1*1e15):($8*1e20) w l t 'metal(below)

p '../../../../nometal/nometal_epdf1e-15/dump/EzdivdAtline.dat' u ($1*1e15):($4*1e20) w l  t 'no metal(above)', 'EzdivdAtline.dat' u ($1*1e15):($4*1e20) w p pt 5 t 'metal(above)
