#!/bin/sh
for i in  $(seq 1 5)  
do
 for j in 1 2 3 
 do
   for k in 0 2 
   do  
     for m in 1 2 3 
     do  
      if [ $i -le 9 ]
      then
           tmp=0${i}${j}${k}${m}
      else
           tmp=${i}${j}${k}${m}
      fi  
        echo "remove"
        rm $tmp/voltage.data $tmp/dip.data  $tmp/deltaq_atom.dat  $tmp/qtd$tmp.hy.o* $tmp/qtd$tmp.hy.po*
        rm $tmp/hamiloverl.data $tmp/dipole.data  $tmp/system.data $tmp/denmatrix.data $tmp/overl.dat
    done
    done
   done
  done
