#!/bin/sh
for i in  $(seq 1 11)  
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
        rm $tmp/eField.dat $tmp/dip_nLoc.1.data  $tmp/timer.data  $tmp/qtd$tmp.hy.o* $tmp/qtd$tmp.hy.po*
    done
    done
   done
  done
