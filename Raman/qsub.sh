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
        echo "qsub"	
	cd $tmp
 	qsub qtd$tmp.hy
	cd ..
    done
    done
   done
  done
