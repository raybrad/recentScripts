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
        echo "sed"
	sed -i "s/ext_sta=250/ext_sta=1/g" $tmp/in.td
	sed -i "s/ref_sta=250/ref_sta=1/g" $tmp/in.td
    done
    done
   done
  done
