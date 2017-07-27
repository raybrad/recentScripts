#!/bin/sh
myarr=($(awk  '{if($1>0.0) print NR}' test.dat))
for i in "${myarr[@]}"
do
         echo $i
done
