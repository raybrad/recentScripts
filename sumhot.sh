#!/bin/sh
for k in $(seq 0 1000 99000 )
do
#    awk '{ sum += $2; n++ } END { if (n > 0) print sum; }' diagMO.${k} >> diag.sum
     awk 'BEGIN {sum = 0} {if ($1>-4.418) {sum +=$2}} END {print  sum}' diagMO.${k} >>hot.sum
done
