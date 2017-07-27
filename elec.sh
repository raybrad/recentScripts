#!/bin/sh
for k in $(seq 0 1000 99000 )
do
    awk '{ sum += $2; n++ } END { if (n > 0) print sum; }' diagMO.${k} >> diag.sum
done
