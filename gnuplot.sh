#!/bin/sh
for ((i=0; i<=10; i++))
do
    var[$i] = i*10
done

gnuplot<<EOF
values="${var[*]}"
do for [j in values] {
val=sprintf("%d", j) ##access part of the array directly
}
EOF

#need to access ${var[j]} somehow
