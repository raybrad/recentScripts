awk 'BEGIN {max = 0;linum=0} {if ($1>max) {max=$1; linum=NR}} END {print "Max=", max ,"linenum",linum}'

awk 'BEGIN {min =10000;linum=0} {if ($1<min) {min=$1; linum=NR}} END {print "Min=", min ,"linenum",linum}'
