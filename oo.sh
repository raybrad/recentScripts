paste time.dat 0.20E2.dat>>temp.dat
awk -F " " '{print $1"  "sqrt($2*5.291772e-11)}' temp.dat>0.20timeE2.dat
rm temp.dat
