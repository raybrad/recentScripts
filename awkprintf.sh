awk '{ printf("% 13.11e  % 13.11e \n", $1, $2*0.001) }' inc_Ex.dat > inc_Ex2.dat
