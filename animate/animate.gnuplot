#!/usr/local/gnuplot4.6/bin/gnuplot

set term gif animate background 'white'

# viewpoint
set view 30,90;

# colorbar option
set colorbox v user size .02,0.65;

# color palette
#set palette rgb 23,13,-31;
set palette defined ( 0.00   'red', \
                      0.333333 'white', \
                      1.00  'blue')

set xlabel "x"; set ylabel "y"; set zlabel "z"; 

# color sacle
set cbrange[-0.1:0.2];
set output "animation.gif"
set samples 100

# start point
i = 0;
load 'looper.inc'

set output;

