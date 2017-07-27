#!/bin/sh
paste absorb.dat inc_int.dat>temp.dat
awk -F " " '{print $1"	"($2/$4)}' temp.dat > abs_cros.dat
paste scatter.dat inc_int.dat>temp.dat
awk -F " " '{print $1"	"($2/$4)}' temp.dat > scat_cros.dat
rm temp.dat

