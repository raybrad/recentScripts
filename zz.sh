#!/bin/sh
awk -F " " '{print $1"	"$2-$4}' temp.dat > abs_cross.dat
rm temp.dat

