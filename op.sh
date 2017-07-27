#!/bin/sh
awk -F " " '{print $3}' voltage.data > vol.dat
awk -F " " '{print $1"	"$2}' dip.data > dip.dat
line=$(cat vol.dat|wc -l)
sed -i "s@line_input=605@line_input=$line@g" fftnamelist
./ndipfft <fftnamelist

