#!/bin/sh
for i in `ls 011*`; do mv -f $i `echo $i | sed 's/011/11/'`;  done  
