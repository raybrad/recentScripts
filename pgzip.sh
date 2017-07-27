#!/bin/sh
sourcepath=`pwd`
targetpath=$sourcepath
for i in `ls ${sourcepath}/*`
do
tar xvf $i -C $targetpath
done
