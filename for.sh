#!/bin/bash
#-----------------------------------------------
# run MonteCarlo simulation with different input
#-----------------------------------------------
inputDIR="./input"
resultDIR="./result"
exeDIR="../"
exe=DoubleWell
ntimes=100

for i in 0.1 1.0 10
do

  for j in $(seq 1 $ntimes)
do
#echo "copying input files of RefT$i"
#cp $inputDIR/PotB0-RefT$i.in  .
echo "running: ../DoubleWell <PotB0-RefT$i.in >PotB0-RefT$i.out"
$exeDIR$exe< $inputDIR/PotB0-RefT$i.in >$resultDIR/PotB0-RefT$i.out$j

echo "grep MeanCoord & MeanPotential"
grep  -iw "MeanCoord" $resultDIR/PotB0-RefT$i.out$j |cat >>$resultDIR/MeanCoord$i.dat
grep  -iw "MeanPotential" $resultDIR/PotB0-RefT$i.out$j |cat >>$resultDIR/MeanPotential$i.dat

echo "rm output data"
rm $resultDIR/PotB0-RefT$i.out$j
echo "finish successfully"

done
done
