inputFile=dftb_in.hsd
bandFile=band_tot.dat
detailFile=detailed.out
eBot=5
eTop=5
eBot2=7
eTop2=7

eFermi=`grep -n "Fermi energy" $detailFile |awk '{print $6}' `
ebb=`echo "$eFermi-$eBot"|bc`
ett=`echo "$eFermi+$eTop"|bc`
echo "Ef=$eFermi"
echo "unset key"
echo "set yr [-$eBot2:$eTop2]"
lineNum=`grep -n Klines $inputFile | cut -f1 -d:`
lineNum=$((lineNum+1))
npoint=`tail -n +$lineNum dftb_in.hsd |awk '{print $1}'` 
label=`tail -n +$lineNum dftb_in.hsd |awk '{print $6}'`
k=0
xtics=0
echo -n "set xtics("
for ll in $label
do
    k=$((k+1))
	npt=`echo $npoint|cut -f $k -d ' '`
	xtics=$((xtics+npt))
	if [ $k -ne 1 ]; then
		echo -n " , "
	fi
	echo -n "\"$ll\" $xtics"
done
echo ")"
echo "set xr [1:$xtics]"
k=0
xtics=0
for ll in $label
do
    k=$((k+1))
	npt=`echo $npoint|cut -f $k -d ' '`
	xtics=$((xtics+npt))
	echo "set arrow from $xtics,graph(0,0) to $xtics,graph(1,1) nohead"
done

#nLevel=`grep KPT band.out |cut -f 1 -d:|sed -n 2p`
#nLevel=$((nLevel-2))
nLevel=`awk  'NR==1{print NF;exit}' $bandFile`
echo "p 0.0 w l lc 5,for [col=2:$nLevel] '$bandFile' using 1:(column(col)-Ef)  lw 2 lc 1 w l"
