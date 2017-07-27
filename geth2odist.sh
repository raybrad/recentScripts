grep -iw 'H' tdCoord.data >temp 
grep -iw 'O' tdCoord.data >o.xyz 
sed -n '1~2p' temp>h1.xyz
sed -n '2~2p' temp>h2.xyz
rm temp
line=$(wc -l <h1.xyz)

~/program/gadgets/getH2Odist  $line $1
