grep -iw 'H' tdCoord.data >coord.dat 
sed -n '1~2p' coord.dat>h1.xyz
sed -n '2~2p' coord.dat>h2.xyz

line=$(wc -l <h1.xyz)

~/program/gadgets/writeCoord  $line $1
