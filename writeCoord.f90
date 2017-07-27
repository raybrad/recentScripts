program writeCoord
implicit none

integer::i,j,nn
real*8 ::dist,iTime
real*8,allocatable::coord1(:,:),coord2(:,:)
character::char
character*20::arg(2)

iTime=10.0
nn=126

CALL getarg(1, arg(1))
CALL getarg(2, arg(2))
read(arg(1),*)nn
read(arg(2),*)iTime
write(6,*) 'nOrbs',nn
write(6,*) 'output dt', iTime
allocate(coord1(3,nn),coord2(3,nn))
open(16,file='h1.xyz')
do i=1,nn
    read(16,*) char,(coord1(j,i),j=1,3)
enddo
close(16)

open(17,file='h2.xyz')
do i=1,nn
    read(17,*) char,(coord2(j,i),j=1,3)
enddo
close(17)

open(18,file='dist.dat')
do i=1,nn
   dist=sqrt( (coord1(1,i)-coord2(1,i))**2+(coord1(2,i)-coord2(2,i))**2+(coord1(3,i)-coord2(3,i))**2)
   write(18,'(F12.5,X,E12.5)') (i-1)*iTime, dist
enddo
close(18)
end program writeCoord
