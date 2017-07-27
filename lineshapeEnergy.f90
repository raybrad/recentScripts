program lineshapeEnergy
  implicit none
integer::i,j,nOrbs
real*8 :: fermi
character*20::arg(2),FileName
real*8,allocatable ::energy0(:),energy(:),occ0(:),occ(:)      


CALL getarg(1, arg(1))      
CALL getarg(2, arg(2))      

nOrbs=2241
fermi=-4.418

read(arg(1),'(I4)') nOrbs
! write(6,*) 'input nOrbs'
! read(5,*) nOrbs

allocate(energy0(nOrbs),energy(nOrbs),occ0(nOrbs),occ(nOrbs))      

FileName='./diagMO.0'
open(15,file=trim(FileName))      
do i=1,nOrbs
     read(15,*) energy0(i),occ0(i)
enddo      
close(15)

FileName='./diagMO.'//trim(arg(2))
open(16,file=trim(FileName))      
do i=1,nOrbs
     read(16,*) energy(i),occ(i)
enddo
close(16)


FileName=trim(arg(2))//'.eline'
open(17,file=trim(FileName))
do i=1,nOrbs
     write(17,*) energy0(i),energy(i), 0.d0
     write(17,*) energy0(i),energy(i),abs(occ(i)-occ0(i))*abs(energy0(i)-fermi)
     write(17,*)
enddo
close(17)


end program lineshapeEnergy
