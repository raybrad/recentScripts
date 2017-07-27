program FermiDistCDF
implicit none

integer::i,j
character*30::arg(2),FileName
real*8 ::ein,eFermi,tElec,energy,dE,nE,fermiDiracCDF,temp
real*8 ::eRange(2)
real*8, parameter :: boltz = 1.3806506d-23
real*8, parameter :: ev2j  = 1.60217646d-19

eFermi=-4.418
tElec=298.0

CALL getarg(1, arg(1))   
CALL getarg(2, arg(2))      

!input must be *.* real format or there will be error
read(arg(1),'(F12.5)') eFermi
read(arg(2),'(F12.5)') tElec

write(6,*) 'fermi energy:',eFermi
write(6,*) 'temperature:',tElec
eRange(1)=-12.0
eRange(2)= 12.0
dE=0.01
nE=int((eRange(2)-eRange(1))/dE)+1
FileName='FermiDiracCDF'//trim(arg(2))//'.dat'
open(16,file=trim(FileName))

do i=1,nE
    energy=eRange(1)+dble(i-1)*dE
!
! fermidirac = 1.d0 / ( 1.d0 + exp(ein / kT) )
!
! ein in eV, tElec in Kelvin
!
ein=energy-eFermi
  temp = dexp(ein * ev2j / (boltz * tElec))
  fermiDiracCDF= temp / (1.d0 + temp)**2
  write(16,'(2(E20.12,3x))') energy, fermiDiracCDF
enddo

close(16)

end program FermiDistCDF

