program fitFermi
implicit none

integer::i,j
character*20::arg(3),FileName
integer:: n
real*8 ::ein,eFermi,occ,x,b,a,y,r,temp
real*8 :: sumx,sumy,sumx2,sumy2,sumxy,deno
real*8,allocatable :: energy0(:),energy(:)
real*8, parameter :: boltz = 1.3806506d-23
real*8, parameter :: ev2j  = 1.60217646d-19

eFermi=-4.418
n=1

CALL getarg(1, arg(1))      
CALL getarg(2, arg(2))      
CALL getarg(3, arg(3))      

!input must be *.* real format or there will be error
read(arg(1),'(F12.5)') eFermi
read(arg(2),'(I5)') n
read(arg(3),*) FileName 

write(6,*) 'fermi energy:',eFermi
write(6,*) 'line No. :',n
write(6,*) 'processing ', FileName

allocate(energy0(n),energy(n))
energy0=0.d0
energy =0.d0
open(14,file='diagMO.base')
do i=1,n
read(14,*) energy0(i),temp
enddo

open(15, file ='temp.dat')
open(16,file=trim(FileName))
sumx  =0.d0   
sumy  =0.d0
sumx2 =0.d0
sumy2 =0.d0
sumxy =0.d0
do i=1,n
    read(16,*) energy(i), occ
    if(occ>2.0) then
         write(6,*) 'occ approx' , occ
!        cycle
         occ=2.d0-1.d-10
    elseif(occ<0.d0) then
         write(6,*) 'occ approx' , occ
         occ=1.d-10
         !cycle
    endif
!
! fermidirac = 1.d0 / ( 1.d0 + exp(ein / kT) )
!
! ein in eV, tElec in Kelvin
!
!y=ax+b
!energy = T * ln(2/f-1)*k/ev2j + e_f
  y = energy0(i)
  x = log(2.d0/occ-1.0)*boltz/ev2j
  b = eFermi
  sumx =sumx+x
  sumy =sumy+y
  sumx2=sumx2+x*x
  sumy2= sumy2 + y * y 
  sumxy=sumxy+x*y
  write(15,'(2(E20.12,3X))') x , y 
enddo
close(15)
close(16)
deno=  n*sumx2-sumx*sumx
a   = (n*sumxy-sumx*sumy)/deno
b   = (sumx2*sumy-sumx*sumxy)/deno
r   = (sumxy - sumx * sumy / n) /sqrt((sumx2 - sumx**2/n) * (sumy2 - sumy**2/n))

write(6,*) 'slope a', a 
write(6,*) 'y-intercept b', b 
write(6,*) 'Correlation  r = ', r
! open(17,file='fitFermi.dat')
! write(17,'(3(E20.12,3X))') x , y , b
! close(17)

end program fitFermi

