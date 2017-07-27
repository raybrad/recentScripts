program distch2
implicit none
integer:: i,j,nAt
real*8 :: xa,yb,c,deltax,deltay,deltac,Ox,Oy,Hx,Hy
character*20::arg(1)
character*1,allocatable:: ele(:)
real*8,allocatable     :: xyz0(:,:)


CALL getarg(1, arg(1))
read(arg(1),*) deltac

nAt=3 
allocate(ele(nAt),xyz0(nAt,3))
open(15, file='h2o.xyz')
do i=1,nAt
read(15,*) ele(i),(xyz0(i,j),j=1,3)
enddo
Ox=xyz0(2,1)
Hx=xyz0(3,1)
Hy=xyz0(3,2)
xa=Hx-Ox
yb=Hy
write(6,*) 'xa',xa,'yb',yb
c=dsqrt(xa**2.d0+yb**2.d0)

deltax=xa/c*deltac
deltay=yb/c*deltac

xyz0(1,1)= Hx+deltax
xyz0(1,2)=-yb-deltay

open(16,file='h2o.xyz2')
do i=1,nAt
write(16,'(A1,3X,3(E15.8,X))') ele(i),(xyz0(i,j),j=1,3)
enddo
close(16)
end program distch2
