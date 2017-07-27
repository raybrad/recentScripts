program Anglech
implicit none
real*8 :: xa,yb,c,deltax,deltay,deltac
real*8 :: theta
real*8, parameter :: pi=3.1415926535897932384626433832795029d0
character*20::arg(1)


CALL getarg(1, arg(1))
read(arg(1),*) deltac
xa=0.57401577d0
yb=0.77848691d0
c=dsqrt(xa**2.d0+yb**2.d0)
theta=datan(yb/xa)
theta=theta+deltac*10.d0/360.d0*2*pi
yb=c*dsin(theta)
xa=c*dcos(theta)

open(16,file='h2o.xyz')
write(16,'(I1)') 3
write(16,*)
write(16,'(A1,3X,3(E12.5,X))') 'O',  0.d0, 0.d0, 0.d0
write(16,'(A1,3X,3(E12.5,X))') 'H', -xa,  yb,0.d0
write(16,'(A1,3X,3(E12.5,X))') 'H', -0.57402E+00,-0.77849E+00, 0.00000E+00
close(16)
end program Anglech      
