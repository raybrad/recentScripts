program Anglech
implicit none
real*8 :: xa,yb,c,deltax,deltay,deltac,Ox,Oy,Hx,Hy
real*8 :: theta
real*8, parameter :: pi=3.1415926535897932384626433832795029d0
character*20::arg(1)


CALL getarg(1, arg(1))
read(arg(1),*) deltac
Ox=10.47633210
Hx=9.91448764
Hy=0.79392405
xa=Ox-Hx
yb=Hy
c=dsqrt(xa**2.d0+yb**2.d0)
theta=datan(yb/xa)
theta=theta+deltac*10.d0/360.d0*2*pi
yb=c*dsin(theta)
xa=Ox-c*dcos(theta)

open(16,file='h2o.xyz')
write(16,'(A1,3X,3(E15.8,X))') 'H', 9.91448764,-0.79392405,0.00000000
write(16,'(A1,3X,3(E15.8,X))') 'O',10.47633210,0.00000000,0.00000000
write(16,'(A1,3X,3(E15.8,X))') 'H', xa,yb, 0.00000
close(16)
end program Anglech      
