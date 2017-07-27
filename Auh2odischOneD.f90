program distch2
implicit none
real*8 :: xa,yb,c,deltax,deltay,deltac,Ox,Oy,Hx,Hy
character*20::arg(1)


CALL getarg(1, arg(1))
read(arg(1),*) deltac

Ox=10.47633210
Hy=9.91448764
xa=Ox-Hx
yb=Hy
c=dsqrt(xa**2.d0+yb**2.d0)

deltax=xa/c*deltac
deltay=yb/c*deltac

open(16,file='h2o.xyz2')
write(16,'(I1)') 3
write(16,*)
write(16,'(A1,3X,3(E15.8,X))') 'H', Ox-deltax,-yb-deltay,0.d0
write(16,'(A1,3X,3(E15.8,X))') 'O',10.47633210,0.00000000,0.00000000
write(16,'(A1,3X,3(E15.8,X))') 'H', 9.91448764,0.79392405,0.00000000
close(16)
end program distch2
