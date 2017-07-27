program distch2
implicit none
real*8 :: xa,yb,c,deltax,deltay,deltac
character*20::arg(1)


CALL getarg(1, arg(1))
read(arg(1),*) deltac
xa=0.57401577d0
yb=0.77848691d0
c=dsqrt(xa**2.d0+yb**2.d0)

deltax=xa/c*deltac
deltay=yb/c*deltac

open(16,file='h2o.xyz')
write(16,'(I1)') 3
write(16,*)
write(16,'(A1,3X,3(E12.5,X))') 'O',  0.d0, 0.d0, 0.d0
write(16,'(A1,3X,3(E12.5,X))') 'H', -xa, yb,0.d0
write(16,'(A1,3X,3(E12.5,X))') 'H', -xa-deltax,-yb-deltay,0.d0
close(16)
end program distch2
