program reducenoise
implicit none
integer::n,i,j,k
real*8 ::tt, diff, tempPol(9),lastPol(9)

n=700
tempPol=0.d0
lastPol=0.d0
!open(16,file='tdPolar.data')      
open(16,file='tdPolar.sphere')      
open(17,file='tdPolar.red')      
do j=1,n
    !read(16,'(10(E20.12,X))') tt, (tempPol(k),k=1,9)
    read(16,*) tt, (tempPol(k),k=1,9)
    if(j>0.5) then
        do k=1,6
           diff=abs((tempPol(k)-lastPol(k))/tempPol(k))
           if (diff > 5.0) then
               tempPol(k) = lastPol(k)
           endif
        enddo
    endif
    lastPol = tempPol
     write(17,'(10(E20.12,X))') tt,(tempPol(k),k=1,9)
enddo
close(16)
close(17)

end program reducenoise      
