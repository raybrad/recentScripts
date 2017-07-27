program getOcc
      implicit none

integer::i,j,k,nn,nOrbs,tmpcountera,tmpcounterb
integer,allocatable::fcountera(:),fcounterb(:)
real*8 ::tg,tb,te,tp
real*8 ::energy,occ,coef1a,coef1b,coef2a,coef2b,tempa,tempb,diffa,diffb,counterb,countera
real*8 ::tmpcoef1,tmpcoef2
real*8,allocatable::bond(:),antibond(:)
character*30::fileName1,fileName2,apd

nOrbs=2243

tb=0.0
te=740.0
tg=10.0

write(6,*) 'input total time'
read(5,*) te
nn=int((te-tb)/tg)+1
allocate(bond(nn),antibond(nn))
allocate(fcountera(nn),fcounterb(nn))
open(16,file='fcounterb.dat')
open(17,file='fcountera.dat')
do i=1,nn
    read(16,*) fcounterb(i)
    read(17,*) fcountera(i)
enddo
close(16)
close(17)

 open(22,file='bondocc.dat')
 open(23,file='antibondocc.dat')
 do i=1,nn
   tp=(i-1)*tg
   write(6,*) tp
   if(tp==0) then
       write(apd,'(I1)') int(tp)
   elseif(tp<10) then
       write(apd,'(I3)') int(tp)*100
   elseif(tp<100) then              
       write(apd,'(I4)') int(tp)*100
   elseif(tp<1000) then             
       write(apd,'(I5)') int(tp)*100
   elseif(tp<10000) then            
       write(apd,'(I6)') int(tp)*100
   endif
   fileName1='diagMO.'//trim(apd)
   write(6,*) fileName1
   open(18,file=fileName1)
   do j=1,nOrbs
       read(18,*) energy,occ
       if(j==fcounterb(i)) then
           write(22,*) tp,energy,occ
       endif
       
       if(j==fcountera(i)) then
           write(23,*) tp,energy,occ
       endif
   enddo

 enddo
    close(22)
    close(23)
end program getOcc  
