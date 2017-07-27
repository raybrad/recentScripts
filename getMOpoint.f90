program getMOpoint
      implicit none

integer::i,j,k,nn,nOrbs,tmpcountera,tmpcounterb
integer,allocatable::fcountera(:),fcounterb(:)
real*8 ::tg,tb,te,tp
real*8 ::energy,occ,coef1a,coef1b,coef2a,coef2b,tempa,tempb,diffa,diffb,counterb,countera
real*8 ::tmpcoef1,tmpcoef2
real*8,allocatable::bond(:),antibond(:)
character*30::fileName1,fileName2,dirName

nOrbs=2243

tb=0.0
te=730.0
tg=10.0

! write(6,*) 'input total time'
! read(5,*) te

nn=int((te-tb)/tg)+1
allocate(bond(nn),antibond(nn))
allocate(fcountera(nn),fcounterb(nn))
bond=0.d0
antibond=0.d0
open(16,file='bond.dat')
do i=1,nn
read(16,*)bond(i)
enddo
close(16)

open(17,file='antibond.dat')
do i=1,nn
read(17,*)antibond(i)
enddo
close(17)

 open(20,file='bondtime.dat')
 open(21,file='antibondtime.dat')
 do i=1,nn
    tp=(i-1)*tg
    write(6,*) tp
    if(tp<10) then
        write(dirName,'(I1)') int(tp)
    elseif(tp<100) then
        write(dirName,'(I2)') int(tp)
    elseif(tp<1000) then
        write(dirName,'(I3)') int(tp)
    elseif(tp<10000) then
        write(dirName,'(I4)') int(tp)
    endif
    fileName1=trim(dirName)//'fs/lode/diagMO.0'
    write(6,*) fileName1
    open(18,file=fileName1)
    tempb=1.d3
    tempa=1.d3
    counterb=0
    countera=0
    do j=1,nOrbs
         read(18,*) energy,occ
         diffb=abs(bond(i)-energy)
         diffa=abs(antibond(i)-energy)
         if(diffb<tempb) then
             tempb=diffb
             counterb=j
         endif
         if(diffa<tempa) then
             tempa=diffa
             countera=j
         endif
     enddo
     write(6,*)'pos', counterb,countera
     close(18)
     fileName2=trim(dirName)//'fs/lode/MOcoeff.dat'
     write(6,*) fileName2
     open(19,file=fileName2)
     tmpcounterb=0
     tmpcountera=0
     coef1a=0.0
     coef1b=0.0
     coef2a=0.0
     coef2b=0.0
     do k=1,nOrbs
          read(19,*) tmpcoef1,tmpcoef2
          if(abs(k-counterb)<=10) then
              tmpcounterb=tmpcounterb+1
             if(abs(coef1b)<abs(tmpcoef1)) then
                coef1b=tmpcoef1
                coef2b=tmpcoef2
                fcounterb(i)=k
             endif
             if(tmpcounterb==21) then
              write(20,*) tp, coef1b,coef2b
              write(6,*) 'fcounterb',fcounterb(i)
             endif
          endif
          
          if(abs(k-countera)<=10) then
              tmpcountera=tmpcountera+1
             if(abs(coef1a)<abs(tmpcoef1)) then
                coef1a=tmpcoef1
                coef2a=tmpcoef2
                fcountera(i)=k
             endif
             if(tmpcountera==21) then
                write(21,*) tp, coef1a,coef2a
                write(6,*) 'fcountera',fcountera(i)
             endif
          endif
   enddo
   close(19)

 enddo
    close(20)
    close(21)
end program getMOpoint  
