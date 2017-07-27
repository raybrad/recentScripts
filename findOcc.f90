program findOcc
 implicit none
integer::i,j,k,nn,nOrbs,counter
real*8 ::tg,tb,te,tp
real*8 ::energy,occ,targetE
real*8 ::temp,diff,tempE,tempOcc
real*8,allocatable::recEnergy(:),recOcc(:)
character*30::fileName,outputName,apd


write(6,*) 'input target energy'
read(5,*) targetE
nOrbs=2243

tb=0.0
te=995.0
tg=5.0

nn=int((te-tb)/tg)+1

allocate(recEnergy(nn),recOcc(nn))
write(apd,'(I1)') int(targetE)
outputName=trim(apd)//'ev.time'
open(22,file=outputName)
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
   elseif(tp<100000) then            
       write(apd,'(I7)') int(tp)*100
   endif
    fileName='diagMO.'//trim(apd)
    write(6,*) fileName
    open(18,file=fileName)
    temp=1.d3
    counter=0
    do j=1,nOrbs
         read(18,*) energy,occ
         diff=abs(targetE-energy)
         if(diff<temp) then
             temp=diff
             counter=j
             tempE=energy
             tempOcc=occ
         endif
    enddo
    recEnergy(i)=tempE
    recOcc(i)=occ
    write(22,*) tp,recEnergy(i),recOcc(i)

    write(6,*)'pos',counter
    close(18)

 enddo
close(22)
end program findOcc
