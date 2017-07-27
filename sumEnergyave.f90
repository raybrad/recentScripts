program sumEnergy
implicit none

integer::i,j,k,nOrbs,nn,counter(5),counter2(5)
real*8 ::fermi,energy,occ
real*8 ::tg,tb,te,tp,eRange(2,5),hotenergy(5),holeenergy(5),maxm(5),temp(5)
real*8,allocatable::energy0(:),occ0(:)
character*30::fileName,outputName,apd


nOrbs=2241
fermi=-4.418

tb=0.0
te=990.0
tg=10.0

write(6,*) 'input tot time'
read(5,*)te

nn=int((te-tb)/tg)+1

do j=1,5
    eRange(1,j)=-2.71d0*dble(j)
    eRange(2,j)= 2.71d0*dble(j)
enddo

allocate(energy0(nOrbs),occ0(nOrbs))
open(16,file='diagMO.0')
do i=1,nOrbs
    read(16,*) energy0(i),occ0(i)
enddo
close(16)

open(22,file='hot.eng')
open(23,file='hole.eng')
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
   elseif(tp<1000000) then
       write(apd,'(I8)') int(tp)*100
   endif
    fileName='diagMO.'//trim(apd)
    write(6,*) fileName
    open(18,file=fileName)
    hotenergy=0.d0
    holeenergy=0.d0
    temp=0.d0
    counter=0
    counter2=0
    do j=1,nOrbs
         read(18,*) energy,occ
         if( energy-fermi>=0 .and. energy-fermi<eRange(2,1) ) then
             hotenergy(1)=hotenergy(1)+abs(occ-occ0(j))*abs(energy-fermi)
             counter(1)=counter(1)+1
             if(occ-occ0(j)>temp(1)) then
                 maxm(1)=occ-occ0(j)
             endif
         elseif(energy-fermi>=eRange(2,1) .and. energy-fermi<eRange(2,2)) then
             hotenergy(2)=hotenergy(2)+abs(occ-occ0(j))*abs(energy-fermi)
             counter(2)=counter(2)+1
             if(occ-occ0(j)>temp(2)) then
                 maxm(2)=occ-occ0(j)
             endif
         elseif(energy-fermi>=eRange(2,2) .and. energy-fermi<eRange(2,3)) then
             hotenergy(3)=hotenergy(3)+abs(occ-occ0(j))*abs(energy-fermi)
             counter(3)=counter(3)+1
             if(occ-occ0(j)>temp(3)) then
                 maxm(3)=occ-occ0(j)
             endif
         elseif(energy-fermi>=eRange(2,3) .and. energy-fermi<eRange(2,4)) then
             hotenergy(4)=hotenergy(4)+abs(occ-occ0(j))*abs(energy-fermi)
             counter(4)=counter(4)+1
             if(occ-occ0(j)>temp(4)) then
                 maxm(4)=occ-occ0(j)
             endif
         elseif(energy-fermi>=eRange(2,4) .and. energy-fermi<eRange(2,5)) then
             hotenergy(5)=hotenergy(5)+abs(occ-occ0(j))*abs(energy-fermi)
             counter(5)=counter(5)+1
             if(occ-occ0(j)>temp(5)) then
                 maxm(5)=occ-occ0(j)
             endif
         endif

         if( energy-fermi<0 .and. energy-fermi>=eRange(1,1) ) then
             holeenergy(1)=holeenergy(1)+abs(occ0(j)-occ)*abs(energy-fermi)
             counter2(1)=counter2(1)+1
         elseif( energy-fermi<eRange(1,1) .and. energy-fermi>=eRange(1,2)) then
             holeenergy(2)=holeenergy(2)+abs(occ0(j)-occ)*abs(energy-fermi)
             counter2(2)=counter2(2)+1
         elseif( energy-fermi<eRange(1,2) .and. energy-fermi>=eRange(1,3)) then
             holeenergy(3)=holeenergy(3)+abs(occ0(j)-occ)*abs(energy-fermi)
             counter2(3)=counter2(3)+1
         elseif( energy-fermi<eRange(1,3) .and. energy-fermi>=eRange(1,4)) then
             holeenergy(4)=holeenergy(4)+abs(occ0(j)-occ)*abs(energy-fermi)
             counter2(4)=counter2(4)+1
         elseif( energy-fermi<eRange(1,4) .and. energy-fermi>=eRange(1,5)) then
             holeenergy(5)=holeenergy(5)+abs(occ0(j)-occ)*abs(energy-fermi)
             counter2(5)=counter2(5)+1
         endif
    enddo
    close(18)
   ! write(22,'(6(F12.5))') tp,(hotenergy(k),k=1,5)
    !write(23,'(6(F12.5))') tp,(holeenergy(k),k=1,5)
    write(22,'(8(E12.5))') tp,(hotenergy(k)/dble(counter(k)),k=1,5),sum(hotenergy)/dble(sum(counter)),(sum(hotenergy)+sum(holeenergy))/(dble(sum(counter))+dble(sum(counter2)))
    write(23,'(7(E12.5))') tp,(holeenergy(k)/dble(counter2(k)),k=1,5),sum(holeenergy)/dble(sum(counter2))
 enddo
close(22)
close(23)
end program sumEnergy
