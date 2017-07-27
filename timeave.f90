program timeave
implicit none
integer::i,j,k,ia,ib,tp,nOrbs,counter
real*8 :: tmp1,tmp2
real*8,allocatable ::energy(:),occ(:)
character*30::apd,fileName

nOrbs=2241
allocate(energy(nOrbs),occ(nOrbs))
energy=0.d0
occ=0.d0
counter=0
do i=0,17
    ia=100*i+10
    ib=100*(i+1)
    ! energy=0.d0
    ! occ=0.d0
    ! counter=0
    do tp=ia,ib,10
        counter=counter+1
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
        fileName='../diagMO.'//trim(apd)
        write(6,*) fileName
        open(18,file=fileName)
        do j=1,nOrbs
             read(18,*) tmp1,tmp2
             energy(j)=energy(j)+tmp1
             occ(j)   = occ(j)  +tmp2
        enddo
        close(18)
        ! if(tp==50) then
        !     open(19,file='5000.ave')
        !     do j=1,nOrbs
        !         write(19,*) energy(j)/dble(counter),occ(j)/dble(counter)
        !     enddo 
        ! endif
     enddo
     fileName=trim(apd)//'.ave'
     open(19,file=fileName)
     do j=1,nOrbs
         write(19,*) energy(j)/dble(counter),occ(j)/dble(counter)
     enddo 
     close(19)
enddo


end program timeave
