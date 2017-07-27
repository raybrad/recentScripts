program trick
      implicit none
      integer::i,j,k,nOrbs,ntime
      character*30::timeName,fileName
      logical::lExist

      integer,allocatable:: eigInd(:)
      real*8,allocatable :: eigVal(:),occ(:),temp(:)

      nOrbs=2241
      allocate(eigInd(nOrbs),eigVal(nOrbs),occ(nOrbs),temp(nOrbs))
      open(16,file='diagMO.base')
      do i=1,nOrbs
          read(16,*) eigVal(i),occ(i),eigInd(i)
      enddo
      close(16)
      
      ntime=14126
      do k=0,ntime,2

            if(k <10) then
                write(timeName,'(I1)') k
            elseif(k<100) then
                write(timeName,'(I2)') k
            elseif(k<1000) then
                write(timeName,'(I3)') k
            elseif(k<10000) then
                write(timeName,'(I4)') k
            elseif(k<100000) then
                write(timeName,'(I5)') k
            endif
            filename='diagMO.'//trim(timeName)
            write(6,*) filename
            inquire(file=trim(fileName), exist=lExist)
            if ( (lExist)) then
                  open(17,file=trim(fileName))
                  do i=1,nOrbs
                      read(17,*) eigVal(i),occ(i)
                  enddo
                  close(17)

                  do i=1,nOrbs
                    temp(eigInd(i))=occ(i)
                  enddo
                    occ=temp

                  open(18,file=trim(fileName))
                  do i=1,nOrbs
                      write(18,'(2(f20.12,x))') eigVal(i),occ(i)
                  enddo
                  close(18)
            end if
      enddo
end program trick
