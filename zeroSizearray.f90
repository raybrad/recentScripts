program main
      implicit none
      integer:: dimen
      real*8,allocatable::mat
      do dimen=0,3
      call zeroSizeArray(dimen,mat)
      enddo
end program main
subroutine zeroSizeArray(dimen,mat)
      implicit none
      integer,intent(in)::dimen
      real*8, intent(in)::mat(dimen)
      
      integer		:: length

      length=size(mat)

      write(6,*) 'length of mat :',length
end subroutine zeroSizeArray
      

