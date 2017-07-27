program checkCoord
implicit none
integer:: i, j, k, nAt,ncount
real*8 :: tt, dt,radius,origin(3),temp
character:: tempChar,tempChar1
real*8,allocatable::coord0(:,:), coord(:,:)

tt=10.d0
dt=0.1d0
write(6,*) 'input total time:'
read(5,*) tt

write(6,*) 'input dt:'
read(5,*) dt


ncount=int(tt/dt)+1

open(15,file='radius.dat')
open(16,file='tdCoord.data')
      ! read(16,*) nAt
      ! allocate(coord0(3,nAt),coord(3,nAt))
      ! read(16,*) tempChar
      ! do j=1,nAt
      !   read(16,*) tempChar1,(coord0(k,j),k=1,3)
      ! enddo
      ! write(15,*) 0.d0,0.d0
do i=1,ncount
      origin=0.d0
      radius=0.d0
      read(16,*) nAt
      if(i==1) then
          allocate(coord0(3,nAt),coord(3,nAt))
      endif
      read(16,*) tempChar
      do j=1,nAt
        read(16,*) tempChar1,(coord(k,j),k=1,3)
        origin(1:3)=origin(1:3)+coord(1:3,j)
      enddo
        origin=origin/dble(nat) 
       do j=1,nAt
         temp=dsqrt((coord(1,j)-origin(1))**2+(coord(2,j)-origin(2))**2+(coord(3,j)-origin(3))**2)
         if(temp>radius) then
             radius=temp
         endif
       enddo    
      
       write(15,*) i*dt,radius
enddo

end program checkCoord
