program sphere
implicit none
integer::nAt
integer:: i, j, k,ncounter
real*8 :: temp,radius,distance
character*20::filename
character*2,allocatable ::ele(:)
integer,allocatable::counter(:)
real*8,allocatable:: pos(:,:)

print *,'input filename'
read(5,*) filename
print *, 'input radius'
read(5,*) temp
radius=temp*10.d0  !A
open(16,file=filename)
read(16,*) nAt
read(16,*)
allocate(pos(3,nAt),ele(nAt),counter(nAt))
pos=0.d0
ele=''
counter=0
ncounter=0
do i=1,nAt
read(16,*) ele(i),(pos(j,i),j=1,3)
    distance=sqrt(pos(1,i)**2+pos(2,i)**2+pos(3,i)**2)
    if (distance>=radius-3.d0) then
        ncounter=ncounter+1
        counter(ncounter)=i
    endif
enddo
close(16)
        
open(20,file='shell.xyz')
write(20,'(I4)') ncounter
write(20,*)
do i=1,ncounter
    k=counter(i)
    write(20,'(A2,x,3(f12.5,x),I4)') ele(k),(pos(j,k),j=1,3),k
enddo
close(20)
end program sphere
