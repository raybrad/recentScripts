program combinePolarSphere
implicit none
integer:: i, j, k, linum,iu, nFiles
real*8,allocatable :: tt(:),pol(:,:)
real*8:: tempPol(9)
character*30:: fileName
character*20::arg(2)

linum=10000
nFiles=10

CALL getarg(1, arg(1))
CALL getarg(2, arg(2))
read(arg(1),*)linum
read(arg(2),*)nFiles

allocate(tt(linum),pol(linum,9))
tt=0.d0
pol=0.d0

do i=1,nFiles
    if(i<100) then 
    write(fileName,'(A5,I2.2,A13)') 'raman',i,'/tdPolar.data'
    else
    write(fileName,'(A5,I3.3,A13)') 'raman',i,'/tdPolar.data'
    endif
    write(6,*) 'processing ',trim(fileName)
    iu=10+i
    open(iu,file=trim(fileName)) 
    do j=1,linum
        read(iu,*) tt(j),(tempPol(k),k=1,9)
        pol(j,1:9)=pol(j,1:9)+tempPol(1:9)
    enddo
    close(iu)
enddo    
     pol=pol/dble(nFiles)

 open(30,file='tdPolar.sphere')    
 do j=1,linum
     write(30,'(10E20.12)') tt(j),(pol(j,k),k=1,9)
 enddo
 close(30)
end program combinePolarSphere
