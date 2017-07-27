program combineDIP
implicit none
integer:: i, j, k, linum,iu, nFiles
real*8,allocatable :: tt(:),dip(:,:)
real*8:: tempPol(3)
character*40:: fileName
character*20::arg(2)

linum=10000
nFiles=10

CALL getarg(1, arg(1))
CALL getarg(2, arg(2))
read(arg(1),*)linum
read(arg(2),*)nFiles

allocate(tt(linum),dip(linum,3))
tt=0.d0
dip=0.d0

do i=1,nFiles
    if(i<100) then 
    write(fileName,'(A5,I2.2,A16)') 'raman',i,'/dip_nLoc.1.data'
    else
    write(fileName,'(A5,I3.3,A16)') 'raman',i,'/dip_nLoc.1.data'
    endif
    write(6,*) 'processing ',trim(fileName)
    iu=10+i
    open(iu,file=trim(fileName)) 
    do j=1,linum
        read(iu,*) tt(j),(tempPol(k),k=1,3)
        dip(j,1:3)=dip(j,1:3)+tempPol(1:3)
    enddo
    close(iu)
enddo    
     dip=dip/dble(nFiles)

 open(30,file='dip_nLoc.1.data')    
 do j=1,linum
     write(30,'(7E20.12)') tt(j),(dip(j,k),k=1,3)
 enddo
 close(30)
end program combineDIP
