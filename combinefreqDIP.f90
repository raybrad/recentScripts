program combinefreqDIP
implicit none
integer:: i, j, k, linum,iu, nFiles
real*8,allocatable :: freq(:),dip(:)
real*8:: tempPol
character*40:: fileName
character*20::arg(2)

linum=10000
nFiles=10

CALL getarg(1, arg(1))
CALL getarg(2, arg(2))
read(arg(1),*)linum
read(arg(2),*)nFiles

allocate(freq(linum),dip(linum))
freq=0.d0
dip=0.d0

do i=1,nFiles
    if(i<100) then 
    write(fileName,'(A5,I2.2,A12)') 'raman',i,'/freqDIP.out'
    else
    write(fileName,'(A5,I3.3,A12)') 'raman',i,'/freqDIP.out'
    endif
    write(6,*) 'processing ',trim(fileName)
    iu=10+i
    open(iu,file=trim(fileName)) 
    do j=1,linum
        read(iu,*) freq(j),tempPol
        dip(j)=dip(j)+tempPol
    enddo
    close(iu)
enddo    
     dip=dip/dble(nFiles)

 open(30,file='freqDIP.comb')    
 do j=1,linum
     write(30,'(2E20.12)') freq(j),dip(j)
 enddo
 close(30)
end program combinefreqDIP
