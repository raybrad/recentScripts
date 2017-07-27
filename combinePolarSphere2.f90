program combinePolarSphere2
implicit none
integer:: i, j, k, linum,iu, nFiles
real*8,allocatable :: tt(:),pol(:,:)
real*8:: tempPol(6)
character*30:: fileName
character*20::arg(2)

linum=10000
nFiles=10

CALL getarg(1, arg(1))
CALL getarg(2, arg(2))
read(arg(1),*)linum
read(arg(2),*)nFiles

allocate(tt(linum),pol(linum,6))
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
        read(iu,*) tt(j),(tempPol(k),k=1,6)
        pol(j,1:6)=pol(j,1:6)+tempPol(1:6)
    enddo
    close(iu)
enddo    
     pol=pol/dble(nFiles)

 open(30,file='tdPolar.sphere')    
 do j=1,linum
     write(30,'(10(E20.12,X))')tt(j),pol(j,1),pol(j,4),pol(j,6),pol(j,3),pol(j,3),pol(j,5),pol(j,5),pol(j,2),pol(j,2)
 enddo
   ! alpha_xx=polar_dq(1,)  alpha_xy=polar_dq(2,)  alpha_xz=polar_dq(3,)
   !                        alpha_yy=polar_dq(4,)  alpha_yz=polar_dq(5,)
   !                                               alpha_zz=polar_dq(6,)
 close(30)
end program combinePolarSphere2
