program sphereCharge
implicit none
integer::shellNat,nAt,nTime
integer:: i, j, k,tmpcounter,intv,point,flag(3),enh
real*8::time,qtot,radius,startTime
character*2,allocatable ::ele(:)
integer,allocatable::counter(:)
real*8,allocatable:: pos(:,:),qAt(:)
character*20::filename,cname

print *,'input Total num of atoms'
read(5,*) nAt
print *,'input Total time steps'
read(5,*) nTime
enh=10.d0
radius=0.8d0
intv=5
startTime=30.d0
open(16,file='shell.xyz')
read(16,*) shellNat
read(16,*)
allocate(ele(shellNat),pos(3,shellNat),counter(shellNat))
ele=''
pos=0.d0
counter=0
do k=1,shellNat
    read(16,*) ele(k),(pos(j,k),j=1,3),counter(k)
enddo
close(16)

allocate(qAt(nAt))
qAt=0.d0
qTot=0.d0
flag=0
open(20,file='deltaq_nLoc.1.data')
do i=1,nTime
    read(20,*) time,(qAt(j),j=1,nAt),qtot
      if (time >startTime .and. mod(i,intv)==0) then
        write(30,'(A10,F12.5)') '#at time:',time
            point=i/intv+1
            if(point<10) then
            write(cname,'(I1)') point
            elseif(point>=10 .and. point<100) then
            write(cname,'(I2)') point
            elseif(point>=100 .and. point<1000) then
            write(cname,'(I3)') point
            elseif(point>=1000 .and. point<10000) then
            write(cname,'(I4)') point
            else
                write(6,*) 'not dealt now'
            endif
        filename='spheres.dat'//trim(cname)
        open(30,file=filename)
        do k=1,shellNat
               tmpcounter=counter(k) 
 !              write(30,'(3(f12.5,x),E20.12)') (pos(j,k),j=1,3),(radius+qAt(tmpcounter)*10.d0)
                if(qAt(tmpcounter)>=0.d0) then
                    flag(1)=0
                    flag(2)=0
                    flag(3)=1
                else
                    flag(1)=1
                    flag(2)=0
                    flag(3)=0
                endif
                if(tmpcounter >1289) then
                    enh=20.d0
                endif
                write(30,'(3(f12.5,x),E20.12,x,3(I1,x))') (pos(j,k),j=1,3),abs(qAt(tmpcounter)*enh),flag(1),flag(2),flag(3)
        enddo
              ! write(30,*)
               close(30)
     endif
enddo
close(20)

end program sphereCharge
