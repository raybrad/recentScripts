program curr2q
implicit none

integer::i,j,n
real*8::t,dt,deltaq,sumdq
real*8, parameter :: el_charge=1.602176462d-19

n=1000
dt=0.01d0
write(6,*) "Please input the line number:"
read(5,*) n
write(6,*) "Please input the dt:"
read(5,*) dt

t=0.d0
deltaq=0.d0
sumdq=0.d0


open(17,file='deltaq.data')
open(18,file='aveq.data')
do i=1,n
	read(17,*) t, deltaq

	sumdq=sumdq+deltaq
	write(18,'(e12.5,2x,e12.5)') t,sumdq/dble(i)
enddo
close(17)
close(18)
end


