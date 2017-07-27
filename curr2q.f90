program curr2q
implicit none

integer::i,j,n
real*8::t,curr(2),temp
real*8::dt,q(2),factor,sumq(2)
real*8, parameter :: el_charge=1.602176462d-19

n=1000
dt=0.01d0
write(6,*) "Please input the line number:"
read(5,*) n
write(6,*) "Please input the dt:"
read(5,*) dt

t=0.d0
curr=0.d0
q=0.d0
sumq=0.d0
temp=0.d0

factor=1.d-9*1.d-15/el_charge

open(16,file='curr.data')
open(17,file='qtot.dat')
do i=1,n
	read(16,*) t, curr(1),curr(2)
	q(1)=q(1)+curr(1)*dt*factor
	q(2)=q(2)+curr(2)*dt*factor

	sumq(:)=sumq(:)+q(:)
	write(17,'(e12.5,2x,4(e12.5,2x))') t,q(1),q(2),sumq(1),sumq(2)
	!write(17,'(e12.5,2x,2(e12.5,2x))') t,q(1),q(2)
enddo
close(16)
close(17)
end


