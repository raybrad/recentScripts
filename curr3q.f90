program curr2q
implicit none

integer::i,j,n
real*8::t,curr(3),temp
real*8::dt,q(3),factor,sumq(3)
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
	read(16,*) t, curr(1),curr(2),curr(3),temp
	q(1)=q(1)+curr(1)*dt*factor
	q(2)=q(2)+curr(2)*dt*factor
	q(3)=q(3)+curr(3)*dt*factor
	
	sumq(:)=sumq(:)+q(:)
	write(17,'(e12.5,2x,6(e12.5,2x))') t,q(1),q(2),q(3),sumq(1)/i,sumq(2)/i,sumq(3)/i
	!write(17,'(e12.5,2x,3(e12.5,2x))') t,q(1),q(2),q(3)
enddo
close(16)
close(17)
end


