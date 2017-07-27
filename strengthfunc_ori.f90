program strengthfunc

implicit none
type spectra
  real::freq
  real::oscistr
end type spectra

character*10::cshape
character*30::datname
real::pi,minev,maxev,cev,increment,osc_int,temp1,temp2,dep,summ,read1,read2
integer::i,j,intno,point,num
character*10::file_uv
type(spectra),allocatable::spectrum(:)
namelist /uv_input/num,file_uv,point,dep,minev,maxev

pi=3.141592653
num=15
summ=0.0

minev=1.0
maxev=60.0
point=6000
dep=0.1
temp2=1.0d0/pi

rewind 5
read(5,uv_input,end=10)
10 continue

allocate(spectrum(num))

!cshape='line'
open(17,file=trim(file_uv))
do i=1,num
read(17,*) read1,read2
spectrum(i)%freq=read1
spectrum(i)%oscistr=read2
enddo

close(17)
	

cshape='lorentz'


datname=trim(cshape)//'_spectrum.dat'
if (trim(cshape)== 'line') then
open(25,file=datname)

do i = 1, ubound(spectrum,1)
	write(25,'(2f10.4)') spectrum(i)%freq,0.d0
	write(25,'(2f10.4)') spectrum(i)%freq, spectrum(i)%oscistr
        write(25,*)
enddo
close(25)
write(25,'(''   Created gnuplot data file  : '',A)')'line_spectrum.dat'
elseif (trim(cshape)=='lorentz') then

open(25,file=datname)
increment=(maxev-minev)/(point-1)
cev=minev
do i=1,point
   osc_int=0.0
   do j=1,ubound(spectrum,1)
	temp1=1.0d0/((cev-spectrum(j)%freq)**2+(dep)**2)
	osc_int=osc_int+spectrum(j)%oscistr*temp1
   enddo
   write(25,'(2f12.6)') 1239.84187/cev, temp2*dep*osc_int   	
   cev=cev+increment
enddo
close(25)
endif
open(26,file='lorentz.plt')
write(26,*) "set xlabel 'excite energy(nm)' "
write(26,*) "set ylabel 'strength function(1/eV)'"
write(26,*) "set xrange [200:10]"
write(26,*) "p 'lorentz_spectrum.dat' u 1:2 w l t 'strength function lorentz spectrum' "
end program strengthfunc
