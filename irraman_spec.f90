program irraman_spec
implicit none
type spectra
 real::freq
 real::ir_int
 real::raman_int
end type spectra

character*10::cshape
character*30::datname,file_raman
real::pi,minwave,maxwave,cwave,increment,ir_int,raman_int,temp1,temp2,dep
integer::i,j,k,intno,point,nmode,num
type(spectra),allocatable::spectrum(:)
namelist /raman_input/nmode,cshape,point,file_raman,minwave,maxwave,dep

pi=3.141592453
nmode=1
cshape='lorentz'
minwave=0.0
maxwave=4000.0
point=8000
dep=(maxwave-minwave)*0.005
temp2=0.5d0/pi
file_raman='irraman.dat'

rewind 5
read(5,raman_input,end=10)
10 continue
allocate(spectrum(nmode))

open(20,file=file_raman)
do i=1,nmode,3
     num=min(2,nmode-i)
     read(20,*)(spectrum(i+k)%freq ,k=0,num)
     read(20,*)(spectrum(i+k)%ir_int ,k=0,num)
     read(20,*)(spectrum(i+k)%raman_int ,k=0,num)
enddo
!do i=1,nmode	
!read(20,*)spectrum(i)%freq,spectrum(i)%ir_int,spectrum(i)%raman_int
!enddo

close(20)	
datname=trim(cshape)//'_irraman_spec.dat'
if (trim(cshape)== 'line') then
open(25,file=datname)

do i = 1, nmode
	write(25,'(3f10.4)') spectrum(i)%freq,0.d0,0.d0
	write(25,'(3f10.4)') spectrum(i)%freq, spectrum(i)%ir_int,spectrum(i)%raman_int
        write(25,*)
enddo
close(25)
write(25,'(''   Created gnuplot data file  : '',A)')'line_spectrum.dat'
elseif (trim(cshape)=='lorentz') then
open(25,file=datname)
increment=(maxwave-minwave)/(point-1)
cwave=minwave
do i=1,point
    ir_int=0.0
    raman_int=0.0
   do j=1,nmode
	temp1=1.0d0/((cwave-spectrum(j)%freq)**2+(dep/2.0)**2)
	ir_int=ir_int+spectrum(j)%ir_int*temp1
	raman_int=raman_int+spectrum(j)%raman_int*temp1
   enddo
   write(25,'(3f10.4)') cwave, temp2*dep*ir_int,temp2*dep*raman_int   	
   cwave=cwave+increment
enddo
close(25)
endif
open(26,file='irraman_spec.plt')
write(26,*) 'set term pngcairo'
write(26,*) "set output 'irraman.png'"
write(26,*) "set xlabel 'wavenumber [1/cm]'" 
write(26,*) 'set xrange [0.0:4000.0]'
!write(26,*) 'set yrange [0.0:0.5]'
!write(26,*) 'set y2range [0.0:0.5]'
write(26,*) 'set y2tics'
write(26,*) "set ylabel 'Absorption coefficient [km/mol]'"
write(26,*) "set y2label 'Raman activity [A^4/amu]'"
write(26,*) "p '",trim(datname),"'  u 1:2 lc 1 lw 3  w l t 'Infrared spectrum'\"
write(26,*) ",'",trim(datname),"'  u 1:3 lc 2 lw 3  w l  axis x1y2 t 'Raman spctrum'"
end program irraman_spec
