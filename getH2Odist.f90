program pyscale
implicit none

integer::i,j,nn
real*8 ::disth1,disth2,iTime
real*8,allocatable::freq(:),intensity(:)
character*20::FileName
character*20::arg(2)

nn=8000

CALL getarg(1, arg(1))
read(arg(1),*)nn
read(arg(2),*)FileName
write(6,*) 'lines',nn
write(6,*) 'FileName',FileName
allocate(freq(nn),intensity(nn))
open(16,file=trim(FileName))
do i=1,nn
    read(16,*) freq(i),intensity(i)
enddo
close(16)

open(19,file='lorentz_spectrum.dat')
do i=1,nn
    if(freq(i) <1400) then
        freq(i)=freq(i) * 0.91
    elseif(freq(i) >=1400) then
        freq(i)=freq(i) * 1.03
    endif
    write(19,'(E20.12,2X,E20.12)') freq(i),intensity(i)
enddo
close(19)
end program pyscale
