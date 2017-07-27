program genEnField
implicit none
!
!
character*8         :: pFile
logical             :: lExist
integer             :: i, j, k, istat
integer             :: XYZ1, XYZmin, XYZmax
real*8              :: ePt, resolution, Fieldat
real*8, allocatable :: pDOS(:), Field(:)
integer		    :: direction

direction=1
rewind 5
namelist /DOS/ direction
read(5, DOS, end=111)
111 continue
!
pFile = 'deltaq_atom.dat'
inquire(file=pFile, exist=lExist)
if (.not. lExist) then
  write(*,*) ' deltaq_atom.dat not found '
  stop
end if
!
!resolution = 0.7d0 ! in angstrom
!resolution = 1.9197949d0*2.d0 ! in angstrom
if (direction==1) then
resolution = 2*0.705d0 ! in angstrom x
elseif(direction==2) then
resolution = 2*1.221d0! in angstrom y
elseif(direction==3) then
resolution = 2*2.0391d0 ! in angstrom z
endif
!resolution = 1.9197949d0 ! in angstrom
XYZmin = 100000
XYZmax =-100000
do i = 1, nAt
  XYZ1 = nint(XYZ(direction,i) / resolution * BOHR)
  XYZmin = min(XYZmin, XYZ1)
  XYZmax = max(XYZmax, XYZ1)
enddo
!
write(6,*) 'XYZmin:',XYZmin,'XYZmax:',XYZmax
allocate(pDOS(nOrbs), Field(XYZmin:XYZmax), STAT=istat)
if (istat /= 0) then
  write(6,*)"allocation error in genField.f90  "
  stop
endif
!
open(30,file=pFile)
if(direction==1) then	
open(31,file='xField.dat')
elseif(direction==2) then	
open(31,file='yField.dat')
elseif(direction==3) then	
open(31,file='zField.dat')
endif

read: do
  read(30,'(999999ES20.12)',iostat=istat) ePt, (pDOS(j), j=1,nOrbs)
  readstatus: if (istat == 0) then
    Field(:)= 0d0
!   Fieldat = 0d0
    do i = 1, nAt
      XYZ1 = nint(XYZ(direction,i) / resolution * BOHR)
!     if (XYZ(1,i) * BOHR > 60d0) cycle
!     if (XYZ(1,i) * BOHR <-60d0) cycle
      do j = ind(i)+1, ind(i+1)
        Field(XYZ1) = Field(XYZ1) + pDOS(j)
!       Fieldat = Fieldat + pDOS(j)
      enddo
    enddo
    do i = XYZmin, XYZmax
      if (dabs(Field(i)) < 1d-9) cycle
      write(31,'(f10.5, 2X, 2(e15.7,2X))') dble(i)*resolution, ePt, Field(i)
    enddo
    write(31,*)
!   write(31,*) ePt, Fieldat
  else readstatus
    exit
  end if readstatus
enddo read
close(30)
close(31)
!
end subroutine genEnField
