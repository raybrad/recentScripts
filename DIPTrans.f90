program DIPTrans
implicit none
!
real*8, parameter :: pi=3.1415926535897932384626433832795029d0
real*8, parameter :: el_mass=9.10938188d-31
real*8, parameter :: planck_h=6.62606876d-34
real*8, parameter :: mu_0=4*pi*1.d-07
real*8, parameter :: speedoflight=2.99792458d08
real*8, parameter :: el_charge=1.602176462d-19
real*8, parameter :: planck_hquer=planck_h/(2*pi)            ! 1.054571596d-34
real*8, parameter :: epsilon_0=1/(speedoflight*speedoflight*mu_0) ! 8.85419d-12
real*8, parameter :: bohr_m=4.0d0*pi*epsilon_0*planck_hquer*planck_hquer/(el_mass*el_charge*el_charge) ! 5.291772083d-11
real*8, parameter :: factor=1.0/(4*pi*epsilon_0*bohr_m/el_charge)!e/V/bohr     !alpha(e*bohr/V) -> alpha(bohr^3)
real*8, parameter :: hbar=0.658211928d0 !eV*fs
real*8, parameter :: avogadro=6.02214199d23
real*8, parameter :: atomic_mass_unit=1/avogadro*1d-03       ! 1.66056d-27 kg =1 a.m.u
real*8, parameter :: boltz   = 1.3806503d-23
real*8, parameter :: conver=2.d0*planck_h*pi*pi/speedoflight*1.d-40/atomic_mass_unit
real*8, parameter :: vbyt=planck_h*speedoflight*100.0d0/boltz
real*8, parameter :: ev2j  = 1.60217646d-19
    !----------------------------------------------------------------
    ! note: lambda_... is a wavenumber
    !----------------------------------------------------------------
    ! vbyt = planck_h*speedoflight*10^2(cm/m) / boltz
    !----------------------------------------------------------------
complex*16 ::s(3)
complex*16,allocatable :: dept_dip(:,:)
complex*16 :: meanV(3)
real*8,allocatable:: tt(:)
real*8::nmomega,omega,omega_start,omega_end,domega,dep 
real*8:: dtt,tt_min,tt_max
real*8:: time,t_dip(3)
real*8:: lambda_z,temperature,qmfactor,tempfactor
real*8:: dipol_dq2mks,prefactor,ir_int
integer*4:: i,j,k,line_output
integer*4:: line_input
real*8:: win_blackman_harris,norm
real*8:: wsum
real*8:: win_Gauss,win_hanning,win_FlatTop
COMMON /glob/ wsum

namelist /freq/omega_start,omega_end,line_input,line_output,dep,temperature

wsum=0.d0
temperature=298.0
rewind 5
read(5,freq,end=102)
102 continue 
PRINT *,'starting freq. in eV',omega_start
PRINT *,'end freq. in eV',omega_end
PRINT *,'number of lines in the file',line_input
PRINT *,'number of output points',line_output
PRINT *,'dephasing',dep

dipol_dq2mks= el_charge/sqrt(atomic_mass_unit)
prefactor = 1/(4*pi*epsilon_0)*(avogadro*pi)/(3.0d0*speedoflight**2)
prefactor = 1.0d-3 * prefactor
!
!
dep=-1d0*dep/hbar
!
allocate(dept_dip(line_input,3),tt(line_input))
!
meanV=(0.0,0.0)
OPEN(60,FILE='dip_nLoc.1.data',STATUS='OLD')
do i=1,line_input
 read(60,*)time,(t_dip(j),j=1,3)
  tt(i)=time
   dept_dip(i,1:3)=(1.,0.)*t_dip(1:3)
   meanV(1:3)=meanV(1:3) + dept_dip(i,1:3)
enddo
close(60)


meanV=meanV/dble(line_input)
!zero
do i=1,line_input
   dept_dip(i,1:3)=(dept_dip(i,1:3)-meanV(1:3))
enddo
!detrend
! do i=1,3
!     call detrend(line_input,dept_dip(:,i))
! enddo
!window
do i=1,line_input
   !dept_dip(i,1:3)=(dept_dip(i,1:3))*dexp(dep*tt(i))
    dept_dip(i,1:3)=(dept_dip(i,1:3))*win_blackman_harris(i,line_input)
   !dept_dip(i,1:3)=(dept_dip(i,1:3))*win_Gauss(i,line_input)
   !dept_dip(i,1:3)=(dept_dip(i,1:3))*win_hanning(i,line_input)
   !dept_dip(i,1:3)=(dept_dip(i,1:3))*win_FlatTop(i,line_input)
enddo

!norm = sqrt(2.0/dble(wsum*line_input));


!
OPEN(72,FILE='freqDIP.out',STATUS='Replace')
OPEN(73,FILE='freqDIP.s',STATUS='Replace')
tt_min=tt(1)
tt_max=tt(line_input)
dtt=tt(2)-tt(1)
domega=(omega_end-omega_start)/(line_output-1.d0)
omega=omega_start-domega
!
do i=1,line_output
  omega=omega+domega
  !new
  !omega=0.d0
  s=(0.d0,0.d0)
  do j=1,line_input-1
     dtt=tt(j+1)-tt(j)
      do k=1,3
    	 s(k)=s(k)+exp(dcmplx(0.d0,omega*tt(j)/hbar))*dept_dip(j,k)*dtt
      enddo
  enddo
  s(k)=s(k)!*norm
  omega=omega
  nmomega=1239.84187/omega
  
  lambda_z=omega*1.d7/1239.84187
  ! if(lambda_z < 400.d0) then
  !       s=(0.d0,0.d0)
  ! endif
  tempfactor = ev2j/(boltz*temperature)*omega
!  qmfactor   = tempfactor/(1.d0-dexp(-tempfactor))
   qmfactor=omega**2*tempfactor/(1.d0-dexp(-tempfactor))
  ir_int= prefactor*dipol_dq2mks**2*(abs(s(1))**2+abs(s(2))**2+abs(s(3))**2)*qmfactor
  WRITE(72,'(2e24.14)')lambda_z,ir_int
  WRITE(73,'(4e24.14)')lambda_z,abs(s(1)),abs(s(2)),abs(s(3))
enddo
close(72)
close(73)

!     auf debye/(angstrom*amu)
!maxir=0.0d0
!do i=1,ubound(spectrum,1)
!  if(spectrum(i)%ir_int.gt.maxir) &
!   & maxir = spectrum(i)%ir_int
!end do
!if(maxir.gt.1.d-8) then
!  do i=1,ubound(spectrum,1)
!    spectrum(i)%rel_ir_int=spectrum(i)%ir_int/maxir
!  end do
!else
!  do i=1,ubound(spectrum,1)
!    spectrum(i)%rel_ir_int=0.0d0
!  end do
!end if

!Blackman–Harris window
!dM=M+1
! for(j=0; j<M/2; j++)
!     {
!      WinCoeff[j] = 0.35875
!      - 0.48829 * cos((double)(j+1) * M_2PI / dM)
!      + 0.14128 * cos((double)(j+1) * M_2PI * 2.0 / dM)
!      - 0.01168 * cos((double)(j+1) * M_2PI * 3.0 / dM);
!     }
end program DIPTrans

function win_blackman_harris(j, n)
integer,intent(in)::j,n
real*8:: win_blackman_harris
real*8, parameter :: pi=3.1415926535897932384626433832795029d0
real*8:: a,w
real*8:: wsum
COMMON /glob/ wsum

    a = 2.0*pi/dble(n-1)

    w = 0.35875 - 0.48829*cos(a*dble(j)) + 0.14128*cos(2.0*a*dble(j)) - 0.01168*cos(3.0*a*dble(j))
    wsum = wsum+w
    win_blackman_harris=w
end function win_blackman_harris

function win_hanning(j, n)
integer,intent(in)::j,n
real*8:: win_hanning
real*8, parameter :: pi=3.1415926535897932384626433832795029d0
real*8:: a,w
real*8:: wsum
COMMON /glob/ wsum
    a = 2.0*pi/dble(n-1)
    
    w = 0.5 - 0.5*cos(a*j)
    wsum = wsum+w;
    win_hanning=w
end function win_hanning

function win_Gauss(j, n)
integer,intent(in)::j,n
real*8::win_Gauss
real*8::w
 w = (dble(j+1) - dble(n)/2.0) / (dble(n)/2.0) * 2.7183
 w=w*w 
 w=exp(-w)
 win_Gauss=w
end function win_Gauss

function win_FlatTop(j, n)
integer,intent(in)::j,n
real*8:: win_FlatTop
real*8, parameter :: pi=3.1415926535897932384626433832795029d0
real*8:: a,w
real*8:: wsum
COMMON /glob/ wsum

    a = 2.0*pi/dble(n-1)

    w = 1.0 &
     - 1.93293488969227 * cos(dble(j) * a) &
     + 1.28349769674027 * cos(dble(j) * a * 2.0 ) &
     - 0.38130801681619 * cos(dble(j) * a * 3.0 ) &
     + 0.02929730258511 * cos(dble(j) * a * 4.0 ) 
    wsum = wsum+w
    win_FlatTop=w
end function win_FlatTop

subroutine detrend(n,c)
integer,intent(in):: n
real*8,intent(inout)::c(n)
real*8:: a, b, tsqsum, ysum, t

a=0.d0
b=0.d0
tsqsum=0.d0
ysum=0.d0
t=0.d0

do i = 1, n
    ysum =ysum+c(i)
enddo

do i = 1, n
    t = i - n/2 - 0.5
    tsqsum =tsqsum+ t*t
    b = b+t*c(i)
enddo
b =b/tsqsum
a = ysum/n - b*(n-1)/2.0

do i = 1, n
c(i) = c(i)- (a + b*i)
enddo
if (b < -0.04 .and. b > 0.04) then
    write(6,*) " (warning) possibly significant trend in input series\n"
endif
end subroutine detrend
