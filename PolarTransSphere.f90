program PolarTransSphere
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
complex*16 ::Umatrix(9,9),UinvMat(9,9),alpha_c,gamma_c
complex*16 ::t_polSphere(9),tmpMat(9),PisoMat(3,3),PanisoMat(3,3),PisoMat2(3,3),PanisoMat2(3,3)
complex*16 ::Piso(9),Paniso(9)
complex*16,allocatable ::t_polw(:,:)
complex*16,allocatable:: alpha2(:),gamma2(:)
real*8,allocatable:: tt(:)
real*8::nmomega,omega,omega_start,omega_end,domega,dep 
real*8:: dtt,tt_min,tt_max
real*8:: time,t_polMean(9)
real*8,allocatable::t_pol(:,:)
real*8:: alpha,alpha_p,gamma_p,depolRatio,tempfactor,qmfactor
real*8:: wvn,lambda_0,lambda_z,lambda_t,temperature,tempfac
real*8:: gscr,ascr,cperu
integer*4:: i,j,k,line_output
integer*4:: line_input
real*8:: win_blackman_harris,norm
real*8:: wsum
real*8:: win_Gauss,win_hanning,win_FlatTop
COMMON /glob/ wsum

namelist /freq/omega_start,omega_end,line_input,line_output,dep,wvn,temperature


wsum=0.d0
wvn=514.d0
temperature=298.d0
rewind 5
read(5,freq,end=102)
102 continue 
PRINT *,'starting freq. in eV',omega_start
PRINT *,'end freq. in eV',omega_end
PRINT *,'number of lines in the file',line_input
PRINT *,'number of output points',line_output

write(6,*) 'incident wave(nm)',wvn
write(6,*) 'temperature: ',temperature
write(6,*) 'conver',conver
write(6,*) 'vbyt',vbyt
!-------------------------------------------------------------------
! set default wave number (Ar ion laser, wavelength =  514.5 nm 
! => wavenumber lambda_0 = 19436.3459669582118561710d0 1/cm)
!------------------------------------------------------------------- thetarad=theta*2.d0*pi/360.d0
lambda_0=1.d0/wvn*1.0d7
write(6,*) 'incident wave(cm^-1)',lambda_0
!
!
!
Umatrix=     dcmplx( 0.d0,0.d0)
Umatrix(1,1)=dcmplx(-1.d0/dsqrt(3.d0),0.d0)
Umatrix(1,2)=dcmplx(-1.d0/dsqrt(3.d0),0.d0)
Umatrix(1,3)=dcmplx(-1.d0/dsqrt(3.d0),0.d0)
Umatrix(2,4)=dcmplx(-1.0/2.0,0.d0)
Umatrix(2,5)=dcmplx( 1.0/2.0,0.d0)
Umatrix(2,6)=dcmplx( 0.d0,-1.0/2.0)
Umatrix(2,7)=dcmplx( 0.d0, 1.0/2.0)
Umatrix(3,8)=dcmplx( 0.d0, 1.0/dsqrt(2.d0))
Umatrix(3,9)=dcmplx( 0.d0,-1.0/dsqrt(2.d0))
Umatrix(4,4)=dcmplx(-1.0/2.0,0.d0)
Umatrix(4,5)=dcmplx( 1.0/2.0,0.d0)
Umatrix(4,6)=dcmplx( 0.d0, 1.0/2.0)
Umatrix(4,7)=dcmplx( 0.d0,-1.0/2.0)
Umatrix(5,1)=dcmplx( 1.0/2.0,0.d0)
Umatrix(5,2)=dcmplx(-1.0/2.0,0.d0)
Umatrix(5,8)=dcmplx( 0.d0,1.0/2.0)
Umatrix(5,9)=dcmplx( 0.d0,1.0/2.0)
Umatrix(6,4)=dcmplx(-1.0/2.0,0.d0)
Umatrix(6,5)=dcmplx(-1.0/2.0,0.d0)
Umatrix(6,6)=dcmplx( 0.d0,-1.0/2.0)
Umatrix(6,7)=dcmplx( 0.d0,-1.0/2.0)
Umatrix(7,1)=dcmplx(-1.d0/dsqrt(6.d0),0.d0)
Umatrix(7,2)=dcmplx(-1.d0/dsqrt(6.d0),0.d0)
Umatrix(7,3)=dcmplx( 2.d0/dsqrt(6.d0),0.d0)
Umatrix(8,4)=dcmplx( 1.0/2.0,0.d0)
Umatrix(8,5)=dcmplx( 1.0/2.0,0.d0)
Umatrix(8,6)=dcmplx( 0.d0,-1.0/2.0)
Umatrix(8,7)=dcmplx( 0.d0,-1.0/2.0)
Umatrix(9,1)=dcmplx( 1.0/2.0,0.d0)
Umatrix(9,2)=dcmplx(-1.0/2.0,0.d0)
Umatrix(9,8)=dcmplx( 0.d0,-1.0/2.0)
Umatrix(9,9)=dcmplx( 0.d0,-1.0/2.0)
UinvMat= conjg(transpose(Umatrix))
allocate(tt(line_input),t_pol(line_input,9))
allocate(t_polw(line_output,9))
allocate(alpha2(line_output),gamma2(line_output))
!
OPEN(60,FILE='tdPolar.sphere',STATUS='OLD')
t_polMean=0.d0
do i=1,line_input
 read(60,*)time,(t_pol(i,j),j=1,9)
 tt(i)=time
 t_polMean(1:9)=t_polMean(1:9)+t_pol(i,1:9)
enddo
t_polMean = t_polMean/dble(line_input)
!zero
do i=1,line_input
    t_pol(i,1:9) = t_pol(i,1:9) - t_polMean(1:9)
enddo
!window
do i=1,line_input
   t_pol(i,1:9)=t_pol(i,1:9)*win_blackman_harris(i,line_input)
enddo

domega=(omega_end-omega_start)/(line_output-1.d0)
omega=omega_start-domega
!
t_polw=(0.d0,0.d0)
do i=1,line_output
  omega=omega+domega
  do j=1,line_input-1
     dtt=tt(j+1)-tt(j)
     do k=1,9
        t_polw(i,k)=t_polw(i,k)+exp(dcmplx(0.d0,omega*tt(j)/hbar))*t_pol(j,k)*dtt
     enddo
  enddo
enddo

do i=1,line_output
   t_polSphere=(0.d0,0.d0)
   do j=1,9
       do k=1,9
           t_polSphere(j)=t_polSphere(j)+Umatrix(j,k)*t_polw(i,k)
       enddo
   enddo
   tmpMat=dcmplx(0.d0,0.d0)
   tmpMat(1)=t_polSphere(1)
   Piso=(0.d0,0.d0)
   do j=1,9
       do k=1,9
           Piso(j)=Piso(j)+UinvMat(j,k)*tmpMat(k)
       enddo
   enddo
     PisoMat(1,1)=Piso(1)   
     PisoMat(2,2)=Piso(2)   
     PisoMat(3,3)=Piso(3)   
     PisoMat(1,3)=Piso(4)   
     PisoMat(3,1)=Piso(5)   
     PisoMat(2,3)=Piso(6)   
     PisoMat(3,2)=Piso(7)   
     PisoMat(1,2)=Piso(8)   
     PisoMat(2,1)=Piso(9)   

     PisoMat2=dcmplx(0.d0,0.d0)
     do j=1,3
        do k=1,3
           PisoMat2(j,k) = PisoMat2(j,k)+PisoMat(j,k)*PisoMat(k,j)
        enddo
     enddo
     alpha2(i)=dcmplx(0.d0,0.d0)
     do j=1,3
        alpha2(i)=alpha2(i)+PisoMat2(j,j)
     enddo
     alpha2(i)=alpha2(i)/3.d0
    !alpha2(i) = Piso(1)*Piso(1)
     tmpMat=dcmplx(0.d0,0.d0)
     tmpMat(5:9)=t_polSphere(5:9)
     Paniso=(0.d0,0.d0)
     do j=1,9
         do k=1,9
             Paniso(j)=Paniso(j)+UinvMat(j,k)*tmpMat(k)
         enddo
     enddo
     PanisoMat(1,1)=Paniso(1)   
     PanisoMat(2,2)=Paniso(2)   
     PanisoMat(3,3)=Paniso(3)   
     PanisoMat(1,3)=Paniso(4)   
     PanisoMat(3,1)=Paniso(5)   
     PanisoMat(2,3)=Paniso(6)   
     PanisoMat(3,2)=Paniso(7)   
     PanisoMat(1,2)=Paniso(8)   
     PanisoMat(2,1)=Paniso(9)   

     PanisoMat2=dcmplx(0.d0,0.d0)
     do j=1,3
        do k=1,3
           PanisoMat2(j,k) = PanisoMat2(j,k)+PanisoMat(j,k)*PanisoMat(k,j)
        enddo
     enddo
     gamma2(i)=dcmplx(0.d0,0.d0)
     do j=1,3
        gamma2(i)=gamma2(i)+PanisoMat2(j,j)
     enddo
     gamma2(i)=3.d0/2.0*gamma2(i)
enddo
close(60)

!
OPEN(72,FILE='freqPolar.sphere',STATUS='Replace')
tt_min=tt(1)
tt_max=tt(line_input)
dtt=tt(2)-tt(1)
domega=(omega_end-omega_start)/(line_output-1.d0)
omega=omega_start-domega
!
do i=1,line_output
  omega=omega+domega
  
  alpha_p=abs(alpha2(i))
  gamma_p=abs(gamma2(i))
  
  nmomega=1239.84187/omega
  lambda_z=omega*1.d7/1239.84187
  ! if(lambda_z < 400.d0) then
  !       s=(0.d0,0.d0)
  ! endif
  
  if (temperature.le.1.0d-10) then
     tempfac=1.0d0
  else
     tempfac=1.0d0/(1.0d0-exp(-vbyt*lambda_z/temperature))
  end if

  lambda_t=tempfac*(lambda_0-lambda_z)**4/lambda_z *1.0d6
  tempfactor = ev2j/(boltz*temperature)*omega
  ! qmfactor   = tempfactor/(1.d0-dexp(-tempfactor))*omega**2
   qmfactor   = 1.d0*omega**2
  !qmfactor   = 1.d0
  gscr = gamma_p/15.0d0*conver*lambda_t*qmfactor
  ascr = (alpha_p+4.0d0/45.0d0*gamma_p)*conver*lambda_t*qmfactor
  depolRatio=gscr/ascr
  cperu=(gscr+ascr)*1d4*1d30 !m2/sr-> 10^-30 cm^2/sr
  WRITE(72,'(7e24.14)')lambda_z,alpha_p,gamma_p,gscr,ascr,cperu,depolRatio
enddo
close(72)
    
       

       


end program PolarTransSphere

function win_blackman_harris(j, n)
integer,intent(in)::j,n
real*8:: win_blackman_harris
real*8, parameter :: pi=3.1415926535897932384626433832795029d0
real*8:: a,w
real*8:: wsum
COMMON /glob/ wsum
    a = 2.0*pi/dble(n-1)

    w = 0.35875 - 0.48829*cos(a*j) + 0.14128*cos(2*a*j) - 0.01168*cos(3*a*j)
    wsum = wsum+w;
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
