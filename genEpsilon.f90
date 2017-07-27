program genEpsilon
implicit none
integer::i,n
real*8::eps_b,omega,omegaStart,omegaEnd,domega
real*8::omega_p,gamma_0
real*8::depsr_0_1,omega_0_1,gamma_0_1
real*8::depsr_0_2,omega_0_2,gamma_0_2
real*8::lightspeed,n_M,eps_M,q,prefac,volume,surface
real*8::Q_ext,Q_sca,Q_abs,k_M,a,pi,hbar,eps_0,IMalpha
complex*16::eps,betaS,alpha
complex*16::eye=(0.d0,1.d0)

!Au
omegaStart=1.0d0
omegaEnd=6.0d0
domega=0.01d0
eps_b=3.559d0
omega_p=8.812d0
gamma_0=0.0752d0

depsr_0_1=2.912d0
omega_0_1=4.693d0
gamma_0_1=1.541d0

depsr_0_2=1.272d0
omega_0_2=3.112d0
gamma_0_2=0.525d0

!Ag
! omegaStart=1.24d0
! omegaEnd=5.0d0
! domega=0.01d0
!
! eps_b=3.189 
! omega_p=9.183 
! gamma_0=0.0179 
!
! depsr_0_1=0.4323 
! omega_0_1=4.668
! gamma_0_1=0.207
!
! depsr_0_2=0.2237
! omega_0_2=4.240 
! gamma_0_2=0.186

!///////////////////////////////////Drude+ Double Lorentz Model /////////////////////////////////////////////////////
!/                                                                                                                   /
!/                    (omega_p)^2               depsr_0_1*(omega_0_1)^2               depsr_0_2*(omega_0_2)^2        /  
!/  eps(w) = eps_b - ------------------   + --------------------------------- + ---------------------------------    /
!/                   w^2 + j w gamma_0      omega_0_1^2 - w^2 - j w  2 gamma_0_1   omega_0_2^2-w^2-j w 2 gamma_0_2   /
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
a=10.d0 !nm
pi=3.141592453d0
lightspeed=2.99792458d8
q=1.602176487d-19
eps_M=1.d0
hbar=0.6582d0
eps_0=8.854d-12
n_M=sqrt(eps_M)  ! 
n=int((omegaEnd-omegaStart)/domega)+1
volume=4.d0/3.d0*pi*a**3.d0
surface=pi*a**2.d0
open(16,file='epsilon.dat')
open(17,file='ext-sca.dat')
do i=1,n
    eps=(0.d0,0.d0) 
    omega=omegaStart+(i-1)*domega !eV
    eps=eps_b-omega_p**2.d0/(omega**2.d0+eye*omega*gamma_0)+ &
              depsr_0_1*omega_0_1**2.d0/(omega_0_1**2.d0-omega**2.d0-eye*2.d0*omega*gamma_0_1)+ &
              depsr_0_2*omega_0_2**2.d0/(omega_0_2**2.d0-omega**2.d0-eye*2.d0*omega*gamma_0_2)
    write(16,'(3(e20.12,x))') omega,dble(eps),aimag(eps)
   
    !\alpha=\epsilon_0*3*V*\frac{\espilon_r- \epsilon_{env}{\espilon_r +2*\epsilon_{env}}
    ! C_sca=k^4/(6\pi)*|\alpha/\epsilon_0|^2, Q_sca=C_sca/S
    ! C_abs=k*Im(\alpha/\epsilon_0), Q_abs=C_abs/S
    ! k=2\pi*n_M/\lambda ,S=\pi r^2
    alpha= eps_0*3.d0*volume*(eps-1.d0*eps_M)/(eps+2.d0*eps_M)  !no dimension
    ! 2*pi n_M  /lambda = 2*pi*n_M/(1239.84187/omega)
    k_M=2*pi*n_M*omega/(1239.84187d0) !1/nm 
    Q_sca=k_M**4.d0/(6.d0*pi)*abs(alpha/eps_0)**2.d0/surface    !   
    Q_abs=k_M*aimag(alpha/eps_0)/surface
    Q_ext=Q_abs+Q_sca
    IMalpha=aimag(alpha/eps_0)/surface
    write(17,'(6(e20.12,x))') omega,Q_abs,Q_sca,Q_ext,Q_abs/Q_ext,IMalpha
enddo
close(16)
close(17)
end
