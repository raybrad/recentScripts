!--------------------------------------------------------------------------------------
!> @brief:  Get the initial Fock_AO and denMat0 from ground state calculation
!>          
!> @date   22-JUL-2015
!> @author CSG
!--------------------------------------------------------------------------------------
!!!declare
!subroutine initDenMat0(nAt, nOrbs, Fock_AO, denMatR, overl, gMat)
subroutine initDenMat0(sys, tdMat,  denMatR)
  use dataType
  use heom_mod
  use mat_mod
  use parameters, only: AUEV

  implicit none
  
  type(system),     intent(inOut) :: sys
  type(tdMatrix),   intent(inOut) :: tdMat
  reaL*8,           intent(Out)   :: denMatR(sys%nOrbs, sys%nOrbs)
!!!end declare
  !------------------------------------------------------------------------------------------------------------------
  integer                 :: i, j, k, l, m, nAt, ndim, nOrbs, nOrb2
  integer                 :: istat
  real*8,  allocatable    :: tmpMat(:),ArrayMat(:,:)
 !==================================== read TAPE.matrix from ground state result  ==================================
  nAt   = sys%nAt 
  nOrbs = sys%nOrbs
  nOrb2 = nOrbs*(nOrbs+1)/2
  allocate(tmpMat(nOrb2), STAT=istat)
  if (istat /= 0) then
    write(6,*) 'allocate tmpMat failed in initDenMat0.f90'
    stop
  endif
  allocate(ArrayMat(nOrbs,nOrbs), STAT=istat)
  if (istat /= 0) then
    write(6,*) 'allocate tmpMat failed in initDenMat0.f90'
    stop
  endif

  open(unit=29, file='TAPE.matrix', status='old', form='unformatted', iostat=istat)
  if (istat /= 0) then
    write(6,*) 'TAPE.matrix open error in initMat when doing TD propagation in isoEpi'
    stop
  endif

  rewind 29
  read(29) ndim 
  if (ndim /= nOrbs) then
    write(6,*) '!!! nOrbs in TAPE.matrix not correct!!!!'
    stop
  endif

  read(29) (tmpMat(i), i=1,nOrb2)
  do i = 1, nOrbs
    do j = 1, i
     m = i*(i - 1)/2 + j
     ArrayMat(j,i) = tmpMat(m)*AUEV
     ArrayMat(i,j) = ArrayMat(j,i)
    enddo
  enddo
  call MatFromArray(tdMat%Fock_AO,ArrayMat,nOrbs,nOrbs)
  write(6,*) 'FockAO read from TAPE.matrix'

  read(29) (tmpMat(i), i=1,nOrb2)
  do i = 1, nOrbs
    do j = 1, i
     m = i*(i - 1)/2 + j
     ArrayMat(j,i) = tmpMat(m)*AUEV
     ArrayMat(i,j) = ArrayMat(j,i)
    enddo
  enddo
  call MatFromArray(tdMat%Fock_0,ArrayMat,nOrbs,nOrbs)
  write(6,*) 'Fock_0 read from TAPE.matrix'

  read(29) (tmpMat(i), i=1,nOrb2)
  do i = 1, nOrbs
    do j = 1, i
     m = i*(i - 1)/2 + j
     ArrayMat(j,i) = tmpMat(m)
     ArrayMat(i,j) = ArrayMat(j,i)
    enddo
  enddo
  call MatFromArray(tdMat%denMatR,ArrayMat,nOrbs,nOrbs)
  write(6,*) 'denMatR read from TAPE.matrix'
  !denMat%mb/nb/is/js are assigned in allocateRSDM(rsdm)

  read(29) (tmpMat(i), i=1,nOrb2)
  do i = 1, nOrbs
    do j = 1, i
     m = i*(i - 1)/2 + j
     ArrayMat(j,i) = tmpMat(m)
     ArrayMat(i,j) = ArrayMat(j,i)
    enddo
  enddo
  call MatFromArray(tdMat%overl,ArrayMat,nOrbs,nOrbs)
  write(6,*) 'Overlap Matrix read from TAPE.matrix'
!??
  read(29)  (sys%gMat(i), i=1, nAt*(nAt+1)/2)
  
  close(29)
  deallocate(tmpMat,ArrayMat)  



  end subroutine initDenMat0

