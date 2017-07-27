!-------------------------------------------------------------------------
!> @brief  calculate transmission coefficient and integrate to obtain current
!>          
!> @date  01-Jun-2013
!> @author CY
!-------------------------------------------------------------------------
!!!declare
#DEFINE __SUB_INTEGRATETC__
subroutine integrateTC(lPDOS, lTCoef, nGrid, tElec, eIntBot, eIntTop, eFermi, eImag, sys, shift, iLead, jLead, voltage, lCondense, lead, sk_data, selfMat, current)
use parameters, only: FSMALL, PI, HBAR, EV2J2, MGAULEG1, dpt, dw
use dataType
use blkMatrix
#ifdef __MPI__
use parallel, only: iRank, nNode
#endif
implicit none
!
! Integrate current to be
! 
! J = (-1) * e * 8 / h * Int_{eIntBot}^{eIntTop} T_{AB}'(E) dE
!
! T_{AB}'(E) = Tr( G^r * DIMAG(SELF^r_B) * G^a * DIMAG(SELF^r_A) )
!
! along (eIntBot, eImag) -> (eIntTop, eImag)
!
logical, intent(in)                   :: lPDOS              !< \b input: output PDOS
logical, intent(in)                   :: lTCoef             !< \b input: output Transmission Coefficients
integer, intent(in)                   :: nGrid              !< \b input: number of intervals for real axis integration
real*8,  intent(in)                   :: tElec              !< \b input: temperature
real*8,  intent(in)                   :: eIntBot            !< \b input: lower bound of energy integration
real*8,  intent(in)                   :: eIntTop            !< \b input: upper bound of energy integration
real*8,  intent(in)                   :: eFermi             !< \b input: Fermi energy
real*8,  intent(in)                   :: eImag              !< \b input: imaginary energy used in energy integration
type(system), intent(in)              :: sys                !< \b input: system information
real*8,  intent(in)                   :: shift(sys%nAt)     !< \b input: Hamiltonian shift
integer, intent(in)                   :: iLead              !< \b input: contact A
integer, intent(in)                   :: jLead              !< \b input: contact B
real*8,  intent(in)                   :: voltage(sys%nLead) !< \b input: bias voltage (actually potential energy)
logical, intent(in)                   :: lCondense          !< \b input: condense contact Hamiltonian
type(system), intent(in)              :: lead(sys%nLead)    !< \b input: contacts information
type(SKData), intent(in)              :: sk_data            !< \b input: SK file data
type(selfEnergyMatrix), intent(inout) :: selfMat(sys%nLead) !< \b in/out: self energy related matrices
real*8,  intent(out)                  :: current            !< \b output: current
!!!end declare
#include "./DFTB/DFTB.interface"
#include "./Steady/Steady.interface"
!
complex*16, parameter   :: cOne = dcmplx(1d0, 0d0)
complex*16, parameter   :: cZero= dcmplx(0d0, 0d0)
!
character*8             :: pFile
character*13            :: tFile
logical                 :: lTemp
integer                 :: i, j, k, istat
integer                 :: nsplit, npoint
integer                 :: count1, count2, count3, count4, countr
real*8                  :: fermiDirac, angStart, angEnd, pDOS(0:sys%nType), fermi
real                    :: etime, elapsed(2), total
real                    :: cpu1, cpu2, cpu3, cpu4
complex*16              :: const, tCoef, ePt, ePtS(sys%nLead)
real*8, allocatable     :: angPt(:), weight(:)
real*8, allocatable     :: qOrb(:)
type(dBlkMat)           :: sMat
type(cBlkMat)           :: gMat
type(dBlkSequence)      :: gamma(sys%nLead)
type(cBlkSequence)      :: gOuter(sys%nLead)
type(cBlkSequence)      :: gCorner
!
current = 0d0
if ((eIntTop - eIntBot) < FSMALL) return
!
total = etime(elapsed)
write(6,*)
write(6,*) ' Enter integrateTC  '
write(6,*) ' start: total=', total, ' user=', elapsed(1),' system=', elapsed(2)
call cpu_time(cpu1)
call system_clock(count=count1, count_rate=countr)
!
if (tElec >= 1d-1) then
  lTemp = .true.
else
  lTemp = .false.
end if
!
#ifdef __MPI__
if (iRank < 10) then
  write(pFile, '(A7,I1)') 'PDOS.00', iRank
  write(tFile, '(A6,I1,A1,I1,A3,I1)') 'TCoef_', iLead, '-', jLead, '.00', iRank
else if (iRank < 100) then
  write(pFile, '(A6,I2)') 'PDOS.0',  iRank
  write(tFile, '(A6,I1,A1,I1,A2,I2)') 'TCoef_', iLead, '-', jLead, '.0',  iRank
else if (iRank < 1000) then
  write(pFile, '(A5,I3)') 'PDOS.',   iRank
  write(tFile, '(A6,I1,A1,I1,A1,I3)') 'TCoef_', iLead, '-', jLead, '.',   iRank
end if
#else
write(pFile, '(A8)') 'PDOS.dat'
write(tFile, '(A6,I1,A1,I1,A4)') 'TCoef_', iLead, '-', jLead, '.dat'
#endif
if ( lTCoef ) then
  open(28, file=tFile, status='unknown')
end if
if ( lPDOS ) then
  open(29, file=pFile, status='unknown')
  open(30, file='PDOS.orb', status='unknown')
  do i = 1, sys%nLead
    allocate(gOuter(i)%seq(sys%nSigma(i),lead(i)%nOrbs), STAT=istat)
    if (istat /= 0) then
      write(6,*)"allocation gOuter error in integrateTC.f90", i
      stop
    else
      gOuter(i)%seq(:,:) = cZero
    end if
  enddo
  allocate(qOrb(sys%nOrbs), STAT=istat)
  if (istat /= 0) then
    write(6,*)"allocation 1 error in integrateTC.f90  "
    stop
  end if
end if
!
nsplit = nGrid
npoint = MGAULEG1
!
allocate(gCorner%seq(sys%nSigma(iLead),sys%nSigma(jLead)), STAT=istat)
if (istat /= 0) then
  write(6,*)"allocation 2 error in integrateTC.f90  "
  stop
end if
allocate(angPt(nsplit*npoint), weight(nsplit*npoint), STAT=istat)
if (istat /= 0) then
  write(6,*)"allocation 3 error in integrateTC.f90  "
  stop
end if
call allocateBlkMat(sys%nBlk, sys%blkSize, gMat)
do i = 1, sys%nLead
  allocate(gamma(i)%seq(sys%nSigma(i),sys%nSigma(i)), STAT=istat)
  if (istat /= 0) then
    write(6,*)"allocation 4 error in integrateTC.f90  "
    stop
  end if
enddo
!
k = 0
do i = 1, nsplit
  angStart = (eIntTop - eIntBot) / dble(nsplit) * dble(i - 1) + eIntBot
  angEnd   = (eIntTop - eIntBot) / dble(nsplit) * dble(i)     + eIntBot
!
  do j = 1, npoint
    k = k + 1
    angPt(k)  = ((angEnd - angStart) * dpt(j) + (angEnd + angStart)) * 5.d-1
    weight(k) = dw(j)
  enddo
enddo
!
! Energy integration starts here
!
write(6,1201)
#ifdef __MPI__
do i = iRank+1, nsplit*npoint, nNode
#else
do i = 1, nsplit*npoint
#endif
  call cpu_time(cpu3)
  call system_clock(count=count3, count_rate=countr)
!
  ePt  = dcmplx(angPt(i), eImag)
  const = weight(i) * (eIntTop - eIntBot) * 5.d-1 / dble(nsplit)
!
! Fermi distribution
!
!-> modify for extra contacts
  ePtS(1:sys%nLead) = dble(ePt) - eFermi - voltage(1:sys%nLead)
  fermi = 1d0
  if (lTemp) then
    fermi = dabs( fermiDirac(ePtS(2), tElec) - fermiDirac(ePtS(1), tElec) )
    const = const * fermi
  end if
!-> modify for extra contacts
!
  ePtS(1:sys%nLead) = ePt ! non-WBL	!overide the above one
  call formKSMatrix(ePt=ePt, sys=sys, shift=shift, sk_data=sk_data, kMat=gMat)
!
! write PDOS(E) for device region 
!
  if ( lPDOS ) then
!
! Note: greenEqBlk here output gMat%Bi as Si instead of Di*Si
!       construct Di*Si before mullikenBlk
!
    if ( sys%nBlk == 1 ) then  ! full matrix
      call greenEq(sys=sys, lead=lead, ePt=ePt, ePtS=ePtS, const=cOne, lCondense=lCondense, selfMat=selfMat, gOuter=gOuter, gMat=gMat, gamma=gamma, iLead=iLead, jLead=jLead, gCorner=gCorner) ! gOuter required for PDOS
    else
      call greenEqBlk(sys=sys, lead=lead, ePt=ePt, ePtS=ePtS, const=cOne, lCondense=lCondense, selfMat=selfMat, gOuter=gOuter, gMat=gMat, gamma=gamma, iLead=iLead, jLead=jLead, gCorner=gCorner) ! gOuter required for PDOS
    end if
    call allocateBlkMat(sys%nBlk, sys%blkSize, sMat)
    call formKSMatrix(sys=sys, sk_data=sk_data, sMat=sMat)
    qOrb(1:sys%nOrbs) = 0d0
    call mullikenBlk(sys%nBlk, sys%blkSize, gMat, sMat, qOrb, 2) ! take imaginary part of G^r
    call deallocateBlkMat(sys%nBlk, sys%blkSize, sMat)
    do k = 1, sys%nLead
      gOuter(k)%seq(:,:) = dcmplx(dimag(gOuter(k)%seq(:,:)),0d0)
      call mullikenOuterBlk(gOuter(k)%seq, sys%nSigma(k), selfMat(k)%sLD, lead(k)%nOrbs, qOrb(lead(k)%nStart))
    enddo
!
    pDOS = 0d0
    do j = 1, sys%nAt
      do k = sys%ind(j)+1, sys%ind(j+1)
        pDOS(0) = pDOS(0) + qOrb(k)
        pDOS(sys%izp(j)) = pDOS(sys%izp(j)) + qOrb(k)
      enddo
    enddo
    pDOS = pDOS * -2d0 / PI
    write(29,'(9ES20.12)')      dble(ePt) - eFermi, (pDOS(j),             j=0,sys%nType)
    write(30,'(999999ES20.12)') dble(ePt) - eFermi, (qOrb(j) * -2d0 / PI, j=1,sys%nOrbs)
  else
    if ( sys%nBlk == 1 ) then  ! full matrix
      call greenEq(sys=sys, lead=lead, ePt=ePt, ePtS=ePtS, const=cOne, lCondense=lCondense, selfMat=selfMat, gMat=gMat, gamma=gamma, iLead=iLead, jLead=jLead, gCorner=gCorner) ! gOuter only needed for PDOS
    else
      call greenEqBlk(sys=sys, lead=lead, ePt=ePt, ePtS=ePtS, const=cOne, lCondense=lCondense, selfMat=selfMat, gMat=gMat, gamma=gamma, iLead=iLead, jLead=jLead, gCorner=gCorner) ! gOuter only needed for PDOS
    end if
  end if
!
! Calculate transmission coefficient spectrum T(E)
!
! T(E) = Tr[ G^r * 2Im(\Sigma_B) * G^a * 2Im(\Sigma_A) ]
!
  call tCoefficient(sys%nSigma(iLead), sys%nSigma(jLead), gamma(iLead)%seq, gamma(jLead)%seq, gCorner%seq, tCoef)
  if ( lTCoef ) then
    write(28,'(3(e15.7,3x))') dble(ePt) - eFermi, dble(tCoef * 4.d0), dble(tCoef * 4.d0) * fermi
  end if
!
  current = current + dble(tCoef * const)
!
  call cpu_time(cpu4)
  call system_clock(count=count4, count_rate=countr)
  write(6,'(9X,I6,3X,4(2X,F10.4),A)') i, ePt, (cpu4 - cpu3)/60.0, dble((count4-count3)/countr)/60d0, ' min'
enddo
!
if ( lTCoef ) then
  close(28)
end if
if ( lPDOS ) then
  close(29)
  close(30)
  do i = 1, sys%nLead
    deallocate(gOuter(i)%seq, STAT=istat)
  enddo
  deallocate(qOrb, STAT=istat)
end if
!
! Analyze current
!
current = current * 4.d0 / HBAR / PI * EV2J2 * 1.d24
!
deallocate(angPt, weight, gCorner%seq, STAT=istat)
call deallocateBlkMat(sys%nBlk, sys%blkSize, gMat)
do i = 1, sys%nLead
  deallocate(gamma(i)%seq, STAT=istat)
enddo
!
npoint = nGrid * MGAULEG1
call cpu_time(cpu2)
call system_clock(count=count2, count_rate=countr)
write(6,'(76("-"))')
write(6,'(/,2X,A,2(2X,F10.4),A,I8,A)') '  Total time : ', (cpu2 - cpu1)/60.0, dble((count2-count1)/countr)/60d0, ' min for ', npoint, ' energy points'
write(6,  '(2X,A,2(2X,F10.4),A,I8,A)') ' Averge time : ', (cpu2 - cpu1)/npoint/60.0, dble((count2-count1)/countr)/npoint/60d0, ' min for each energy point'
!
total = etime(elapsed)
write(6,*)
write(6,*) ' Leave integrateTC  '
write(6,*) ' start: total=', total, ' user=', elapsed(1),' system=', elapsed(2)
!
1201 format(/,10X,'Energy',10X,'Real',7X,'Imag.',5X,'CPU Time',5X,'System Time')
end subroutine integrateTC
