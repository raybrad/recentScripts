!-------------------------------------------------------------------------
!> @brief  construct hamiltonian or overlap matrix in CSC format
!>          
!> @date 05-Aug-2015
!> @author ZHOU Yi
!-------------------------------------------------------------------------
!!!declare
#DEFINE __SUB_FORMKSMATRIX_MT__
subroutine formKSMatrix_mt(sys, shift, sk_data, hMat, sMat, crit)
    use parameters, only: AUEV, LDIM
    use dataType
    use cscMatrix_op

    implicit none
    !
    type(system), intent(in) :: sys !< system information
    real*8, optional, intent(in) :: shift(sys%nAt) !< Hamiltonian shift
    type(SKData), intent(in) :: sk_data !< SK file data
    type(zCSCMat), optional, intent(inout) :: hMat
    type(zCSCMat), intent(inout) :: sMat
    real*8, intent(in) :: crit
    real*8, external :: getDistance
!!!end declare
#include "./mt/mt.interface"
#include "./DFTB/DFTB.interface"
    !
    real*8, allocatable :: hamMatOneCol(:, :)
    real*8, allocatable :: hamMatOneBlk(:, :)
    real*8, allocatable :: overMatOneCol(:, :)
    real*8, allocatable :: overMatOneBlk(:, :)
    integer, allocatable :: cscMatRowInd(:)
    real*8 :: distance, skCut 
    integer :: matOneRowIndex, matOneColIndex
    integer :: nAtom = 1
    integer :: info
    integer :: iAtom, jAtom, iRow, jCol, i, j, inOrbs, jnOrbs
    integer :: pHamil, pOverl

    !
    skCut = sk_data%slkcut
    allocate(cscMatRowInd(sys%nOrbs), stat = info)
    if (info /= 0) then
        print*, ' Allocation fails in formKSMatrix for cscMatRowInd with code:', info
        stop
    endif

    allocate(overMatOneCol(sys%nOrbs, LDIM), overMatOneBlk(LDIM, LDIM), stat = info)
    if (info /= 0) then
        print*, ' Allocation fails in formKSMatrix for overMat with code:', info
        stop
    endif

    if (present(hMat)) then
        allocate(hamMatOneCol(sys%nOrbs, LDIM), hamMatOneBlk(LDIM, LDIM), stat = info)
        if (info /= 0) then
            print*, ' Allocation fails in formKSMatrix for hamMat with code:', info
            stop
        endif
    endif
    matOneColIndex = 1
    pHamil = 1
    pOverl = 1
    sMat%colPtr(1) = 1
    if (present(hMat)) hMat%colPtr(1) = 1
    do jAtom = 1, sys%nAt
        jnOrbs = sys%ind(jAtom + 1) - sys%ind(jAtom)
        matOneRowIndex = 1
        
        overMatOneCol = 0.d0
        if (present(hMat)) hamMatOneCol = 0.d0

        do iAtom = jAtom, sys%nAt
            ! Calculate distance between two atoms
            distance = getDistance(sys%xyz, iAtom, jAtom)

            if (distance > skCut) cycle
            ! Only fill-in elements of two atoms within skCut
            overMatOneBlk = 0.d0
            if (present(hMat)) hamMatOneBlk = 0.d0
            
            ! In order to form H, S, and K in CSC format, every time formHS for one pair of atoms, i.e. nAt1 = nAt2 = 1 in formHS
            if (present(hMat)) then
                call formHS(nAt1 = nAtom, nAt2 = nAtom, nOrbs1 = LDIM, nOrbs2 = LDIM, &
                        &   ind1 = sys%ind(iAtom), ind2 = sys%ind(jAtom), izp1 = sys%izp(iAtom), izp2 = sys%izp(jAtom), &
                        &   xyz1 = sys%xyz(1, iAtom), xyz2 = sys%xyz(1, jAtom), shift1 = shift(iAtom), shift2 = shift(jAtom), &
                        &   sk_data = sk_data, overMat = overMatOneBlk, hamMat = hamMatOneBlk)
            else ! sMat
                call formHS(nAt1 = nAtom, nAt2 = nAtom, nOrbs1 = LDIM, nOrbs2 = LDIM, &
                        &   ind1 = sys%ind(iAtom), ind2 = sys%ind(jAtom), izp1 = sys%izp(iAtom), izp2 = sys%izp(jAtom), &
                        &   xyz1 = sys%xyz(1, iAtom), xyz2 = sys%xyz(1, jAtom), sk_data = sk_data, overMat = overMatOneBlk)
            endif ! present(hMat)
            
            inOrbs = sys%ind(iAtom + 1) - sys%ind(iAtom)
            j = sys%Ind(iAtom) + 1
            do i = matOneRowIndex, matOneRowIndex + inOrbs - 1
                cscMatRowInd(i) = j 
                j = j + 1
            enddo !i
            overMatOneCol(matOneRowIndex : matOneRowIndex + inOrbs - 1, 1 : jnOrbs) = overMatOneBlk(1 : inOrbs, 1 : jnOrbs)
            if (present(hMat)) then
                hamMatOneCol(matOneRowIndex : matOneRowIndex + inOrbs - 1, 1 : jnOrbs) = hamMatOneBlk(1 : inOrbs, 1 : jnOrbs)
            endif

            matOneRowIndex = matOneRowIndex + inOrbs
        enddo !iAtom = jAtom, sys%nAt
        
        ! Fill-in matrix elements of column(jAtom) in CSC format
        j = 1
        do jCol = matOneColIndex, matOneColIndex + jnOrbs - 1
            do iRow = 1, matOneRowIndex - 1
                if (dabs(overMatOneCol(iRow, j)) >= crit) then
                    sMat%val(pOverl) = dcmplx(overMatOneCol(iRow, j))
                    sMat%rowInd(pOverl) = cscMatRowInd(iRow)
                    pOverl = pOverl + 1
                endif ! fill-in sMat

                if (present(hMat)) then 
                    if (dabs(hamMatOneCol(iRow, j)) >= crit) then
                        hMat%val(pHamil) = dcmplx(hamMatOneCol(iRow, j))
                        hMat%rowInd(pHamil) = cscMatRowInd(iRow)
                        pHamil = pHamil + 1
                    endif
                endif
            enddo !iRow
            sMat%colPtr(jCol + 1) = pOverl
            if (present(hMat)) hMat%colPtr(jCol + 1) = pHamil

            j = j + 1
        enddo !jCol

        matOneColIndex = matOneColIndex + jnOrbs
    enddo !jAtom = 1, sys%nAt

    !sMat%val = sMat%val * AUEV
    if (present(hMat)) hMat%val = hMat%val * AUEV
end subroutine formKSMatrix_mt

