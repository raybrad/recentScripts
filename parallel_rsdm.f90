module parallel_rsdm
  use heom_mod
  use parallel
  implicit none
  private
  !------------------------
  ! public functions
  !------------------------
  public startMPI_RSDM
  public inCol
  public inRow
  public Tier2_triangular   !< only calculate 2nd Tier for ia2 <= ia 
  public Tier1_AllGather    !< allGather 1st Tier, used in columnwise scheme
  public Tier1_Bcast
  public Mat_AllReduce
  public Mat_Reduce
  public stopMPI_RSDM
  public mpi_comm_col, mpi_comm_row, mpi_comm_diag
  !
  integer :: iCol, iRow, nCol, nRow
  !
  integer,allocatable,dimension(:) :: cntCol  !< count   
  integer,allocatable,dimension(:) :: cntRow  !< count 
  integer,allocatable,dimension(:) :: dispCol !< displacement 
  integer,allocatable,dimension(:) :: dispRow !< displacement
  !
  character*30 :: mpi_scheme
  !
  integer :: mpi_type_mat
  !
  integer :: mpi_comm_col
  integer :: mpi_comm_row
  integer :: mpi_comm_diag


  contains
  subroutine startMPI_RSDM(rsdm, scheme)
    type(rsdm_T)    ,intent(in) :: rsdm
    character(len=*),intent(in) :: scheme
    !
    integer :: nLead, nk, nTotal, blkSize, i, iq, nBlk, iDiag
    nLead  = rsdm%nLead
    nk     = rsdm%nk
    nTotal = nLead*nk
    if(trim(scheme)=="column")then
      iCol = iRank
      nCol = nNode
      iRow = 0
      nRow = 1
    else !if(trim(scheme)=="checkerboard")then

      ! BCM scheme ref: http://naturalunits.blogspot.hk/2013/07/on-parallelism-of-symmetric-sparse.html
      nCol = nint( ( dsqrt(8d0*nNode+1)-1 )/2 )
      nRow = nCol
      if( nCol*(nCol+1)/2 /= nNode)then
        write(6,*)'nNode = ',nNode
        write(6,*)'nNode has to be n(n+1)/2 to use checkerboard scheme (e.g. 6, 10, 15, etc)'
        stop
      endif
      call BCMBoard(nCol,nRow,iRank,iCol,iRow)
      ! note: iCol, iRow start from 0
      write(6,*)'Col:',iCol+1,'/',nCol
      write(6,*)'Row:',iRow+1,'/',nRow
      call PrintBoard(6,nNode,nCol,nRow)
    endif
    !
    allocate(cntCol(nCol), dispCol(nCol))
    allocate(cntRow(nRow), dispRow(nRow))
    cntCol(1:nCol) = nTotal/nCol
    cntCol(1:mod(nTotal,nCol)) =  cntCol(1:mod(nTotal,nCol)) + 1
    cntRow(1:nRow) = nTotal/nRow
    cntRow(1:mod(nTotal,nRow)) =  cntRow(1:mod(nTotal,nRow)) + 1
    dispCol = 0
    dispRow = 0
    do i=2,nCol
      dispCol(i) = dispCol(i-1) + cntcol(i-1)
    enddo
    do i=2,nRow
      dispRow(i) = dispRow(i-1) + cntRow(i-1)
    enddo
    !write(6,*)'cntCol : ',cntCol
    !write(6,*)'dispCol: ',dispCol
    !
    ! Create column, row, diag communicator
    iDiag = (iRow - iCol)
    call mpi_comm_split(mpi_comm_world, iRow,  abs(iDiag), mpi_comm_row,  ierr)
    call mpi_comm_split(mpi_comm_world, iCol,  abs(iDiag), mpi_comm_col,  ierr)
    call mpi_comm_split(mpi_comm_world, iDiag, abs(iDiag), mpi_comm_diag, ierr)
    !
  end subroutine
  subroutine stopMPI_RSDM
    deallocate(cntCol, cntRow, dispCol, dispRow)
  end subroutine

  function inCol(ia)
    integer :: ia
    logical :: inCol
    !inCol = ( mod(ia,nCol) == iCol )
    if( dispCol(iCol+1) < ia .and. ia <= dispCol(iCol+1)+cntCol(iCol+1))then
      inCol = .true.
    else
      inCol = .false.
    endif
  end function
  function inRow(ia)
    integer :: ia
    logical :: inRow
    !inRow = ( mod(ia,nRow) == iRow )
    if( dispRow(iRow+1) < ia .and. ia <= dispRow(iRow+1)+cntRow(iRow+1))then
      inRow = .true.
    else
      inRow = .false.
    endif
  end function
  subroutine BCMBoard(nCol,nRow,iRank,iCol,iRow)
    integer,intent(in) ::  nCol,nRow,iRank
    integer,intent(out) :: iCol, iRow
    !
    integer :: nBlk,nBlk2, iq  
    nBlk  = ceiling( (nCol+1)*0.5d0 )
    nBlk2 = (nCol+1)/2 
    if( iRank <= nBlk*(nCol+1-nBlk) )then
      iCol = iRank/nBlk + 1
      iRow = mod(iRank,nBlk) + iCol
    else
      iq = iRank - nBlk*(nCol + 1 - nBlk) 
      iCol = nCol + 2 - nBlk + iq/ nBlk2
      iRow = mod(mod(iq,nBlk2)+iCol-1,nCol) + 1
    endif
    iCol = iCol - 1
    iRow = iRow - 1
  end subroutine
  subroutine PrintBoard(hFile, nNode, nCol, nRow)
    integer,intent(in) ::  hFile, nNode, nCol, nRow
    !
    integer :: i,j,ic,ir,jCol(nNode), jRow(nNode)
    integer,parameter :: nDigit = 4
    integer  :: rank(nCol,nRow)
    character(len=nDigit) :: id 
    character,parameter :: char_corner = '+'
    character,parameter :: char_hori = '-'
    character,parameter :: char_vert = '|'
    rank = -1
    do i=1,nNode
      call BCMBoard(nCol,nRow, i-1, jCol(i), jRow(i))
      jCol(i) = jCol(i) + 1
      jRow(i) = jRow(i) + 1
      rank(jCol(i),jRow(i))=i 
    enddo
    do ir = 1, nRow
      do ic = 1, nCol
        write(hFile,'(A\)'),char_corner
        do j=1,nDigit+2
          write(hFile,'(A\)'),char_hori
        enddo
      enddo
      write(hFile,'(A)'),char_corner
      !
      do ic = 1, nCol
        write(hFile,'(A,A\)'),char_vert,' '
        if(rank(ic,ir) >= 0)then
          write(hFile,'(I4,A\)'),rank(ic,ir),' '
        else
          write(hFile,'(4A,A\)'),(' ',i=1,nDigit),' '
        endif
      enddo
      write(hFile,'(A)'),char_vert
    enddo
    !
    do ic = 1, nCol
      write(hFile,'(A\)'),char_corner
      do j=1,nDigit+2
        write(hFile,'(A\)'),char_hori
      enddo
    enddo
    write(hFile,'(A)'),char_corner
    !write(6,*)'+------+------+------+------+-------'
    !write(6,*)'|  20  |  01  | 1234 |--------------'
    !write(6,*)'+------+------+------.--------------'
  end subroutine
  subroutine printLine()
  end subroutine

  function MatToMPI(mat)
    type(mat_T),intent(in) :: mat
    integer :: MatToMPI
    !
    integer,parameter :: nVar = 4
    integer :: typelist(nVar)
    integer :: block_lengths(nVar)
    integer(kind=mpi_address_kind) :: address(0:nVar), displacements(nVar)
    integer :: i

    typelist(1:4) = mpi_integer
    !
    block_lengths(1:4) = 1
    call mpi_get_address(mat   ,address(0),ierr)
    call mpi_get_address(mat%is,address(1),ierr)
    call mpi_get_address(mat%js,address(2),ierr)
    call mpi_get_address(mat%mb,address(3),ierr)
    call mpi_get_address(mat%nb,address(4),ierr)
    !
    do i=1,nVar
      displacements(i) = address(i)-address(0)
    enddo
    !write(6,*)'displacements',displacements
    !
    call mpi_type_struct(nVar, block_lengths, displacements, typelist, MatToMPI, ierr)

  end function

  function Tier2_triangular
    logical :: Tier2_triangular
    if(trim(mpi_scheme) == "column")then
      Tier2_triangular = .false.
    else
      Tier2_triangular = .true.
    endif 
  end function
  subroutine Tier1_AllGather(rsdm)
    type(RSDM_T) :: rsdm
    !
    integer,parameter :: nVar = 4
    integer :: intBuff(nVar,rsdm%nk * rsdm%nLead)
    integer :: i,n,ic,mb,nb,is,ie,nTotal
    complex*16,pointer :: matBuff(:), matPointer(:,:) 
    integer :: ind(rsdm%nk * rsdm%nLead + 1)
    integer :: cntMat(nCol)
    integer :: dispMat(nCol)
    !________________________________________________________
    ! 1. gather scalar information such as is,js,mb,nb       |

    ! no idea why this doesn't work
    ! call mpi_allgatherv(mpi_in_place, 0, mpi_datatype_null, rsdm%phi01, cntCol, dispCol, mpi_type_mat, mpi_comm_world, ierr)

    nTotal = rsdm%nk * rsdm%nLead
    is = dispCol(iCol+1)+1
    ie = dispCol(iCol+1)+cntCol(iCol+1)
    do i=is,ie
      intBuff(1,i) = rsdm%phi01(i)%is
      intBuff(2,i) = rsdm%phi01(i)%js
      intBuff(3,i) = rsdm%phi01(i)%mb
      intBuff(4,i) = rsdm%phi01(i)%nb
    enddo
    call mpi_allgatherv(mpi_in_place, 0, mpi_datatype_null, intBuff, cntCol*nVar, dispCol*nVar, mpi_integer, mpi_comm_world, ierr)
    !________________________________________________________|
    ! 2. count the sizes of matrices &  pack into 1D array   |

    ! the difficulty is that each phi01 matrix can be of different size =_=
    ! so the simplest way is to pack the matrices into a 1D array first
    ind(1) = 1
    do ic=1,nCol
      cntMat(ic) = 0
      do i=dispCol(ic)+1,dispCol(ic)+cntCol(ic)
        n = intBuff(3,i)*intBuff(4,i)
        cntMat(ic) = cntMat(ic) + n 
        ind(i+1)   = ind(i) + n
      enddo
    enddo
    dispMat = 0
    do i=2,nCol
      dispMat(i) = dispMat(i-1) + cntMat(i-1)
    enddo
    !write(6,*)'cntMat : ',cntMat
    !write(6,*)'dispMat: ',dispMat
    !write(6,*)'ind    : ',ind
    allocate(matBuff(ind(nTotal+1)-1))
    do i=is,ie
      n = ind(i+1) - ind(i)
      matBuff(ind(i):ind(i+1)-1) = reshape( rsdm%phi01(i)%mat, (/n/) ) 
    enddo
    !________________________________________________________|
    ! 3. allgather the 1D array matBuff                      |
    call mpi_allgatherv(mpi_in_place, 0, mpi_datatype_null, matBuff, cntMat, dispMat, mpi_double_complex, mpi_comm_world, ierr)
    !________________________________________________________|
    ! 4. unpack 1D array to phi01                            |
    do ic=1,nCol
      if(ic==iCol+1)cycle !< skip those own by myself
      do i=dispCol(ic)+1,dispCol(ic)+cntCol(ic)
        mb = intBuff(3,i)
        nb = intBuff(4,i)
        matPointer(1:mb,1:nb) => matBuff(ind(i):ind(i+1)-1) ! reshape 1D back to 2D by pointer
        call MatFromArray(rsdm%phi01(i), matPointer, mb, nb, intBuff(1,i), intBuff(2,i))
      enddo
    enddo
    !________________________________________________________|
    deallocate(matBuff)
  end subroutine


  !> @brief Broadcast tier1 
  subroutine Tier1_Bcast(rsdm, group)
    type(RSDM_T) :: rsdm
    character(len=*),intent(in) :: group

    integer ::  comm, root
    !
    integer,parameter :: nVar = 4
    integer :: i,n,ic,mb,nb, dc, cc, ierr
    complex*16,pointer :: matBuff(:), matPointer(:,:) 
    integer,allocatable :: ind(:), intBuff(:,:)
    !
    logical :: debug = .true.
    !________________________________________________________
    ! 1. bcast  scalar information such as is,js,mb,nb       |

    ! no idea why this doesn't work
    ! call mpi_allgatherv(mpi_in_place, 0, mpi_datatype_null, rsdm%phi01, cntCol, dispCol, mpi_type_mat, mpi_comm_world, ierr)

    if(trim(group) == "row")then
      comm = mpi_comm_row
      root = 0
      dc = dispRow(iRow+1)
      cc = cntRow(iRow+1)
    else 
      comm = mpi_comm_col
      root = 0
      dc = dispCol(iCol+1)
      cc = cntCol(iCol+1)
    endif

    if(debug) write(6,*)'dc, cc = ',dc,cc
    allocate(intBuff(nVar, cc), ind(cc+1))
    do i=1,cc
      intBuff(1,i) = rsdm%phi01(i+dc)%is
      intBuff(2,i) = rsdm%phi01(i+dc)%js
      intBuff(3,i) = rsdm%phi01(i+dc)%mb
      intBuff(4,i) = rsdm%phi01(i+dc)%nb
    enddo
    if(debug)then
      do i=1,cc
        write(6,'(A,5I0)')'intBuff:',i,intBuff(1:4,i)
      enddo
    endif
    call mpi_bcast(intBuff, cc*nVar, mpi_integer, root, comm, ierr)
    if(debug)then
      do i=1,cc
        write(6,'(A,5I0)')'intBuff:',i,intBuff(1:4,i)
      enddo
    endif
    !________________________________________________________|
    ! 2. pack matrices into 1D array                         |

    ! the difficulty is that each phi01 matrix can be of different size =_=
    ! so the simplest way is to pack the matrices into a 1D array first
    ind(1) = 1
    do i=1,cc
      n = intBuff(3,i)*intBuff(4,i)
      ind(i+1)   = ind(i) + n
    enddo
    if(debug)write(6,*)'ind    : ',ind
    allocate(matBuff(ind(cc+1)-1))
    do i=1,cc
      n = ind(i+1) - ind(i)
      matBuff(ind(i):ind(i+1)-1) = reshape( rsdm%phi01(i+dc)%mat, (/n/) ) 
    enddo
    !________________________________________________________|
    ! 3. bcast the 1D array matBuff                          |
    call mpi_bcast(matBuff, ind(cc+1)-1, mpi_double_complex, root, comm, ierr)
    !________________________________________________________|
    ! 4. unpack 1D array to phi01                            |
    do i=1,cc
      mb = intBuff(3,i)
      nb = intBuff(4,i)
      matPointer(1:mb,1:nb) => matBuff(ind(i):ind(i+1)-1) ! reshape 1D back to 2D by pointer
      call MatFromArray(rsdm%phi01(i+dc), matPointer, mb, nb, intBuff(1,i), intBuff(2,i))
    enddo
    !________________________________________________________|
    deallocate(matBuff, intBuff, ind)
    
  end subroutine

  !> @brief AllReduce a mat (e.g. comm = mpi_comm_diagonal)
  subroutine Mat_AllReduce(mat, comm)
    type(mat_T) :: mat
    integer,intent(in) ::  comm
    !
    integer,parameter :: nVar = 4
    integer :: intBuff(nVar)
    complex*16,allocatable :: matBuff(:,:)
    integer :: is,js,mb,nb
    intBuff(1) = mat%is 
    intBuff(2) = mat%js 
    intBuff(3) = mat%is+mat%mb 
    intBuff(4) = mat%js+mat%nb
    ! find the matrix size after addition (min and max of is,js,ie,je)
    call mpi_allreduce(mpi_in_place, intBuff(1:2), 2, mpi_integer, mpi_min, comm, ierr)
    call mpi_allreduce(mpi_in_place, intBuff(3:4), 2, mpi_integer, mpi_max, comm, ierr)
    is = intBuff(1)
    js = intBuff(2)
    mb = intBuff(3)-intBuff(1)
    nb = intBuff(4)-intBuff(2)

    allocate(matBuff(mb,nb))
    call MatToArray(mat,matBuff,mb,nb,is,js)
    call mpi_allreduce(mpi_in_place, matBuff, mb*nb, mpi_double_complex, mpi_sum, comm, ierr)
    call MatFromArray(mat, matBuff, mb, nb, is, js)
    deallocate(matBuff)
  end subroutine 

  !> @brief Reduce a mat
  subroutine Mat_Reduce(mat, comm, root)
    type(mat_T) :: mat
    integer,intent(in) ::  comm, root
    !
    integer,parameter :: nVar = 4
    integer :: intBuff(nVar), intSend(nVar)
    complex*16,allocatable :: matBuff(:,:), matSend(:,:)
    integer :: is,js,mb,nb,ierr
    intSend(1) = mat%is 
    intSend(2) = mat%js 
    intSend(3) = mat%is+mat%mb 
    intSend(4) = mat%js+mat%nb
    ! find the matrix size after addition (min and max of is,js,ie,je)
    call mpi_reduce(intSend(1:2), intBuff(1:2), 2, mpi_integer, mpi_min, root, comm, ierr)
    call mpi_reduce(intSend(3:4), intBuff(3:4), 2, mpi_integer, mpi_max, root, comm, ierr)
    is = intBuff(1)
    js = intBuff(2)
    mb = intBuff(3)-intBuff(1)
    nb = intBuff(4)-intBuff(2)

    write(6,*),'Mat_Reduce: mb, nb = ',mb,nb

    allocate(matSend(mb,nb),matBuff(mb,nb))
    call MatToArray(mat,matSend,mb,nb,is,js)
    call mpi_reduce(matSend, matBuff, mb*nb, mpi_double_complex, mpi_sum, root, comm, ier)
    call MatFromArray(mat, matBuff, mb, nb, is, js)
    deallocate(matBuff,matSend)
  end subroutine 
end module 
