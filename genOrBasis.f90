
!---------------------------------------------------------------- 
!> @brief  Generate the transformation matrix for AO to OR 
!> @author CSG 
!>         Modified on ZQ's EhMD code
!> @note   iOrth = 1, Cholesky decomposition 
!>               = 2, Lowdin's method 
!> 15-JUL-2015        
!---------------------------------------------------------------- 
!!! declare
subroutine genOrBasis(nOrbs, iOrth, overl, Fock_AO, denMat_AO, invOverl, overlU, invOverlU)
  implicit none
  !
  integer,            intent(in)     :: nOrbs
  integer,            intent(in)     :: iOrth
  real*8,             intent(in)     :: overl(nOrbs, nOrbs)
  real*8,             intent(in)     :: Fock_AO(nOrbs, nOrbs)
  real*8,             intent(in)     :: denMat_AO(nOrbs, nOrbs)
  real*8,             intent(out)    :: invOverl(nOrbs, nOrbs)
  real*8,             intent(out)    :: overlU(nOrbs, nOrbs)
  real*8,             intent(out)    :: invOverlU(nOrbs, nOrbs)
!!!end declare
!=============================================================================================
  real*8,  allocatable ::  iOrthMat(:,:)
  integer              ::  i, j, k, l, ii, jj, kk, ll, lwork, istat, info
  integer, allocatable ::  ipiv(:)
  real*8,  allocatable ::  eigval(:), work(:), dtmpmat1(:,:), dtmpmat2(:,:)
 
  real*8 :: fermiDirac
  external fermiDirac 
!---------------------------------------------------------------------------------------------

!__________________________________________________________________________________________|
!      invOverl is the matrix transforming density matrix from AO to OR                   |
!       overlU  is the matrix transforming Fock    matrix from Ao to OR
!     invOverlU is the matrix transforming density matrix from Or to AO
! ========================  Cholesky decompostion Method  ===================================
  if ( iOrth == 0 .or. iOrth == 3) then   ! Cholesky decomposition for orthogonalization

    call dcopy(nOrbs*nOrbs, overl, 1, invOverl, 1)

    call dpotrf('U', nOrbs, invOverl, nOrbs, info)       ! matrix U 
    if (info/=0) then
      write(6,*) 'Cholesky factorization failed in initializeMat.f90'
      stop
    end if
    do ii=2, nOrbs
      do jj=1, ii-1
        invOverl(ii,jj) = 0.d0
      enddo
    enddo

    call dcopy(nOrbs*nOrbs, invOverl, 1, overlU, 1)

    call dtrtri('U', 'N', nOrbs, invOverl, nOrbs, info)   ! matrix U^(-1)
    if (info/=0) then
       write(6,*) 'matrix inversion failed in initializeMat.f90!'
       stop
    end if
    do ii=2, nOrbs
      do jj=1, ii-1
        invOverl(ii,jj) = 0.d0
      enddo
    enddo

    if (iOrth == 0) then
      write(6,*) 'Cholesky factorization finished.'
    else
!======================= MOlecular Orbital Orthogonalization =============================


!----------------------------------------------------------------------------------------
!       V: eigenVectors of Fock_OR;                                                      |
!       Fock_MO   =  V^T * Fock_OR * V    =  V^T * X^T * Fock_AO * X * V                 | 
!         Define: X*V = Q; then                                                          |
!         Fock_MO  = Q^T  * Fock_AO * Q                                                  | 
!  -----------------------------------------------------------------------------------   |
!       denMat_MO =  V^T * denMat_OR * V  =  V^T * U * denMat_AO * U^T * V               | 
!         Define: V^T * U  = K; then                                                     |
!         denMat_MO  = K * denMat_AO * K^T                                               | 
!                                                                                        |
!         Define: {K}^-1  = H; then    i                                                 |
!         denMat_AO  = H * denMat_OR * H^T                                               | 
!  -----------------------------------------------------------------------------------   |
!            Q: invOverl;        K:  overlU;        H:   invOverlU                       |
!----------------------------------------------------------------------------------------
    
  !_______________________________________________________________________________|
  !                 1. Got the eigenVectors of Fock_OR :  V                       |         
      lwork = 3*nOrbs
      allocate(eigval(nOrbs), work(lwork), ipiv(nOrbs), &
               dtmpmat1(nOrbs, nOrbs), dtmpMat2(nOrbs, nOrbs), STAT=istat) 
      if(istat /=0 ) then
        write(6,*) 'error when allocating matrice in MO orthogonalization '
        stop
      endif
    
      dtmpMat1(1:nOrbs, 1:nOrbs) = Fock_AO(1:nOrbs, 1:nOrbs)
      call dtrmm('r', 'u', 'n', 'n', nOrbs, nOrbs, 1.d0, invOverl, nOrbs, dtmpMat1, nOrbs)
      call dtrmm('l', 'u', 't', 'n', nOrbs, nOrbs, 1.d0, invOverl, nOrbs, dtmpMat1, nOrbs)
      !dtmpMat1: Fock_OR

      ! Get the eigenVecters of Fock_OR
      call dsyev('V', 'U', nOrbs, dtmpMat1, nOrbs, eigval,  work, lwork, info)
      if ( info /= 0 ) then
        write(6,*) 'DSYEV fails in MO orthogonalization'
        write(6,*) 'Error info:', info
        stop
      endif

      write(6,*) 'MO eigen energy:'
      call fermi()
      do i=1,nOrbs
        write(6,*) i, eigval(i),fermiDirac(eigVal(i), 300d0)
      enddo

      ! dtmpMat1: V 
  !_______________________________________________________________________________|
  !                 2. Got the tranform matrix {X * V}                            | 

      call dgemm ('N', 'N', nOrbs, nOrbs, nOrbs, 1.d0, invOverl, nOrbs, dtmpMat1, nOrbs, 0.d0, dtmpMat2, nOrbs ) 
      invOverl(1:nOrbs, 1:nOrbs) = dtmpMat2(1:nOrbs, 1:nOrbs)
      ! invOverl = Q : {X * V}
  !_______________________________________________________________________________|
  !                 3. Got the tranform matrix {V^T * U} and {V^T * U}^-1         | 
      call dgemm ('T', 'N', nOrbs, nOrbs, nOrbs, 1.d0, dtmpMat1, nOrbs, overlU, nOrbs, 0.d0, dtmpMat2, nOrbs ) 
      overlU(1:nOrbs, 1:nOrbs) = dtmpMat2(1:nOrbs, 1:nOrbs)
      ! overlU = K : V^T * U
      invOverlU(1:nOrbs, 1:nOrbs) = dtmpMat2(1:nOrbs, 1:nOrbs)
      ! invOverlU : K

      call dgetrf(nOrbs, nOrbs, invOverlU, nOrbs, ipiv,  info)  
      if (info/=0) then
        write(6,*) 'L-U factorization failed in MO orthogonalization'
        write(6,*) 'Error info:', info
        stop
      end if

      call dgetri(nOrbs, invOverlU, nOrbs, ipiv, work, lwork, info)
      if (info/=0) then
        write(6,*) 'inversion of overl failed in MO orthogonalization'
        write(6,*) 'Error info:', info
        stop
      end if
!      do ii=2, nOrbs
!        do jj=1, ii-1
!          invOverlU(ii,jj) = invOverlU(jj,ii)
!        enddo
!      enddo
      ! invOverlU = H :  {V^T * U}^-1
  !_______________________________________________________________________________|

      write(6,*) 'Molecular Orbital Orthogonalization finished.'

!      !debug 
!      write(6,*) '=======  Check Orthogonalization:  ====='
!      dtmpMat1(1:nOrbs, 1:nOrbs) = Fock_AO(1:nOrbs, 1:nOrbs)
!      write(6,*) 'Fock_AO'
!      do i=1, nOrbs
!        do j=1, nOrbs
!          write(6,*) i, j,  dtmpMat1(i,j)
!        enddo
!      enddo
!      call dgemm ('N', 'N', nOrbs, nOrbs, nOrbs, 1.d0, dtmpMat1, nOrbs, invOverl, nOrbs, 0.d0, dtmpMat2, nOrbs ) 
!      call dgemm ('T', 'N', nOrbs, nOrbs, nOrbs, 1.d0, invOverl, nOrbs, dtmpMat2, nOrbs, 0.d0, dtmpMat1, nOrbs ) 
!      write(6,*) 'Fock_MO'
!      do i=1, nOrbs
!        do j=1, nOrbs
!          write(6,*) i, j,  dtmpMat1(i,j)
!        enddo
!      enddo
!       
!      dtmpMat1(1:nOrbs, 1:nOrbs) = denMat_AO(1:nOrbs, 1:nOrbs)
!      write(6,*) 'denMat_AO'
!      do i=1, nOrbs
!        do j=1, nOrbs
!          write(6,*) i, j,  dtmpMat1(i,j)
!        enddo
!      enddo
!
!      call dgemm ('N', 'T', nOrbs, nOrbs, nOrbs, 1.d0, dtmpMat1, nOrbs, overlU, nOrbs, 0.d0, dtmpMat2, nOrbs ) 
!      call dgemm ('N', 'N', nOrbs, nOrbs, nOrbs, 1.d0, overlU, nOrbs, dtmpMat2, nOrbs, 0.d0, dtmpMat1, nOrbs ) 
!      write(6,*) 'denMat_MO'
!      do i=1, nOrbs
!        do j=1, nOrbs
!          write(6,*) i, j,  dtmpMat1(i,j)
!        enddo
!      enddo
!
!      call dgemm ('N', 'T', nOrbs, nOrbs, nOrbs, 1.d0, dtmpMat1, nOrbs, invOverlU, nOrbs, 0.d0, dtmpMat2, nOrbs ) 
!      call dgemm ('N', 'N', nOrbs, nOrbs, nOrbs, 1.d0, invOverlU, nOrbs, dtmpMat2, nOrbs, 0.d0, dtmpMat1, nOrbs ) 
!      write(6,*) 'denMat_AO'
!      do i=1, nOrbs
!        do j=1, nOrbs
!          write(6,*) i, j,  dtmpMat1(i,j)
!        enddo
!      enddo


    endif

! ==================================  Lowdin's Method  ===================================
  else if (iOrth == 2 ) then  ! Lowdin's method for orthogonalization
    write(6,*) 'Lowdin method not available right now'
    stop

    lwork  = 4*nOrbs**2
    allocate( eigval(nOrbs), work(lwork), dtmpmat1(nOrbs,nOrbs), &
              dtmpmat2(nOrbs,nOrbs), STAT=istat )
    if(istat /=0 ) then
      write(6,*) 'error when allocating matrice in Lowdin orthogonalization '
      stop
    endif

    call dcopy(nOrbs*nOrbs, overl, 1, dtmpmat2, 1)
    call DSYEV( 'V', 'U', nOrbs, dtmpmat2, nOrbs, eigval, work, lwork, info )
    if ( info /= 0 ) then
      write(6,*) info
      write(6,*) 'DSYEV fails in Lowdin orthogonalization'
      stop
    endif
  
    do ii=1, nOrbs
      dtmpmat1(1:nOrbs, ii) = dtmpmat2(1:nOrbs, ii) / dsqrt(eigval(ii))
    enddo
    call DGEMM('n', 't', nOrbs, nOrbs, nOrbs, 1.d0, dtmpmat1, nOrbs, dtmpmat2, nOrbs,&
               0.d0, invOverl, nOrbs)
  
    deallocate(eigval, work, dtmpmat1, dtmpmat2) 

  else 
    write(6,*) 'Error! Wrong iOrth method input !'
    stop
  endif


end subroutine genOrBasis
