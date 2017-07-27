module mat_mod
#ifdef __MPI__ 
use parallel
#endif
  implicit none
  !________________________________________________________|
  ! Types                                                  |

  !Some notes:
  !  1. Matrices with zero dimension are assumed to be zero matrices
  !     i.e. mb=nb=0
  !  2. Zero matrices will be automatically deallocated
  !  3. Matrices deallocated manually will be assigned mb=nb=-1
  !  4. The only case we may want to have allocated zero matrix is initialization.
  !     In this case, call MatSetZero(mat) or MatFromScalar(mat,0)
  !  5. MatAllocated=F => MatIsZero=T but reverse is not true

  type :: mat_T
!   complex*16,allocatable,dimension(:,:) :: mat
    integer :: mb,nb  !dimension of row and col
    integer :: is,js
    !for mpi
    complex*16,allocatable :: mem(:)
    integer,allocatable    :: desc(:)
    integer :: nsize  !mem size
    integer :: mrow !num of row calculated from NUMROC
    integer :: mcol !num of col calculated from NUMROC
    integer :: iserial,jserial
  end type mat_T
  
  !Allocation is unnecessary unless you want to construct Mat yourself. Use constructors if possible
  public allocateMat,deallocateMat
  public MatCopy                     !Copy one Mat to another Mat
  public MatFromArray!,MatFromScalar   !Construct Matrix from 2D Array, Scalar 
  public MatZero,MatSetZero           !Set Mat to Zero: 
  public MatWrite
  
  !  -  MatZero will deallocate Mat (save memory)  <- normally, use this
  !  -  MatSetZero will construct a zero-filled Mat

  !Setters
  public MatStartAt           !Set mat_T%is, mat_T%js

  !Getters
!  public MatTr                !Matrix Trace
  public MatToArray           !Convert Mat to 2D Array
  public MatMemUsed           !Calculate the Memory used by the Matrix
  !Printers 
!  public MatPrint             !Print Matrix

  !Matrix Operations
!  public MatCTranspose
  public MatMultiply,MatAdd,MatCTransposeAdd!,MatDivide,MatAddTo

!  public MatMaxNorm

  !Internal use
  private MatResize 

  interface MatFromArray  !Construct Mat(MPI form) from 2D Array
    module procedure MatFromArrayD,MatFromArrayZ!,MatFromArrayD2,MatFromArrayZ2
  end interface
  interface MatToArray    !Convert Mat(MPI form) to 2D Array
    module procedure MatToArrayD,MatToArrayZ
  end interface
  ! interface MatFromScalar !Construct Mat from scalar (scalar times identity matrix)
  !   module procedure MatFromScalarI,MatFromScalarD,MatFromScalarZ,MatFromScalarZ0,MatFromScalarZ1
  ! end interface
  interface MatAdd        !Add two Matrices to give another matrix
    module procedure MatAdd,MatAddD,MatAddZ
  end interface
  ! interface MatAddTo      !Add one matrices onto another
  !   module procedure MatAddTo,MatAddTo2,MatAddToI
  ! end interface
  interface MatResize     !Resize Mat, only used internally
    module procedure MatResize,MatResizeMat
  end interface
  interface MatMultiply   !Multiply two matrices
    module procedure MatMultiply,MatMultiplyScalarD,MatMultiplyScalar!,MatMultiplyScalarI
  end interface
  ! interface MatDivide     !Divide one matrix by another
  !   module procedure MatDivideScalar
  ! end interface
  ! interface MatCTranspose !Matrix Conjugate Transpose
  !   module procedure MatCTranspose,MatCTranspose1
  ! end interface
  ! interface MatWrite !Write matrix
  !   module procedure MatWrite,MatWriteUnformat
  ! end interface


  !Methods
  contains
  !========================================================|
  ! Basic Utilities                                        |
  subroutine allocateMat(mat,mb,nb,istat)
    implicit none
    type(mat_T) :: mat
    integer,intent(in) :: mb,nb 
    integer,intent(out) :: istat
    !
    complex*16 :: czero=(0d0,0d0)
    integer :: mp,info
    integer :: mrow,mcol,rDPN,cDPN
    integer :: i,j,k
    integer :: NUMROC
    external:: NUMROC

   if(allocated(mat%mem))then
     write(6,*)'mat already allocated!',allocated(mat%mem)
     stop
   endif
   if(mb==0 .or. nb==0)then
     write(6,*)'mb or nb=0 in allocateMat',mb,nb
     stop
   endif
    
    !calculate memory needed
    !------------------------------
    call getDPN(mb,nRow,rDPN)
    call getDPN(nb,nCol,cDPN)
    rDPN=min(rDPN,cDPN)  !ensure MB=NB
    cDPN=rDPN
    mrow = max(1,NUMROC(mb, rDPN, iRow, 0, nRow ))
    mcol = max(1,NUMROC(nb, cDPN, iCol, 0, nCol ))
    mat%nsize=mcol*mrow
    mat%mrow=mrow
    mat%mcol=mcol
    mat%mb=mb
    mat%nb=nb
    mat%iserial=0
    mat%jserial=0
    call MatStartAt(mat,0,0)
    !------------------------------
    !allocate
    !------------------------------
    istat=0
    allocate(mat%desc(9), STAT=istat)
    if(istat/=0)then
      write(6,*)'allocate mat%desc failed',istat
      stop
    endif
    mat%desc(1:9)=0

    istat=0
    allocate(mat%mem(mat%nsize), STAT=istat)
    if(istat/=0)then
      write(6,*)'allocate mat%mem failed',istat
      stop
    endif
    mat%mem(1:mat%nsize)=czero
    !------------------------------
    !init descriptor
    !------------------------------
     ! write(6,*) 'mb',mat%mb,'nb',mat%nb
     ! write(6,*) 'rDPN',rDPN,'cDPN',cDPN
     ! write(6,*) 'ictxt',myICTXT
     ! write(6,*) 'mrow',mat%mrow,'mcol',mat%mcol
    call DESCINIT(mat%desc, mat%mb, mat%nb, rDPN, cDPN, 0, 0, myICTXT, mat%mrow, info )
    ! write(6,*)'desc',(mat%desc(i),i=1,9)
    !  write(6,*) 'mb',mat%mb,'nb',mat%nb
    !  write(6,*) 'rDPN',rDPN,'cDPN',cDPN
    !  write(6,*) 'ictxt',myICTXT
    !  write(6,*) 'mrow',mat%mrow,'mcol',mat%mcol
  end subroutine
  
  subroutine deallocateMat(mat,istat)
    implicit none
    type(mat_T) :: mat
    !
    integer :: istat

    istat=0
    if(allocated(mat%mem))deallocate(mat%mem,STAT=istat)
    if(istat/=0)then
      write(6,*)'deallocate mat%mem failed'
    endif
    istat=0
    if(allocated(mat%desc))deallocate(mat%desc,STAT=istat)
    if(istat/=0)then
      write(6,*)'deallocate mat%desc failed'
    endif
    mat%mb=0
    mat%nb=0
  end subroutine
  function MatAllocated(mat)
    implicit none
    type(mat_T) :: mat
    logical :: MatAllocated
    !
    logical ltmp1
    ltmp1 = allocated(mat%mem) .and. (mat%mb>0 .and. mat%nb>0) .and. (mat%mrow>0 .and. mat%mcol >0)
    MatAllocated = ltmp1
  end function
  subroutine MatStartAt(mat,is,js)
    implicit none
    type(mat_T) :: mat
    integer,intent(in) :: is,js
    mat%is=is
    mat%js=js
  end subroutine 
  subroutine MatCopy(mat,matSrc)
    implicit none
    type(mat_T) :: mat,matSrc
    !
    call MatCopy2(mat,(1d0,0d0),matSrc)
  end subroutine
  subroutine MatCopy2(mat,alpha,matSrc)
    implicit none
    type(mat_T) :: mat,matSrc
    complex*16 :: alpha
    integer :: is,js,mb,nb,istat
    mb=matSrc%mb
    nb=matSrc%nb
    is=matSrc%is
    js=matSrc%js
    if( .not. MatAllocated(mat))then
      call allocateMat(mat,mb,nb,istat)
      call MatStartAt(mat,is,js)
      call MatSetZero(mat)
    else
      call deallocateMat(mat,istat)
      call allocateMat(mat,mb,nb,istat)
      call MatStartAt(mat,is,js)
      call MatSetZero(mat)
    endif
   !!sub( Y ) := sub( X )
   ! call pzcopy(matSrc%mb*matSrc%nb, matSrc%mem, 1, 1, matSrc%desc, 1, mat%mem, 1, 1, mat%desc, 1 )
   !!sub( Y ) := sub( Y ) + alpha * sub( X ),
   !call pzaxpy(matSrc%mb*matSrc%nb,alpha,matSrc%mem,1,1,matSrc%desc,1,mat%mem,1,1,mat%desc,1)
    CALL PZLACPY( 'All', mb, nb, matSrc%mem, 1, 1, matSrc%desc,mat%mem, 1, 1, mat%desc) 
    mat%mem=mat%mem*alpha
   ! mat%desc=matSrc%desc
    mat%mrow=matSrc%mrow
    mat%mcol=matSrc%mcol
    mat%iserial=matSrc%iserial
    mat%jserial=matSrc%jserial
  end subroutine
  subroutine MatResize(mat,is,js,ie,je)
    implicit none
    type(mat_T) :: mat
    integer,intent(in) :: is,js,ie,je
    !
    logical :: ltmp1
    integer :: i,j,mb,nb,mb2,nb2,is21,js21,ii,jj,istat
    integer :: mrow,mcol,nsize
    complex*16,allocatable,dimension(:,:) :: zmat1
    mb=mat%mb
    nb=mat%nb
    mb2=ie-is!matSrc%mb
    nb2=je-js!matSrc%nb
    is21=mat%is-is
    js21=mat%js-js
    if(mb2==0 .or. nb2==0)then
      call MatZero(mat)
      call MatStartAt(mat,is,js)
    else if( .not. MatAllocated(mat))then
      call deallocateMat(mat,istat)
      call allocateMat(mat,mb2,nb2,istat)
      call MatStartAt(mat,is,js)
      call MatSetZero(mat)
    else
      ! not the case in block
      ! mrow=mat%mrow
      ! mcol=mat%mcol
      ! nsize=mrow*mcol
      ! allocate(zmat1(nsize),STAT=istat)
      ! zmat1(1:nsize)=mat%mat(1:nsize)
      ! call deallocateMat(mat,istat)
      ! call allocateMat(mat,mb2,nb2,istat)
      ! call MatStartAt(mat,is,js)
      call MatSetZero(mat)
      ! do j=1,nb
      !   do i=1,mb
      !     ii=i+is21 
      !     jj=j+js21 
      !     if(ii>0 .and. jj>0 .and. ii<=mb2 .and. jj<=nb2)then
      !       mat%mat(ii,jj)=zmat1(i,j)
      !       !mat%mat(i+ii,j+jj)=zmat1(i,j)
      !     endif
      !   enddo
      ! enddo
      ! deallocate(zmat1,STAT=istat)
    endif
  end subroutine
  subroutine MatResizeMat(mat,matRef)
    implicit none
    type(mat_T) :: mat,matRef
    integer :: is,js,ie,je
    !
    is=matRef%is
    js=matRef%js
    ie=matRef%is+matRef%mb
    je=matRef%js+matRef%nb
    call MatResize(mat,is,js,ie,je)
  end subroutine

  ! function MatEle(mat,i,j)
  !   implicit none
  !   type(mat_T) :: mat
  !   integer :: i,j
  !   complex*16 :: MatEle
  !   !
  !   complex*16 :: czero=(0d0,0d0)
  !   integer :: is,js,ie,je
  !   if(MatIsZero(mat))then
  !     MatEle=czero
  !   else if(MatAllocated(mat))then
  !     is=mat%is
  !     js=mat%js
  !     ie=mat%is+mat%mb 
  !     je=mat%js+mat%nb
  !     if(i>is .and. i<=ie .and. j>js .and. j<=je)then
  !       MatEle=mat%mat(i-is,j-js)
  !     else
  !       MatEle=czero
  !     endif
  !   else
  !     write(6,*)'MatEle ERROR: Mat not allocated'
  !     MatEle=czero
  !   endif
  ! end function
  ! function MatLastEle(mat,i,j)
  !   implicit none
  !   type(mat_T) :: mat
  !   integer :: i,j
  !   complex*16 :: MatLastEle
  !   !
  !   MatLastEle=MatEle(mat,mat%mb-i+1,mat%nb-j+1)
  ! end function

  ! subroutine MatFromScalarI(mat,Ivalue,matRef)
  !   implicit none
  !   type(mat_T) :: mat
  !   type(mat_T) :: matRef
  !   integer :: Ivalue
  !   !
  !   call MatFromScalarZ(mat,dcmplx(1d0*Ivalue,0d0),matRef)
  ! end subroutine
  ! subroutine MatFromScalarD(mat,Rvalue,matRef)
  !   implicit none
  !   type(mat_T) :: mat
  !   type(mat_T) :: matRef
  !   real*8 :: Rvalue
  !   !
  !   call MatFromScalarZ(mat,dcmplx(Rvalue,0d0),matRef)
  ! end subroutine
  ! subroutine MatFromScalarZ0(mat,value,mb,nb,is,js)
  !   implicit none
  !   type(mat_T) :: mat
  !   complex*16,intent(in) :: value
  !   integer,intent(in)    :: mb,nb,is,js
  !   !
  !   complex*16,parameter :: czero=(0d0,0d0)
  !   complex*16,allocatable::zmat(:,:)
  !   integer :: ishift,i,j,n
  !   call MatResize(mat,is,js,is+mb,js+nb)
  !   allocate(zmat(mb,nb))
  !   zmat=czero  
  !   if(value/=czero)then
  !     n=min(mb,nb)
  !     do i=1,n
  !       zmat(i,i)=value
  !     enddo
  !   endif
  !   MatFromArrayZ(mat,zmat,mb,nb)
  ! end subroutine
  ! subroutine MatFromScalarZ(mat,value,matRef)
  !   implicit none
  !   type(mat_T) :: mat
  !   type(mat_T) :: matRef
  !   complex*16,intent(in) :: value
  !   !
  !   complex*16,parameter :: czero=(0d0,0d0)
  !   integer :: ishift,i,j,n
  !   call MatFromScalarZ0(mat,value,matRef%mb,matRef%nb,matRef%is,matRef%js)
  ! end subroutine
  ! subroutine MatFromScalarZ1(mat,value,mb,nb)
  !   implicit none
  !   type(mat_T) :: mat
  !   complex*16,intent(in) :: value
  !   integer,intent(in)    :: mb,nb
  !   !
  !   call MatFromScalarZ0(mat,value,mb,nb,0,0)
  ! end subroutine

  subroutine MatFromArrayD(mat,dmat,mb,nb)
    implicit none
    type(mat_T) :: mat
    real*8,intent(in) :: dmat(mb,nb)
    complex*16,allocatable:: zmat(:,:)
    integer,intent(in) ::  mb,nb
    !working variables
    integer :: ndim,i,j,istat
    complex*8,allocatable::zwork(:) 
    
    if( .not. MatAllocated(mat))then
        call allocateMat(mat,mb,nb,istat)
        call MatStartAt(mat,0,0)
    else
        call deallocateMat(mat,istat)
        call allocateMat(mat,mb,nb,istat)
        call MatStartAt(mat,0,0)
    endif
    !------------------------------
    !call MatFromArrayD2(mat,zmat,mb,nb,0,0)
    !------------------------------
    ndim =nb
    allocate(zwork(ndim**2),zmat(mb,nb))
    zmat(1:mb,1:nb)=dcmplx(dmat(1:mb,1:nb),0.d0)
    call mpiDistributeZ(ndim,zmat,mat%mem,mat%desc,mat%iserial,mat%jserial,zwork)
    deallocate(zwork,zmat)
  end subroutine
  subroutine MatFromArrayZ(mat,zmat,mb,nb)
    implicit none
    type(mat_T) :: mat
    complex*16,intent(in) :: zmat(mb,nb)
    integer,intent(in) ::  mb,nb
    
    !working variables
    integer :: ndim,istat
    complex*16,allocatable::zwork(:) 
    
    if( .not. MatAllocated(mat))then
        call allocateMat(mat,mb,nb,istat)
        call MatStartAt(mat,0,0)
    else
        call deallocateMat(mat,istat)
        call allocateMat(mat,mb,nb,istat)
        call MatStartAt(mat,0,0)
    endif
    !------------------------------
    !call MatFromArrayZ2(mat,zmat,mb,nb,0,0)
    !------------------------------
    ndim =nb
    allocate(zwork(ndim**2))
    call mpiDistributeZ(ndim,zmat,mat%mem,mat%desc,mat%iserial,mat%jserial,zwork)
    deallocate(zwork)
  end subroutine
  ! subroutine MatFromArrayD2(mat,dmat,mb,nb,is,js)
  !   implicit none
  !   type(mat_T) :: mat
  !   real*8,intent(in) :: dmat(mb,nb)
  !   integer,intent(in) ::  mb,nb,is,js
  !   !
  !   integer :: istat
  !   call deallocateMat(mat,istat)
  !   call allocateMat(mat,mb,nb,istat)
  !   call MatStartAt(mat,is,js)
  !   mat%mat(1:mb,1:nb)=dcmplx(dmat(1:mb,1:nb),0d0)
  ! end subroutine
  ! subroutine MatFromArrayZ2(mat,zmat,mb,nb,is,js)
  !   implicit none
  !   type(mat_T) :: mat
  !   complex*16,intent(in) :: zmat(mb,nb)
  !   integer,intent(in) ::  mb,nb,is,js
  !   !
  !   integer :: istat
  !   call deallocateMat(mat,istat)
  !   call allocateMat(mat,mb,nb,istat)
  !   call MatStartAt(mat,is,js)
  !   mat%mat(1:mb,1:nb)=zmat(1:mb,1:nb)
  ! end subroutine

  ! subroutine MatToArrayD(mat,dmat,mdim,ndim)
  !   implicit none
  !   type(mat_T) :: mat
  !   real*8,intent(out) :: dmat(mdim,ndim)
  !   integer,intent(in) ::  mdim,ndim
  !   !
  !   integer :: mb,nb,is,js
  !   integer :: istat
  !   mb=mat%mb 
  !   nb=mat%nb 
  !   is=mat%is 
  !   js=mat%js 
  !   if(mb==mdim .and. nb==ndim)then
  !     dmat(1:mb,1:nb)=dble( mat%mat(1:mb,1:nb) )
  !   else if(is+mb<=mdim .and. js+nb<=ndim)then
  !     dmat=0d0
  !     dmat(is+1:is+mb,js+1:js+nb)=dble( mat%mat(1:mb,1:nb) )
  !   else  
  !     write(6,*)'MatToArrayD: dimension not match:',mb,nb,'vs',mat%mb,mat%nb
  !     stop
  !   endif
  ! end subroutine
  subroutine MatToArrayD(mat,dmat,mdim,ndim)
    implicit none
    type(mat_T) :: mat
    real*8,intent(out) :: dmat(mdim,ndim)
    complex*16,allocatable::zmat(:,:)
    integer,intent(in) ::  mdim,ndim
    !
    integer :: mb,nb,is,js
    integer :: istat
    ! integer :: iserial,jserial
    complex*16,allocatable::zwork(:) 
    mb=mat%mb 
    nb=mat%nb 
    is=mat%is 
    js=mat%js 
   allocate(zwork(ndim**2),zmat(mdim,ndim))
    zmat(1:mdim,1:ndim)=(0.d0,0.d0)
   call mpiGatherZ(ndim,ndim,mat%mem(1),1,1,mat%desc,mat%iserial,mat%jserial,zmat,zwork)
     ! call mpiGatherZ2(mat%mem,mat%desc, zmat, ndim, mat%iserial,mat%jserial)
   dmat(1:mdim,1:ndim)=dble(zmat(1:mdim,1:ndim))
   deallocate(zwork,zmat)
  end subroutine
  
  ! subroutine MatToArrayZ(mat,zmat,mdim,ndim)
  !   implicit none
  !   type(mat_T) :: mat
  !   complex*16,intent(out) :: zmat(mdim,ndim)
  !   integer,intent(in) ::  mdim,ndim
  !   !
  !   integer :: mb,nb,is,js
  !   integer :: istat
  !   mb=mat%mb 
  !   nb=mat%nb 
  !   is=mat%is 
  !   js=mat%js 
  !   if(mb==mdim .and. nb==ndim)then
  !     zmat(1:mb,1:nb)= mat%mat(1:mb,1:nb) 
  !   else if(is+mb<=mdim .and. js+nb<=ndim)then
  !     zmat=(0d0,0d0)
  !     zmat(is+1:is+mb,js+1:js+nb)= mat%mat(1:mb,1:nb) 
  !   else  
  !     write(6,*)'MatToArrayZ: dimension not match:',mb,nb,'vs',mat%mb,mat%nb
  !     stop
  !   endif
  ! end subroutine
  
  subroutine MatToArrayZ(mat,zmat,mdim,ndim)
    implicit none
    type(mat_T) :: mat
    complex*16,intent(inout) :: zmat(mdim,ndim)
    integer,intent(in) ::  mdim,ndim
    !
    integer :: mb,nb,is,js,i,j
    integer :: istat
    complex*16,allocatable::zwork(:) 
   mb=mat%mb 
   nb=mat%nb 
   is=mat%is 
   js=mat%js 
   allocate(zwork(ndim**2))
     write(6,*) 'ndim',ndim,'iserial',mat%iserial,'jserial',mat%jserial
     write(6,*) 'desc in MatToArrayZ',mat%desc
    ! write(6,*) 'mem in MatToArrayZ',mat%mem
    zmat(1:mdim,1:ndim)=(0.d0,0.d0)
   call mpiGatherZ(ndim,ndim,mat%mem(1),1,1,mat%desc,mat%iserial,mat%jserial,zmat,zwork)
     ! call mpiGatherZ2(mat%mem,mat%desc, zmat, ndim, mat%iserial,mat%jserial)
    write(6,*) 'after gather'
  ! do i=1,ndim
   do j=1,20
     write(6,*) 1,j,dble(zmat(1,j))!,dimag(zmat(i,j))
   enddo
   ! enddo
   deallocate(zwork)
  end subroutine
  subroutine MatZero(mat)
    implicit none
    type(mat_T) :: mat
    integer :: istat
    istat=0
    call deallocateMat(mat,istat)
    mat%mb=0
    mat%nb=0
  end subroutine
  function MatIsZero(mat)
    implicit none
    type(mat_T) :: mat
    logical :: MatIsZero
    !
    logical ltmp1
    ltmp1 = (mat%mb==0 .or. mat%nb==0)
    MatIsZero = ltmp1
  end function
  subroutine MatSetZero(mat)
    implicit none
    type(mat_T) :: mat
    if(MatAllocated(mat))then
      mat%mem=(0d0,0d0)
    else
        write(6,*) 'mat not allocated before set to zero'
        stop
    endif
  end subroutine
!  subroutine MatPrint(id,mode,mat)
!    implicit none
!    type(mat_T) :: mat
!    integer :: id
!    character :: mode
!    !
!    integer :: is,js,mb,nb,i,j,nn,kk,k,k2
!    is=mat%is
!    js=mat%js
!    mb=mat%mb
!    nb=mat%nb
!    write(id,*)'is,js,mb,nb:',is,js,mb,nb
!    if(MatIsZero(mat))return
!    if(MatAllocated(mat))then
!      if(mode/='I')then
!        write(id,*)'Re Part'
!        call MatPrintRe(id,mode,6,mat)
!      endif
!      if(mode/='R')then
!        write(id,*)'Im Part'
!        call MatPrintIm(id,mode,6,mat)
!      endif
!    endif
!    write(id,*)'------------------------'
!  end subroutine
!  subroutine MatPrintRe(id,mode,nprint,mat)
!    implicit none
!    type(mat_T) :: mat
!    integer :: id
!    character :: mode
!    integer :: nprint
!    !
!    integer :: is,js,mb,nb,i,j,nn,kk,k,k2
!    is=mat%is
!    js=mat%js
!    mb=mat%mb
!    nb=mat%nb
!    nn=min(nb,nprint)
!    kk=nb/nn+1
!    do k=1,kk
!      k2=nn*(k-1)
!      if(1+k2>nb)exit
!      write(id,'(A5,20I15)')'#',(i,i=1+k2,min(nn+k2,nb) )
!      do j=1,mb
!        write(id,'(I5,20E15.5)'),j,(dble(mat%mat(j,i)),i=1+k2,min(nn+k2,nb) )
!      enddo
!      write(id,'')
!    enddo
!  end subroutine
!  subroutine MatPrintIm(id,mode,nprint,mat)
!    implicit none
!    type(mat_T) :: mat
!    integer :: id
!    character :: mode
!    integer :: nprint
!    !
!    integer :: is,js,mb,nb,i,j,nn,kk,k,k2
!    is=mat%is
!    js=mat%js
!    mb=mat%mb
!    nb=mat%nb
!    nn=min(nb,nprint)
!    kk=nb/nn+1
!    do k=1,kk
!      k2=nn*(k-1)
!      if(1+k2>nb)exit
!      write(id,'(A5,20I15)')'#',(i,i=1+k2,min(nn+k2,nb) )
!      do j=1,mb
!        write(id,'(I5,20E15.5)'),j,(dimag(mat%mat(j,i)),i=1+k2,min(nn+k2,nb) )
!      enddo
!      write(id,'')
!    enddo
!  end subroutine
!  !========================================================|
!  ! Basic Operations                                       |
!  function MatTr(mat)
!    implicit none
!    type(mat_T) :: mat
!    complex*16 :: MatTr
!    !
!    integer :: i,ii,jj
!
!    MatTr=(0d0,0d0)
!    if(mat%is==mat%js)then
!      do i=1,min( mat%mb,mat%nb )
!        matTr=matTr + mat%mat(i,i)
!      enddo
!      return
!    else if(mat%is>mat%js)then
!      ii=0
!      jj=mat%is-mat%js
!    else
!      ii=mat%js-mat%is
!      jj=0
!    endif
!
!    do i=1,min( mat%mb,mat%nb )
!      if(i+ii<=mat%mb .and. i+jj<=mat%nb)then
!        matTr=matTr + mat%mat(i+ii,i+jj)
!      endif
!    enddo
!    !
!
!  end function
!
!  subroutine MatAddTo(mat1,mat2) !Add mat2 to mat1, i.e. mat1+=mat2
!    implicit none
!    type(mat_T) :: mat1,mat2
!    !
!    call MatAddTo2(mat1,(1d0,0d0),mat2)
!  end subroutine
!  subroutine MatAddToI(mat1,alpha,mat2)
!    implicit none
!    type(mat_T) :: mat1,mat2
!    integer :: alpha
!    call MatAddTo2(mat1,(1d0,0d0)*alpha,mat2)
!  end subroutine 
!  subroutine MatAddTo2(mat1,alpha,mat2) !Add mat2 to mat1, i.e. mat1+=alpha*mat2
!    !mat1 will only be enlarged or remain same dimension
!    implicit none
!    type(mat_T) :: mat1,mat2
!    complex*16 :: alpha
!    !
!    complex*16 :: czero=(0d0,0d0)
!    integer :: i,j,ii,jj,is,js,is2,js2,ie,je,ie2,je2
!    logical :: ltmp1,ltmp2,ltmp3
!    if(MatIsZero(mat2) .or. alpha==czero)return
!    is =mat1%is 
!    js =mat1%js
!    is2=mat2%is 
!    js2=mat2%js
!    ie =mat1%is+mat1%mb 
!    je =mat1%js+mat1%nb
!    ie2=mat2%is+mat2%mb 
!    je2=mat2%js+mat2%nb
!
!    ltmp1 = MatAllocated(mat1)!allocated(mat1%mat) .and. (mat1%mb>0 .and. mat1%nb>0)
!    ltmp2 = (is<=is2 .and. js<=js2) .and. (ie>=ie2 .and. je>=je2)
!    ltmp3 = alpha==(0d0,0d0)
!    if(.not. MatAllocated(mat1))then
!      call MatResize(mat1,mat2)
!    else if(.not. ltmp2)then
!      call MatResize(mat1,min(is,is2),min(js,js2),max(ie,ie2),max(je,je2))
!    endif
!    do j=1,mat2%nb
!      do i=1,mat2%mb
!        ii=i+mat2%is-mat1%is
!        jj=j+mat2%js-mat1%js
!        mat1%mat(ii,jj)=mat1%mat(ii,jj)+alpha*mat2%mat(i,j)
!      enddo
!    enddo 
!  end subroutine
!  subroutine MatAdd(mat1,mat2,mat3)
!    implicit none
!    type(mat_T) :: mat1,mat2,mat3
!    call MatCopy(mat3,mat1)
!    call MatAddTo(mat3,mat2)
!  end subroutine
 subroutine MatMultiplyScalar(mat,value)
   implicit none
   type(mat_T):: mat
   complex*16 :: value
   complex*16 :: czero=(0.d0, 0.d0)
   if(.not. MatIsZero(mat)) then 
       mat%mem=mat%mem*value
    !sub( Y ) := sub( Y ) + alpha * sub( X ),
    !call pzaxpy(matSrc%mb*matSrc%nb,alpha,matSrc%mem,1,1,matSrc%desc,1,mat%mem,1,1,mat%desc,1)

   endif
 end subroutine
 subroutine MatMultiplyScalarD(mat,value)
   implicit none
   type(mat_T) :: mat
   real*8 :: value
   call MatMultiplyScalar(mat,dcmplx(value,0d0))
 end subroutine
! subroutine MatMultiplyScalarI(mat,value)
!   implicit none
!   type(mat_T) :: mat
!   integer :: value
!   call MatMultiplyScalar(mat,dcmplx(value*1d0,0d0))
! end subroutine
! subroutine MatDivideScalar(mat,value)
!   implicit none
!   type(mat_T) :: mat
!   complex*16 :: value
!   if(.not. MatIsZero(mat))mat%mat=mat%mat/value
! end subroutine
  ! subroutine MatMultiply(mat1,mat2,mat3)
  !   ! If (mat1==0 .or. mat2==0)  => mat3 is deallocated
  !   ! Else mat3 dimension always = m1 x n2
  !   implicit none
  !   type(mat_T) :: mat1,mat2,mat3
  !   !
  !   logical :: debug=.false.
  !   complex*16 :: czero=(0d0,0d0)
  !   complex*16 :: cunity=(1d0,0d0)
  !   integer :: i,j,ii,jj
  !   integer :: is,js,is2,js2,ie,je,ie2,je2,ns,ne,istat
  !   integer :: m1,n1,m2,n2,m3,n3,nn
  !   complex*16,allocatable,dimension(:,:) ::  zmat1,zmat2,zmat3
  !   is =mat1%is 
  !   js =mat1%js
  !   is2=mat2%is 
  !   js2=mat2%js
  !   ie =mat1%is+mat1%mb
  !   je =mat1%js+mat1%nb
  !   ie2=mat2%is+mat2%mb 
  !   je2=mat2%js+mat2%nb
  !   m1=mat1%mb
  !   n1=mat1%nb
  !   m2=mat2%mb
  !   n2=mat2%nb
  !   m3=m1
  !   n3=n2
  !   ne=min(je,ie2)
  !   ns=max(js,is2)
  !   nn=ne-ns
  !   if(debug)print*,'  - MatMultiply: je,ie2,js,is2=',je,ie2,js,is2
  !   if(debug)print*,'  - MatMultiply: nn=',nn
  !   if(nn<=0 .or. MatIsZero(mat1) .or. MatIsZero(mat2) )then!result is zero!
  !     if(debug)print*,'  - MatMultiply: mat3 is zero',nn,MatIsZero(mat1),MatIsZero(mat2)
  !     call MatZero(mat3)
  !     !call MatResize(mat3,is,js2,ie,je2)
  !     !mat3%mat=czero
  !     return
  !   else
  !     allocate(zmat1(m3,nn),zmat2(nn,n3),zmat3(m3,n3))
  !     do j=1,nn
  !       do i=1,m3
  !         jj=j+ns-js
  !         if(jj<=n1)zmat1(i,j)=mat1%mat(i,jj)
  !       enddo
  !     enddo
  !     do j=1,n3
  !       do i=1,nn
  !         ii=i+ns-is2
  !         if(ii<=m2)zmat2(i,j)=mat2%mat(ii,j)
  !       enddo
  !     enddo
  !     call zgemm('N','N',m3,n3,nn,cunity,zmat1,m3,zmat2,nn,czero,zmat3,m3)
  !     call MatResize(mat3,is,js2,ie,je2)
  !     mat3%mat(1:m3,1:n3)=zmat3(1:m3,1:n3)
  !     deallocate(zmat1,zmat2,zmat3)
  !   endif
  ! end subroutine
 
  !MatMultiplyMPI for matrix multiplication MPI,added by HY
  subroutine MatMultiply(mat1,mat2,mat3)
    ! If (mat1==0 .or. mat2==0)  => mat3 is deallocated
    ! Else mat3 dimension always = m1 x n2
    implicit none
    type(mat_T) :: mat1,mat2,mat3
    !
    logical :: debug=.false.
    complex*16 :: czero=(0d0,0d0)
    complex*16 :: cunity=(1d0,0d0)
    integer :: i,j,ii,jj
    integer :: is,js,is2,js2,ie,je,ie2,je2,ns,ne,istat
    integer :: m1,n1,m2,n2,m3,n3,nn
    complex*16,allocatable,dimension(:,:) ::  zmat1,zmat2,zmat3

      if( .not. MatAllocated(mat3))then
         call allocateMat(mat3,mat1%mb,mat2%nb,istat)
      endif

      call pzgemm('N','N',mat1%mrow, mat2%mrow, mat3%mrow, cunity, &
               mat1%mem , 1,1, mat1%desc,       &
                       mat2%mem, 1,1, mat2%desc,      &
                               czero, mat3%mem,1,1,mat3%desc)
  end subroutine
   
  subroutine MatCTransposeAdd(beta,mat1,alpha,mat2)
    implicit none
    complex*16 ::beta, alpha
    type(mat_T),intent(in)   :: mat1
    type(mat_T),intent(inOut):: mat2

  !i\dot{\sigma}=[h,\sigma]
  !sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
  !?mat1 = mat2 = h\sigma           
  !                                              alpha      A                        beta        C 
    call pzgeadd("C", mat1%mrow, mat1%mcol, alpha, mat1%mem, 1, 1, mat1%desc,  beta, mat2%mem, 1, 1, mat2%desc) 
  end subroutine
  
  subroutine MatAdd(mat1,mat2)
    implicit none
    type(mat_T) :: mat1,mat2
    complex*16 :: alpha,beta
    complex*16 :: cunity=(1d0,0d0)
    if( MatAllocated(mat1) .and. MatAllocated(mat2)) then
    call pzgeadd("N", mat2%mrow, mat2%mcol,cunity, mat2%mem, 1, 1, mat2%desc,cunity , mat1%mem, 1, 1, mat1%desc) 
    !mat1%mem=mat1%mem+mat2%mem
    else 
        write(6,*) 'mat1/mat2 not allcated before addition'
        stop
    endif
  end subroutine
  subroutine MatAddD(alpha,mat1,beta,mat2,mat3)
    implicit none
    type(mat_T) :: mat1,mat2,mat3
    real*8 :: alpha,beta
    complex*16 :: zalpha,zbeta
    integer::istat
  !i\dot{\sigma}=[h,\sigma]
  !sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
  !?mat1 = mat2 = h\sigma
    !mat3=alpha*mat1+beta*mat2
    zalpha=dcmplx(alpha,0d0)
    zbeta=dcmplx(beta,0d0)
    ! call MatCopy2(mat3,zalpha,mat1)
    ! call MatAddZ(zalpha,mat1,zbeta,mat2,mat3)
    if( .not. MatAllocated(mat3))then
       call allocateMat(mat3,mat1%mb,mat2%nb,istat)
    endif
    mat3%mem=zalpha*mat1%mem+zbeta*mat2%mem
  end subroutine
  subroutine MatAddZ(alpha,mat1,beta,mat2,mat3)
    implicit none
    type(mat_T) :: mat1,mat2,mat3
    complex*16 :: alpha,beta
    complex*16 :: cunity=(1d0,0d0)
    integer::istat
    ! call MatCopy2(mat3,alpha,mat1)
    ! call pzgeadd("N", mat2%mrow, mat2%mcol,cunity, mat2%mem, 1, 1, mat2%desc,cunity , mat3%mem, 1, 1, mat3%desc) 
    if( .not. MatAllocated(mat3))then
       call allocateMat(mat3,mat1%mb,mat2%nb,istat)
    endif
    mat3%mem=alpha*mat1%mem+beta*mat2%mem
  end subroutine
  !
!  subroutine MatCTranspose(mat2,mat1)!transpose mat1 to mat2
!    ! If (mat1==0 .or. mat2==0)  => mat3 is deallocated
!    ! Else mat3 dimension always = m1 x n2
!    implicit none
!    type(mat_T) :: mat1,mat2
!    !
!    integer :: i,j,istat
!    integer :: is,js,mb,nb
!    is =mat1%is 
!    js =mat1%js
!    mb =mat1%mb
!    nb =mat1%nb
!    if(MatIsZero(mat1))then
!      call MatZero(mat2)
!      return  
!    else if(.not. MatAllocated(mat1))then
!      write(6,*),'Mat1 not allocated in MatCTranspose',is,js,mb,nb
!      stop  
!    endif
!    call deallocateMat(mat2,istat)
!    call allocateMat(mat2,nb,mb,istat)
!    do j=1,mb
!      do i=1,nb
!        mat2%mat(i,j)=dconjg(mat1%mat(j,i))
!      enddo
!    enddo
!    mat2%is=js
!    mat2%js=is
!  end subroutine 
!  subroutine MatCTranspose1(mat1)
!    ! If (mat1==0 .or. mat2==0)  => mat3 is deallocated
!    ! Else mat3 dimension always = m1 x n2
!    implicit none
!    type(mat_T) :: mat1
!    !
!    type(mat_T) :: mat2
!    call MatCTranspose(mat2,mat1)
!    call MatCopy(mat1,mat2)
!  end subroutine 
!  subroutine MatMaxNorm(mat,norm)
!    ! Compute the Max Norm (Entrywise infinity norm) of Matrix mat 
!    implicit none
!    type(mat_T) :: mat
!    real*8,intent(out) :: norm
!    !
!    integer :: i,j
!    norm=0d0
!    if(MatIsZero(mat))then
!      return
!    else if(MatAllocated(mat))then
!      do j=1,mat%nb
!        do i=1,mat%mb
!          norm=max( norm,cdabs(mat%mat(i,j)) )
!        enddo
!      enddo
!    else
!      write(6,*)'MatMaxNorm ERROR: Mat not allocated'
!    endif
!  end subroutine 
 subroutine MatWrite(fileID,cFormat,mat)
   ! Write Matrix mat in different format
   integer,intent(in) ::  fileID
   character :: cFormat
   type(mat_T) :: mat
   complex*16,allocatable::zmat(:,:)
   !
   integer :: i,j
   !
   allocate(zmat(mat%mb,mat%nb))
    call MatToArray(mat,zmat,mat%mb,mat%nb)
    ! write(6,*) 'in MatWrite'
    ! do j=1,10
    !    write(6,*) 1,j,dble(zmat(1,j)),dimag(zmat(1,j))
    !  enddo
     ! write(6,*) 'in FockAO.dat'
   if(cFormat /= 'u')then 
     ! sparse format
     write(fileID,*),mat%is,mat%js,mat%mb,mat%nb
     do i=1,mat%mb
       do j=1,mat%nb
        ! if( cdabs(zmat(i,j)) > 0d0 )then
           write(fileID,*),i,j,dble(zmat(i,j))!,dimag(zmat(i,j))
           ! write(6,*),i,j,dble(zmat(i,j)),dimag(zmat(i,j))
        ! endif
       enddo
     enddo
   ! else if(cFormat == 'u')then
   !   ! unformatted
   !   call MatWriteUnformat(fileID,mat)
   endif
   deallocate(zmat)
 end subroutine
!  subroutine MatWriteUnformat(fileID,mat)
!    ! Print out Matrix mat unformattedly
!    integer,intent(in) ::  fileID
!    type(mat_T) :: mat
!    !
!    integer :: i,j
!    !
!    write(fileID),mat%is,mat%js,mat%mb,mat%nb
!    write(fileID),((mat%mat(i,j),i=1,mat%mb),j=1,mat%nb)
!  end subroutine
!  subroutine MatRead(fileID,mat)
!    ! Print out Matrix mat unformattedly
!    integer,intent(in) ::  fileID
!    type(mat_T) :: mat
!    !
!    integer :: i,j
!    !
!    read(fileID),mat%is,mat%js,mat%mb,mat%nb
!    read(fileID),((mat%mat(i,j),i=1,mat%mb),j=1,mat%nb)
!  end subroutine
  subroutine MatMemUsed(mat,memory)
    ! Calculate memory used in storing the Matrix (in bytes)
    type(mat_T) :: mat
    integer :: memory
    !
    if(MatAllocated(mat))then
      memory = mat%mb * mat%nb * 16
    else
      memory = 0
    endif
  end subroutine
  
  subroutine getDPN(ndim,np,DPN)
  !DPN:dimension per node
    integer,intent(in)  :: ndim,np
    integer,intent(out) :: DPN
    integer :: i
    i=mod(ndim,np)
    if(i==0)then
      DPN=(ndim)/(np)
    else
      DPN=(ndim+np-mod(ndim,np))/np
    endif
  end subroutine
end module mat_mod
