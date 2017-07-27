program ConvertMat
implicit none
integer::ndim,nRowIn,nColIn,nRowOut,nColOut
character*20::infileName,outfileName
complex*16,allocatable:: mat(:,:)

ndim=126
nRowIn=2
nColIn=4
nRowOut=2
nColOut=2
allocate(mat(ndim,ndim))
infileName='TAPE.td.restart'
call MatLoad_MPI(infileName,ndim,nRowIn,nColIn,mat)
outfileName='oTAPE.td.restart'
call MatSave_MPI(outfileName,ndim,nRowOut,nColOut,mat)

end program ConvertMat

subroutine MatLoad_MPI(fileName,ndim,nRow,nCol,mat)
  character(*),intent(in) ::fileName
  integer,intent(in)      :: ndim,nRow,nCol
  complex*16,intent(inout):: mat(ndim,ndim)
  character*30		  :: fileNameRank
  character*3 		  :: RankName
  real*8      		  :: tt
  complex*16,allocatable  ::work(:)
  integer      		  :: i,j,k,fileID,icounter,jcounter
  integer      		  :: mrow,mcol,rDPN,cDPN,iRow,iCol,nsize
  integer 		  :: NUMROC
  external		  :: NUMROC
  !
icounter=0
jcounter=0
iRank=-1
do iRow=0,nRow-1
  do iCol=0,nCol-1
    call getDPN(ndim,nRow,rDPN)
    call getDPN(ndim,nCol,cDPN)

    rDPN=min(rDPN,cDPN)  !ensure MB=NB
    cDPN=rDPN
    mrow = max(1,NUMROC(ndim, rDPN, iRow, 0, nRow ))
    mcol = max(1,NUMROC(ndim, cDPN, iCol, 0, nCol ))
    nsize=mcol*mrow
    mcol=mcol
    write(6,*) 'load'
    write(6,*) 'iRow',iRow,'iCol',iCol
    write(6,*) 'mrow',mrow,'mcol',mcol,'nsize',nsize
    iRank=iRank+1
    fileID=233+iRank
    if(iRank>=0 .and. iRank < 10) then
        write(RankName,'(I1)') iRank
    elseif(iRank>=10 .and. iRank <100) then
        write(RankName,'(I2)') iRank
    else
        write(6,*) 'situation not defined'
    endif
    write(6,*)'RankName:', RankName
    fileNameRank=trim(fileName)//trim(RankName)
     write(6,*) 'open ',trim(fileNameRank)
    write(6,*) 'jcounter',jcounter
    write(6,*) 'icounter',icounter
    jcounter = jcounter+mcol

    icounter = icounter+mrow
  enddo
enddo
end subroutine

subroutine MatSave_MPI(fileName,ndim,nRow,nCol,mat)
  ! Write Matrix mat in different format
  character(*),intent(in)	:: fileName
  integer,intent(in)     	:: ndim,nRow,nCol
  complex*16,intent(in)		:: mat(ndim,ndim)
  character*30 			:: fileNameRank
  character*3  			:: RankName
  real*8       	                :: tt
  complex*16,allocatable        :: work(:)
  integer  		        :: mrow,mcol,rDPN,cDPN,iRow,iCol,nsize
  integer  		        :: NUMROC
  external 		        :: NUMROC
  integer  		        :: i,j,k,fileID,icounter,jcounter
  !
icounter=0
jcounter=0
do iRow=0,nRow-1
  do iCol=0,nCol-1
    call getDPN(ndim,nRow,rDPN)
    call getDPN(ndim,nCol,cDPN)

    rDPN=min(rDPN,cDPN)  !ensure MB=NB
    cDPN=rDPN
    mrow = max(1,NUMROC(ndim, rDPN, iRow, 0, nRow ))
    mcol = max(1,NUMROC(ndim, cDPN, iCol, 0, nCol ))
    nsize=mcol*mrow
    mcol=mcol

    write(6,*) 'save'
    write(6,*) 'iRow',iRow,'iCol',iCol
    write(6,*) 'mrow',mrow,'mcol',mcol,'nsize',nsize

    iRank=iCol+ iRow*nRow
    fileID=243+iRank
    if(iRank>=0 .and. iRank < 10) then
        write(RankName,'(I1)') iRank
    elseif(iRank>=10 .and. iRank <100) then
        write(RankName,'(I2)') iRank
    else
        write(6,*) 'situation not defined'
    endif
    write(6,*)'RankName: ',RankName

     fileNameRank=trim(fileName)//trim(RankName)
     write(6,*) 'open ',trim(fileNameRank)
     jcounter = iCol*mcol 
     icounter = iRow*mrow 

 enddo
enddo
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
INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
  INTEGER              IPROC, ISRCPROC, N, NB, NPROCS
  INTEGER              EXTRABLKS, MYDIST, NBLOCKS
  INTRINSIC            MOD
  MYDIST = MOD( NPROCS+IPROC-ISRCPROC, NPROCS )
  NBLOCKS = N / NB
  NUMROC = (NBLOCKS/NPROCS) * NB
  EXTRABLKS = MOD( NBLOCKS, NPROCS )
  IF( MYDIST.LT.EXTRABLKS ) THEN
      NUMROC = NUMROC + NB
  ELSE IF( MYDIST.EQ.EXTRABLKS ) THEN
      NUMROC = NUMROC + MOD( N, NB )
  END IF
  RETURN
  END
