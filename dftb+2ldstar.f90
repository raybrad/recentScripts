program convert_dftb_plus_lead_2_ldstar_lead
   !=====================================================================!
   !                                                                     !
   ! DFTB+ outputs the shift vectors and charge/orb. Lodestar needs      !
   ! shift/atom and charge/atom. Here we construct                       !
   !                                                                     !
   ! **  shift/atom = shiftPerL(1,:)                                     !
   ! **  charge/atom = sum(chargesSt(:,i)                                !
   !                                                                     !
   ! To complete the construction of the lead needed by lodestar, one    !
   ! needs to input by hand the second line of the file that describes   !
   ! the fermi level and the bot energy after this program.              !
   !                                                                     !
   !=====================================================================!
   real*8, allocatable :: nOrbAtom(:), shiftPerL(:,:),chargesSt(:,:),chargesAtom(:)
   open(11,file='shiftcont_source.dat') ! dftb+ output 
   read(11,*) nAtomSt, mShellSt, mOrbSt, nSpinSt 
   ALLOCATE(nOrbAtom(nAtomSt))
   read(11, *) nOrbAtom(:)
   ALLOCATE(shiftPerL(mShellSt, nAtomSt))
   read(11, *) shiftPerL(:,:)
   ALLOCATE(chargesSt(mOrbSt, nAtomSt))
   read(11, *) chargesSt(:,:)
   close(11)
   write(6,*)  nAtomSt ! how many atoms are included in this lead
   write(6,*)  shiftPerL(1,:) ! shift per atom
   allocate(chargesAtom(nAtomSt))
   do i = 1,nAtomSt
      chargesAtom(i)=sum(chargesSt(:,i))
   end do
   write(6,*) chargesAtom(:)  ! charge per atom
   end

