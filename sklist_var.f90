!!!declare
subroutine SKlist(i,j,skfile)
! THIS SUBROUTINE IS TO POINT TO THE CORRECT SK FILES
! NOT ALL THE ELEMENTS HAVE SK FILES
! i,j  THE PROTON MUMBER OF FOR ATOM PAIR
! skfile THE DIRECTORY AND FILENAMES OF SKFILE
  implicit none

  integer, intent(in) :: i, j
  character*250, intent(out)  :: skfile
!!!end declare
character*80 :: match
integer :: len,istat

call get_environment_variable ('HOME', match, len, istat, .true.)
!!only H C N is changed to 3ob-3-1 for now
  select case(i)
  case(1)
    select case(j)
    case(1)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/3ob-3-1/H-H.skf'
    case(6)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/3ob-3-1/H-C.skf'
    case(7)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/3ob-3-1/H-N.skf'
    case(8)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/H-O.skf'
    case(13)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO_new/H-Al.skf'
    case(14)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/Si/H-Si.skf'
!		case(14)
!			skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO/H-Si.skf'
!   case(15)
!     skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/ph'
    case(16)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/H-S.skf'
    case(31)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/H-Ga.skf'
    case(33)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/H-As.skf'
    case(47)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/H-Ag.skf'
    case(79)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/newskfile.jacky/newSK/H-Au.skf'
    case default
      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
    end select
  case(6)
    select case(j)
    case(1)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/3ob-3-1/C-H.skf'
    case(6)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/3ob-3-1/C-C.skf'
    case(7)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/3ob-3-1/C-N.skf'
    case(8)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/C-O.skf'
    case(13)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO_new/C-Al.skf'
    case(14)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/Si/C-Si.skf'
!	case(14)
!			skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO/C-Si.skf'
!   case(15)
!     skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/cp'
    case(16)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/C-S.skf'
    case(31)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/C-Ga.skf'
    case(47)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/C-Ag.skf'
    case(79)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/newskfile.jacky/newSK/C-Au.skf'
    case default
      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
      stop
    end select
  case(7)
    select case(j)
    case(1)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/3ob-3-1/N-H.skf'
    case(6)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/3ob-3-1/N-C.skf'
    case(7)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/3ob-3-1/N-N.skf'
    case(8)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/N-O.skf'
    case(13)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO_new/N-Al.skf'
    case(14)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/Si/N-Si.skf'
!   case(14)
!			skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO/N-Si.skf'
!   case(15)
!     skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/np'
    case(16)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/N-S.skf'
    case(79)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/newskfile.jacky/newSK/N-Au.skf'
    case default
      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
      stop
    end select
  case(8)
    select case(j)
    case(1)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/O-H.skf'
    case(6)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/O-C.skf'
    case(7)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/O-N.skf'
    case(8)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/O-O.skf'
    case(13)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO_new/O-Al.skf'
    case(14)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/Si/O-Si.skf'
!    case(14)
!			skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO/O-Si.skf'
!   case(15)
!     skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/op'
    case(16)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/O-S.skf'
    case(31)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/O-Ga.skf'
    case(47)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/O-Ag.skf'
    case(79)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/newskfile.jacky/newSK/O-Au.skf'
    case default
      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
      stop
    end select
  case(13)
    select case(j)
    case(1)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO_new/Al-H.skf'
    case(6)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO_new/Al-C.skf'
    case(7)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO_new/Al-N.skf'
    case(8)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO_new/Al-O.skf'
    case(13)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO_new/Al-Al.skf'
    case(14)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO_new/Al-Si.skf'
!		case(14)
!		  skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO/Al-Si.skf'
    case default
      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
      stop
    end select
  case(14)
    select case(j)
    case(1)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/Si/Si-H.skf'
    case(6)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/Si/Si-C.skf'
    case(7)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/Si/Si-N.skf'
    case(8)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/Si/Si-O.skf'
    case(13)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO_new/Si-Al.skf'
    case(14)
!      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/Si/Si-Si.skf'
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Si-Si.skf'
    case(16)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Si-S.skf'
    case(31)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Si-Ga.skf'
    case(33)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Si-As.skf'
    case(47)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Si-Ag.skf'
!    case(1)
!      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO/Si-H.skf'
!    case(6)
!      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO/Si-C.skf'
!    case(7)
!      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO/Si-N.skf'
!    case(8)
!      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO/Si-O.skf'
!    case(13)
!      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO/Si-Al.skf'
!		case(14)
!			skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/AlO/Si-Si.skf'
    case default
      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
      stop
    end select

!  case(15)
!    select case(j)
!   case(1)
!     skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/ph'
!   case(6)
!     skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/pc'
!   case(7)
!     skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/pn'
!   case(8)
!     skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/po'
!   case(15)
!     skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/pp'
!   case(16)
!     skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/ps'
!    case default
!      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
!    end select
  case(16)
    select case(j)
    case(1)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/S-H.skf'
    case(6)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/S-C.skf'
    case(7)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/S-N.skf'
    case(8)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/S-O.skf'
!   case(15)
!     skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/sp'
    case(14)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/S-Si.skf'
    case(16)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/mio-0-1/S-S.skf'
    case(31)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-1/S-Ga.skf'
    case(33)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-1/S-As.skf'
    case(79)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/newskfile.jacky/newSK/S-Au.skf'
    case default
      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
    end select
  case(31)
    select case(j)
    case(1)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ga-H.skf'
    case(6)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ga-C.skf'
    case(8)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ga-O.skf'
    case(14)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ga-Si.skf'
    case(16)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ga-S.skf'
    case(31)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ga-Ga.skf'
    case(33)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ga-As.skf'
    case default
      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
      stop
    end select
  case(33)
    select case(j)
    case(1)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/As-H.skf'
    case(14)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/As-Si.skf'
    case(16)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/As-S.skf'
    case(31)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/As-Ga.skf'
    case(33)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/As-As.skf'
    case default
      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
      stop
    end select
  case(47)
    select case(j)
    case(1)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ag-H.skf'
    case(6)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ag-C.skf'
    case(8)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ag-O.skf'
    case(14)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ag-Si.skf'
    case(47)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/hyb-0-2/Ag-Ag.skf'
    case default
      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
      stop
    end select
  case(79)
    select case(j)
    case(1)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/newskfile.jacky/newSK/Au-H.skf'
    case(6)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/newskfile.jacky/newSK/Au-C.skf'
    case(7)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/newskfile.jacky/newSK/Au-N.skf'
    case(16)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/newskfile.jacky/newSK/Au-S.skf'
    case(79)
      skfile = match(1:len)//'/lodestar/DFTB/trunk/sk/newskfile.jacky/newSK/Au-Au.skf'
    case default
      write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
    end select
  case default
    write(6,*) 'SORRY! ONLY FEW SKFILES ARE AVAILABLE AT THE MOMENT'
  end select
  
end subroutine SKlist
