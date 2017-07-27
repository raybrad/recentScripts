!> @file cscMatrix_op.f90
!> @brief Definition of matrix in CSC format
!>
!> @details This file contain definition of matrix in CSC format,
!>          and related subroutines.
!> 
!> @date 28-Dec-2015
!> @author Yi ZHOU
!------------------------------------------------------------------------
!!!declare
module cscMatrix_op
    implicit none

    ! Public members
    integer, parameter :: DP = kind(1.d0)
    real(kind = DP), parameter :: dOne = 1.d0
    real(kind = DP), parameter :: dZero = 0.d0
    real(kind = DP), parameter :: dmOne = -1.d0
    complex(kind = DP), parameter :: cOne = dcmplx(1.d0)
    complex(kind = DP), parameter :: cZero = dcmplx(0.d0)
    complex(kind = DP), parameter :: cmOne = dcmplx(-1.d0)
    complex(kind = DP), parameter :: cIOne = dcmplx(0, 1.d0)

    type :: dCSCMat
        integer :: row, col, nnz
        real(kind = DP), allocatable :: val(:)
        integer, allocatable :: rowInd(:)
        integer, allocatable :: colPtr(:)
    end type dCSCMat

    type :: zCSCMat
        integer :: row, col, nnz
        complex(kind = DP), allocatable :: val(:)
        integer, allocatable :: rowInd(:)
        integer, allocatable :: colPtr(:)
    end type zCSCMat

    ! Operator overload
    interface assignment (=)
        module procedure assignCSCMat_d_d
        module procedure assignCSCMat_d_z
        module procedure assignCSCMat_z_z
    end interface ! =

    interface operator (+)
        module procedure addCSCMat_d_d
        module procedure addCSCMat_d_z
        module procedure addCSCMat_z_d
        module procedure addCSCMat_z_z
    end interface ! +

    interface operator (-)
        module procedure subCSCMat_d_d
        module procedure subCSCMat_d_z
        module procedure subCSCMat_z_d
        module procedure subCSCMat_z_z
        module procedure negCSCMat_d
        module procedure negCSCMat_z
    end interface ! -

    interface operator (*)
        module procedure xmulCSCMat_d_d
        module procedure xmulCSCMat_d_z
        module procedure xmulCSCMat_z_z
    end interface ! *

    interface operator (.dot.)
        module procedure dotCSCMat_d_d
        module procedure dotCSCMat_z_z
    end interface ! .dot.

    ! Member functions
    interface allocated
        module procedure allocated_d
        module procedure allocated_z
    end interface allocated

    interface initCSCMat
        module procedure initCSCMat_d
        module procedure initCSCMat_z
    end interface initCSCMat

    interface allocateCSCMat
        module procedure allocateCSCMat_d
        module procedure allocateCSCMat_z
    end interface allocateCSCMat

    interface deallocateCSCMat
        module procedure deallocateCSCMat_d
        module procedure deallocateCSCMat_z
    end interface deallocateCSCMat

    interface reallocateCSCMat
        module procedure reallocateCSCMat_d
        module procedure reallocateCSCMat_z
    end interface reallocateCSCMat

    interface outputCSCMat
        module procedure outputCSCMat_d
        module procedure outputCSCMat_z
    end interface outputCSCMat

    interface printCSCMat
        module procedure printCSCMat_d
        module procedure printCSCMat_z
    end interface printCSCMat
    
    interface dns2csc
        module procedure dns2csc_d_d
        module procedure dns2csc_d_z
        module procedure dns2csc_z_z
        module procedure dns2csc_op_z_z
    end interface dns2csc

    interface transpose
        module procedure transCSCMat_d
        module procedure transCSCMat_z
    end interface transpose

    interface conjugate
        module procedure conjCSCMat_z
    end interface conjugate

    interface dble
        module procedure dbleCSCMat_z
    end interface dble

    interface dimag
        module procedure dimagCSCMat_z
    end interface dimag

    interface diag
        module procedure diagCSCMat_d
        module procedure diagCSCMat_z
        module procedure diag_csccsc_d_d
        module procedure diag_csccsc_d_z
        module procedure diag_csccsc_z_d
        module procedure diag_csccsc_z_z
        module procedure diag_cscdns_d_d
        module procedure diag_cscdns_d_z
        module procedure diag_cscdns_z_d
        module procedure diag_cscdns_z_z
    end interface diag

    interface cscmm
        module procedure cscmm_cscdns_d_d
        module procedure cscmm_cscdns_d_z
        module procedure cscmm_cscdns_z_d
        module procedure cscmm_cscdns_z_z
        module procedure cscmm_cscdns_d_d_op
        module procedure cscmm_cscdns_z_z_op

        module procedure cscmm_csccsc_d_d
    end interface cscmm

    ! Implementation of member functions
    contains
    include 'cscMatrix_op/assignCSCMat.inc'
    include 'cscMatrix_op/addCSCMat.inc'
    include 'cscMatrix_op/subCSCMat.inc'
    include 'cscMatrix_op/xmulCSCMat.inc'
    include 'cscMatrix_op/dotCSCMat.inc'

    include 'cscMatrix_op/allocatedCSCMat.inc'
    include 'cscMatrix_op/initCSCMat.inc'
    include 'cscMatrix_op/allocateCSCMat.inc'
    include 'cscMatrix_op/deallocateCSCMat.inc'
    include 'cscMatrix_op/reallocateCSCMat.inc'
    include 'cscMatrix_op/outputCSCMat.inc'
    include 'cscMatrix_op/printCSCMat.inc'

    include 'cscMatrix_op/dns2csc.inc'
    include 'cscMatrix_op/transCSCMat.inc'
    include 'cscMatrix_op/conjCSCMat.inc'
    include 'cscMatrix_op/dbleCSCMat.inc'
    include 'cscMatrix_op/dimagCSCMat.inc'

    include 'cscMatrix_op/diagCSCMat.inc'
    include 'cscMatrix_op/cscmm.inc'

end module cscMatrix_op
!!!end declare

