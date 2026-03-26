!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Computes the left-hand side matrix for a Galerkin projection into the
!!        scalar-form of a vector-valued function space.
!> @details A Galerkin projection solves the equation
!!          gamma_i \sum_j \int( phi_i * phi_j  dx) * c_j
!!          = gamma_i int( phi_i * f  dx)
!!          for the coefficients c_i of the projected field in the target space,
!!          where gamma are the test functions in the target space, phi_i are
!!          the basis functions in the target space, and f is the field to be
!!          projected.
!!          This kernel computes the left-hand side integral, for the case that
!!          phi_i = \sum_j |psi_j|, where psi_j are the vector basis functions
!!          in the target space. This allows a scalar field to be projected into
!!          a scalar form of a vector-valued function space.

module sci_mass_matrix_scalarvec_kernel_mod

  use argument_mod,             only: arg_type, func_type,                     &
                                      GH_REAL, GH_READ, GH_WRITE,              &
                                      GH_FIELD, GH_OPERATOR,                   &
                                      GH_BASIS, GH_DIFF_BASIS,                 &
                                      GH_QUADRATURE_XYoZ, GH_OPERATOR,         &
                                      ANY_SPACE_2, ANY_DISCONTINUOUS_SPACE_3,  &
                                      CELL_COLUMN
  use constants_mod,           only: i_def, r_def
  use fs_continuity_mod,       only: Wchi
  use kernel_mod,              only: kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  type, public, extends(kernel_type) :: mass_matrix_scalarvec_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
        arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, ANY_SPACE_2, ANY_SPACE_2),    &
        arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  Wchi),                        &
        arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)    &
    /)
    type(func_type) :: meta_funcs(2) = (/                                      &
        func_type(ANY_SPACE_2, GH_BASIS),                                      &
        func_type(Wchi,        GH_BASIS, GH_DIFF_BASIS)                        &
    /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: mass_matrix_scalarvec_code
  end type mass_matrix_scalarvec_kernel_type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: mass_matrix_scalarvec_code

contains

  !> @brief Computes the LHS matrix for projecting a scalar into a scalar
  !!        form of a vector-valued function space.
  !> @param[in]     col_idx               Column index
  !> @param[in]     nlayers               Number of layers in mesh
  !> @param[in]     ncell_3d              Num of cells in mesh (this partition)
  !> @param[in,out] projection_lhs        LHS operator to calculate
  !> @param[in]     chi1                  First coordinate field
  !> @param[in]     chi2                  Second coordinate field
  !> @param[in]     chi3                  Third coordinate field
  !> @param[in]     panel_id              Field containing the mesh panel ID
  !> @param[in]     ndf_wvec              Number of DOFs per cell for vector
  !!                                      function space
  !> @param[in]     basis_wvec            Basis functions for the vector space
  !!                                      evaluated at quadrature points
  !> @param[in]     ndf_wchi              Number of DOFs per cell for Wchi
  !> @param[in]     undf_wchi             Num of DOFs in this partition for Wchi
  !> @param[in]     map_wchi              Dofmap for Wchi
  !> @param[in]     basis_wchi            Basis functions for Wchi evaluated at
  !!                                      quadrature points
  !> @param[in]     diff_basis_wchi       Derivative of basis functions for Wchi
  !!                                      evaluated at quadrature points
  !> @param[in]     ndf_pid               Number of DOFs per cell for panel_id
  !> @param[in]     undf_pid              Num of panel_id DOFs in this partition
  !> @param[in]     map_pid               Dofmap for panel_id
  !> @param[in]     nqp_h                 Number of horizontal quadrature points
  !> @param[in]     nqp_v                 Number of vertical quadrature points
  !> @param[in]     wqp_h                 Quadrature weights horizontal
  !> @param[in]     wqp_v                 Quadrature weights vertical
  subroutine mass_matrix_scalarvec_code(col_idx, nlayers, ncell_3d,            &
                                        projection_lhs,                        &
                                        chi1, chi2, chi3,                      &
                                        panel_id,                              &
                                        ndf_wvec,                              &
                                        basis_wvec,                            &
                                        ndf_wchi, undf_wchi, map_wchi,         &
                                        basis_wchi, diff_basis_wchi,           &
                                        ndf_pid, undf_pid, map_pid,            &
                                        nqp_h, nqp_v, wqp_h, wqp_v)

    use sci_coordinate_jacobian_mod, only : coordinate_jacobian

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_3d, col_idx
    integer(kind=i_def), intent(in) :: nqp_h, nqp_v
    integer(kind=i_def), intent(in) :: ndf_wvec, ndf_wchi, ndf_pid
    integer(kind=i_def), intent(in) :: undf_wchi, undf_pid

    integer(kind=i_def), intent(in) :: map_wchi(ndf_wchi)
    integer(kind=i_def), intent(in) :: map_pid(ndf_pid)

    real(kind=r_def),    intent(in) :: basis_wvec(3,ndf_wvec,nqp_h,nqp_v)
    real(kind=r_def),    intent(in) :: basis_wchi(1,ndf_wchi,nqp_h,nqp_v)
    real(kind=r_def),    intent(in) :: diff_basis_wchi(3,ndf_wchi,nqp_h,nqp_v)
    real(kind=r_def),    intent(in) :: wqp_h(nqp_h)
    real(kind=r_def),    intent(in) :: wqp_v(nqp_v)

    real(kind=r_def),    intent(inout) :: projection_lhs(ncell_3d,ndf_wvec,ndf_wvec)
    real(kind=r_def),    intent(in)    :: chi1(undf_wchi)
    real(kind=r_def),    intent(in)    :: chi2(undf_wchi)
    real(kind=r_def),    intent(in)    :: chi3(undf_wchi)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)

    ! Internal variables
    integer(kind=i_def) :: df, ipanel, ik, i, k
    integer(kind=i_def) :: df1, df2, qph, qpv
    real(kind=r_def)    :: chi1_e(ndf_wchi), chi2_e(ndf_wchi), chi3_e(ndf_wchi)
    real(kind=r_def)    :: jac(3,3,nqp_h,nqp_v), dj(nqp_h,nqp_v)
    real(kind=r_def)    :: basis_trial(ndf_wvec,nqp_h,nqp_v)

    ipanel = INT(panel_id(map_pid(1)))

    ! Make scalar wvec basis functions
    basis_trial(:,:,:) = 0.0_r_def
    do i = 1, 3
      ! Basis functions are linear and shouldn't be negative, but in case there
      ! is a direction, take the absolute value
      basis_trial(:,:,:) = basis_trial(:,:,:) + basis_wvec(i,:,:,:)
    end do

    do k = 1, nlayers

      ! Fill out coordinates and compute Jacobian for cell at face quad points
      chi1_e(:) = 0.0_r_def
      chi2_e(:) = 0.0_r_def
      chi3_e(:) = 0.0_r_def
      do df = 1, ndf_wchi
        chi1_e(df) = chi1(map_wchi(df)+k-1)
        chi2_e(df) = chi2(map_wchi(df)+k-1)
        chi3_e(df) = chi3(map_wchi(df)+k-1)
      end do

      call coordinate_jacobian(                                                &
              ndf_wchi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,                  &
              ipanel, basis_wchi, diff_basis_wchi, jac, dj                     &
      )

      ! Compute operator index
      ik = k + (col_idx-1)*nlayers

      ! Set matrix to zero
      projection_lhs(ik, :, :) = 0.0_r_def

      do df1 = 1, ndf_wvec
        do df2 = 1, ndf_wvec
          do qph = 1, nqp_h
            do qpv = 1, nqp_v
              projection_lhs(ik, df1, df2) = projection_lhs(ik, df1, df2)      &
                  + wqp_h(qph) * wqp_v(qpv) * basis_trial(df1,qph,qpv)         &
                  * basis_trial(df2,qph,qpv) * dj(qph,qpv)
            end do
          end do
        end do
      end do
    end do ! k

  end subroutine mass_matrix_scalarvec_code

end module sci_mass_matrix_scalarvec_kernel_mod
