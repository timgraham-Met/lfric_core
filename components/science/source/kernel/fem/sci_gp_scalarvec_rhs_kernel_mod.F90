!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Computes the right-hand side for a Galerkin projection of a scalar
!!        field into the scalar-form of a vector-valued function space.
!> @details A Galerkin projection solves the equation
!!          gamma_i \sum_j \int( phi_i * phi_j  dx) * c_j
!!          = gamma_i int( phi_i * f  dx)
!!          for the coefficients c_i of the projected field in the target space,
!!          where gamma are the test functions in the target space, phi_i are
!!          the basis functions in the target space, and f is the field to be
!!          projected.
!!          This kernel computes the right-hand side integral, for the case that
!!          phi_i = \sum_j |psi_j|, where psi_j are the vector basis functions
!!          in the target space. This allows a scalar field to be projected into
!!          a scalar form of a vector-valued function space.

module sci_gp_scalarvec_rhs_kernel_mod

  use argument_mod,             only : arg_type, func_type,                    &
                                       GH_FIELD, GH_REAL, GH_READ, GH_INC,     &
                                       GH_BASIS, GH_DIFF_BASIS,                &
                                       GH_QUADRATURE_XYoZ, CELL_COLUMN,        &
                                       ANY_SPACE_2, ANY_SPACE_3,               &
                                       ANY_DISCONTINUOUS_SPACE_3
  use constants_mod,            only : r_def, i_def
  use fs_continuity_mod,        only : Wchi
  use kernel_mod,               only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: gp_scalarvec_rhs_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                        &
        arg_type(GH_FIELD,   GH_REAL, GH_INC,  ANY_SPACE_2),                   &
        arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_SPACE_3),                   &
        arg_type(GH_FIELD*3, GH_REAL, GH_READ, Wchi),                          &
        arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3)      &
    /)
    type(func_type) :: meta_funcs(3) = (/                                      &
        func_type(ANY_SPACE_2, GH_BASIS),                                      &
        func_type(ANY_SPACE_3, GH_BASIS),                                      &
        func_type(Wchi, GH_BASIS, GH_DIFF_BASIS)                               &
    /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: gp_scalarvec_rhs_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: gp_scalarvec_rhs_code

contains

  !> @brief Computes the RHS for projecting orography into wvec function space
  !> @param[in]     nlayers               Number of layers in mesh
  !> @param[in,out] projection_rhs        RHS vector in wvec to calculate
  !> @param[in]     surface_altitude_w3   Surface altitude field in W3
  !> @param[in]     chi1                  First coordinate field
  !> @param[in]     chi2                  Second coordinate field
  !> @param[in]     chi3                  Third coordinate field
  !> @param[in]     panel_id              Field containing the mesh panel ID
  !> @param[in]     ndf_wvec              Number of DOFs per cell for vector
  !!                                      function space
  !> @param[in]     undf_wvec             Num of DOFs in this partition for
  !!                                      vector function space
  !> @param[in]     map_wvec              Dofmap for vector function space
  !> @param[in]     basis_wvec            Basis functions for vector function
  !!                                      space evaluated at quadrature points
  !> @param[in]     ndf_w3                Number of DOFs per cell for W3
  !> @param[in]     undf_w3               Num of DOFs in this partition for W3
  !> @param[in]     map_w3                Dofmap for W3
  !> @param[in]     basis_w3              Basis functions for W3 evaluated at
  !!                                      quadrature points
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
  subroutine gp_scalarvec_rhs_code(nlayers,                                    &
                                   projection_rhs,                             &
                                   surface_altitude_w3,                        &
                                   chi1, chi2, chi3,                           &
                                   panel_id,                                   &
                                   ndf_wvec, undf_wvec, map_wvec,              &
                                   basis_wvec,                                 &
                                   ndf_w3, undf_w3, map_w3,                    &
                                   basis_w3,                                   &
                                   ndf_wchi, undf_wchi, map_wchi,              &
                                   basis_wchi, diff_basis_wchi,                &
                                   ndf_pid, undf_pid, map_pid,                 &
                                   nqp_h, nqp_v, wqp_h, wqp_v)

    use sci_coordinate_jacobian_mod, only : coordinate_jacobian

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: nqp_h, nqp_v
    integer(kind=i_def), intent(in) :: ndf_wvec, ndf_w3, ndf_wchi, ndf_pid
    integer(kind=i_def), intent(in) :: undf_wvec, undf_w3, undf_wchi, undf_pid

    integer(kind=i_def), intent(in) :: map_wvec(ndf_wvec)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_wchi(ndf_wchi)
    integer(kind=i_def), intent(in) :: map_pid(ndf_pid)

    real(kind=r_def),    intent(in) :: basis_wvec(3,ndf_wvec,nqp_h,nqp_v)
    real(kind=r_def),    intent(in) :: basis_w3(1,ndf_w3,nqp_h,nqp_v)
    real(kind=r_def),    intent(in) :: basis_wchi(1,ndf_wchi,nqp_h,nqp_v)
    real(kind=r_def),    intent(in) :: diff_basis_wchi(3,ndf_wchi,nqp_h,nqp_v)
    real(kind=r_def),    intent(in) :: wqp_h(nqp_h)
    real(kind=r_def),    intent(in) :: wqp_v(nqp_v)

    real(kind=r_def),    intent(inout) :: projection_rhs(undf_wvec)
    real(kind=r_def),    intent(in)    :: surface_altitude_w3(undf_w3)
    real(kind=r_def),    intent(in)    :: chi1(undf_wchi)
    real(kind=r_def),    intent(in)    :: chi2(undf_wchi)
    real(kind=r_def),    intent(in)    :: chi3(undf_wchi)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)

    ! Internal variables
    integer(kind=i_def) :: df, df_w3, qph, qpv, ipanel, i, k
    real(kind=r_def)    :: chi1_e(ndf_wchi), chi2_e(ndf_wchi), chi3_e(ndf_wchi)
    real(kind=r_def)    :: jac(3,3,nqp_h,nqp_v), dj(nqp_h,nqp_v)
    real(kind=r_def)    :: basis_trial(ndf_wvec,nqp_h,nqp_v)

    ipanel = INT(panel_id(map_pid(1)))

    ! Make scalar wvec basis functions
    basis_trial(:,:,:) = 0.0_r_def
    do i = 1, 3
      ! Basis functions are linear and shouldn't be negative, but in case there
      ! is a direction, take the absolute value
      basis_trial(:,:,:) = basis_trial(:,:,:) + ABS(basis_wvec(i,:,:,:))
    end do

    do k = 0, nlayers-1

      ! Fill out coordinates and compute Jacobian for cell at face quad points
      chi1_e(:) = 0.0_r_def
      chi2_e(:) = 0.0_r_def
      chi3_e(:) = 0.0_r_def
      do df = 1, ndf_wchi
        chi1_e(df) = chi1(map_wchi(df)+k)
        chi2_e(df) = chi2(map_wchi(df)+k)
        chi3_e(df) = chi3(map_wchi(df)+k)
      end do

      call coordinate_jacobian(                                                &
              ndf_wchi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,                  &
              ipanel, basis_wchi, diff_basis_wchi, jac, dj                     &
      )

      do df = 1, ndf_wvec
        do df_w3 = 1, ndf_w3
          do qph = 1, nqp_h
            do qpv = 1, nqp_v
              projection_rhs(map_wvec(df)+k) = projection_rhs(map_wvec(df)+k)  &
                  + wqp_h(qph) * wqp_v(qpv)                                    &
                  * basis_trial(df,qph,qpv) * basis_w3(1,df_w3,qph,qpv)        &
                  * surface_altitude_w3(map_w3(df_w3)+k) * dj(qph,qpv)
            end do
          end do
        end do
      end do
    end do ! k


  end subroutine gp_scalarvec_rhs_code

end module sci_gp_scalarvec_rhs_kernel_mod
