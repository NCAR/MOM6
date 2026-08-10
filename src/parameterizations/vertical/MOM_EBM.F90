! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains the MOM6 version of the estuary box model parameterization,
!! which is based on the algorithm described in Sun, Q., Whitney, M. M., Bryan, F. O.,
!! & Tseng, Y. (2017). A box model for representing estuarine physical processes in
!! Earth system models. Ocean Modelling, 112, 139-153.
!! https://doi.org/10.1016/j.ocemod.2017.03.004
module MOM_EBM

use MOM_error_handler, only : MOM_error, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type

implicit none ; private

public EBM_init, calculate_EBM, EBM_is_used

!> Control structure including parameters for the estuary box model.
type, public :: EBM_cs ; private

  real :: tide_amp !< Averaged tidal amplitude at estuary mouth [m]
  real :: W_h      !< Estuary head width [m]
  real :: H        !< Estuary averaged depth [m]
  real :: a2       !< A constant of tidal diffusion [nondim]
  real :: a1       !< A constant of estuarine mixing length [nondim]
  real :: h0       !< A constant of ratio of geometry: h_l/H [nondim]
  real :: g        !< Gravitational acceleration [m s-2]
  real :: rho_ref    !< Reference density for linear equation of state [kg m-3]
  real :: beta_S   !< Saline contraction coefficient [ppt-1]
  real :: Sc       !< Schmidt number [nondim]

end type EBM_cs

character(len=40) :: mdl = "MOM_EBM"  !< This module's name.

contains

!> Initializes the estuary box model parameterization.
!! Returns .true. if the parameterization is enabled.
logical function EBM_init(param_file, CS)

  type(param_file_type), intent(in)    :: param_file !< Run-time parameter file handle
  type(EBM_cs),          intent(inout) :: CS         !< EBM control structure

  ! This include declares and sets the variable "version".
# include "version_variable.h"

  call get_param(param_file, mdl, "USE_EBM", EBM_init, default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, &
       "Estuary box model parameterization", all_default=.not.EBM_init)
  call get_param(param_file, mdl, "USE_EBM", EBM_init, &
                 "If true, enables the estuary box model (EBM) parameterization. ", &
                 default=.false.)

  if (.not. EBM_init) return

  call get_param(param_file, mdl, "EBM_TIDAL_AMPLITUDE", CS%tide_amp, &
                 "Averaged tidal amplitude at the estuary mouth.", &
                 units="m", default=1.0)

  call get_param(param_file, mdl, "EBM_HEAD_WIDTH", CS%W_h, &
                 "Estuary head width.", &
                 units="m", default=2000.0)

  call get_param(param_file, mdl, "EBM_DEPTH", CS%H, &
                 "Estuary averaged depth.", &
                 units="m", default=10.0)

  call get_param(param_file, mdl, "EBM_A1", CS%a1, &
                 "A constant of estuarine mixing length.", &
                 units="nondim", default=0.876)

  call get_param(param_file, mdl, "EBM_A2", CS%a2, &
                 "A constant of tidal diffusion.", &
                 units="nondim", default=0.0)

  call get_param(param_file, mdl, "EBM_H0", CS%h0, &
                 "A constant of ratio of geometry: h_l/H.", &
                 units="nondim", default=0.5)

  call get_param(param_file, mdl, "EBM_G", CS%g, &
                 "Gravitational acceleration used in the EBM.", &
                 units="m s-2", default=9.8)

  call get_param(param_file, mdl, "EBM_RHO_REF", CS%rho_ref, &
                 "Reference density for the linear EoS used in the EBM.", &
                 units="kg m-3", default=1000.0)

  call get_param(param_file, mdl, "EBM_BETA_S", CS%beta_S, &
                 "Saline contraction coefficient used in the EBM.", &
                 units="ppt-1", default=7.7e-4)

  call get_param(param_file, mdl, "EBM_SC", CS%Sc, &
                 "Schmidt number used in the EBM.", &
                 units="nondim", default=2.2)

end function EBM_init

!> Calculates estuary box model exchange and distributes river runoff over the
!! first 4 layers, updating layer thickness, temperature, salinity, and net mass flux.
!! The estuary exchange fluxes (Q_u, Q_l, S_u) are computed by estuary_box_model
!! using the EBM parameters in CS together with the provided river discharge and
!! lower-layer salinity.
subroutine calculate_EBM(CS, lrunoff, S_l, EnthalpyConst, netMassIn, T2d_col, S_col, h2d_col)

  type(EBM_cs), intent(in)    :: CS            !< EBM control structure
  real,         intent(in)    :: lrunoff       !< River runoff for this column [H ~> m or kg m-2]
  real,         intent(in)    :: S_l           !< Salinity at the estuary lower layer [S ~> ppt]
  real,         intent(in)    :: EnthalpyConst !< Enthalpy constant [nondim]
  real,         intent(inout) :: netMassIn     !< Net mass entering this column [H ~> m or kg m-2]
  real,         intent(inout) :: T2d_col(:)    !< Temperature in first 4 layers [C ~> degC]
  real,         intent(inout) :: S_col(:)      !< Salinity in first 4 layers [S ~> ppt]
  real,         intent(inout) :: h2d_col(:)    !< Layer thickness in first 4 layers [H ~> m or kg m-2]

  ! local variables
  real :: Q_u        ! EBM upper layer volume flux [m3 s-1]
  real :: Q_l        ! EBM lower layer volume flux [m3 s-1]
  real :: S_u        ! EBM upper layer salinity [ppt]
  real :: dThickness ! Change in layer thickness [H ~> m or kg m-2]
  real :: dTemp      ! Integrated change in layer temperature [C H ~> degC m or degC kg m-2]
  real :: dSalt      ! Integrated change in layer salinity [S H ~> ppt m or ppt kg m-2]
  real :: Temp_in    ! Temperature of the incoming mass flux [C ~> degC]
  real :: Salin_in   ! Salinity of the incoming mass flux [S ~> ppt]
  real :: hOld       ! Original layer thickness before update [H ~> m or kg m-2]
  real :: Ithickness ! Inverse of the updated layer thickness [H-1 ~> m-1 or m2 kg-1]
  integer :: k       ! Layer index [nondim]

  ! GMM, TODO: need to specify Hu, instead of hard code thickness...

  ! Distribute river runoff over the first 4 layers
  do k = 1, size(T2d_col)
    dThickness = lrunoff * 0.25  ! Each of the 4 layers receives 1/4 of the runoff
    dTemp = 0.
    dSalt = 0.

    netMassIn = netMassIn - dThickness
    Temp_in  = T2d_col(k)
    Salin_in = 0.0
    dTemp = dTemp + dThickness * Temp_in * EnthalpyConst

    hOld = h2d_col(k)
    h2d_col(k) = h2d_col(k) + dThickness
    if (h2d_col(k) > 0.0) then
      Ithickness = 1.0 / h2d_col(k)
      if (dThickness /= 0. .or. dTemp /= 0.) T2d_col(k) = (hOld * T2d_col(k) + dTemp) * Ithickness
      if (dThickness /= 0. .or. dSalt /= 0.) S_col(k)   = (hOld * S_col(k)   + dSalt) * Ithickness
    end if
  end do

  ! GMM, todo
  call estuary_box_model(CS, lrunoff, S_l, Q_u, Q_l, S_u)

end subroutine calculate_EBM

!> Reads the parameter "USE_EBM" and returns state.
!! This function allows other modules to know whether this parameterization will
!! be used without needing to duplicate the log entry.
logical function EBM_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  call get_param(param_file, mdl, "USE_EBM", EBM_is_used, &
                 default=.false., do_not_log=.true.)

end function EBM_is_used

!> Calculate the estuary box model (EBM) exchange.
!!
!! The EBM is assumed steady state, flat bottom and flat surface with a straight
!! channel and rectangular cross-section. It is built on three global conservation
!! laws: water mass, water volume, and potential energy conservation.
!!
!! Estuary Box geometry:
!!              ^ z
!!              |____________________         _______
!!              |                    |               |
!!  Q_u,S_u <--+--  upper layer   <-+-- Q_r         |
!!              |--------------------|        -+-    | H
!!     Q_l,S_l -+->  lower layer     |         | h_l |
!!           ___|____________________|        _|_____|
!!          x   0                   -LE
!!
!! Note: the negative lower layer volume flux Q_l leaves the ocean; the positive
!! upper layer volume flux Q_u flows into the ocean.
!!
!! Reference: Sun, Q., Whitney, M. M., Bryan, F. O., & Tseng, Y. (2017).
!! A box model for representing estuarine physical processes in
!! Earth system models. Ocean Modelling, 112, 139-153.
!! https://doi.org/10.1016/j.ocemod.2017.03.004
subroutine estuary_box_model(CS, Q_r, S_l, Q_u, Q_l, S_u)

  type(EBM_cs), intent(in)  :: CS   !< EBM control structure
  real,         intent(in)  :: Q_r  !< River discharge [m3 s-1]
  real,         intent(in)  :: S_l  !< Salinity at estuary lower layer [ppt]
  real,         intent(out) :: Q_u  !< Upper layer volume flux [m3 s-1]
  real,         intent(out) :: Q_l  !< Lower layer volume flux [m3 s-1]
  real,         intent(out) :: S_u  !< Salinity at estuary upper layer [ppt]

  ! local variables
  real :: ERR_EBM     ! Closure error of the EBM potential energy budget [kg m s-3]
  real :: rho_r       ! River water density [kg m-3]
  real :: rho_l       ! Lower layer inflow density [kg m-3]
  real :: rho_u       ! Upper layer outflow density [kg m-3]
  real :: u_t         ! Tidal current amplitude near bottom [m s-1]
  real :: u_r         ! Riverine velocity at head of EBM [m s-1]
  real :: u_l         ! Estuarine inflow velocity at mouth of EBM [m s-1]
  real :: u_u         ! Estuarine outflow velocity at mouth of EBM [m s-1]
  real :: u_bar       ! Layer-averaged net velocity in EBM [m s-1]
  real :: c_wave      ! Densimetric wave phase speed [m s-1]
  real :: ur0         ! Densimetric riverine Froude number [nondim]
  real :: ut0         ! Densimetric tidal current Froude number [nondim]
  real :: ul0         ! Dimensionless lower layer inflow Froude number [nondim]
  real :: uu0         ! Dimensionless upper layer outflow Froude number [nondim]
  real :: R0          ! Layer densimetric riverine Froude number [nondim]
  real :: T0          ! Layer densimetric tidal Froude number [nondim]
  real :: h_l         ! Lower layer water depth of EBM [m]
  real :: a, b, c, d  ! Normalized coefficients of the cubic equation for ul0 [nondim]
  real :: AD          ! Advective potential energy flux term for EBM closure check [kg m s-3]
  real :: HD          ! Horizontal diffusive potential energy flux term [kg m s-3]
  real :: VD          ! Vertical diffusive potential energy flux term [kg m s-3]
  real :: LF          ! Lateral friction potential energy flux term [kg m s-3]
  real, dimension(3,2) :: roots   ! Roots of the 3rd-order polynomial [nondim];
                                  !! roots(:,1) = real part, roots(:,2) = imaginary part
  integer, dimension(3) :: mask   ! Selection mask for physically valid roots [nondim]
  integer :: i, n                 ! Loop and counter indices [nondim]

  real, parameter :: PI = 4.0 * atan(1.0)  !< Ratio of circumference to diameter [nondim]

  ! Skip computations if S_l <= 0, using limit as S_l->0
  if (S_l <= 0.0) then
    Q_u = Q_r
    Q_l = 0.0
    S_u = 0.0
    return
  end if

  ! Water densities
  rho_r = CS%rho_ref
  rho_l = CS%rho_ref * (1.0 + CS%beta_S * S_l)

  ! River and tidal velocities
  u_t    = -CS%tide_amp * sqrt(CS%g / CS%H)        ! Tidal velocity (toward river)
  u_r    = Q_r / (CS%W_h * CS%H * (1.0 - CS%h0))  ! Riverine velocity at head of upper layer
  c_wave = sqrt(CS%beta_S * S_l * CS%g * CS%H)     ! Densimetric wave phase speed

  ! Dimensionless parameters
  ur0 = u_r / c_wave
  ut0 = u_t / c_wave
  R0  = ur0 * (1.0 - CS%h0)
  T0  = ut0 * (1.0 - CS%h0) / PI

  ! Coefficients of the cubic equation for dimensionless lower layer inflow (ul0)
  a = -CS%h0**3.0

  b = 2.0 * CS%h0**2.0 * ((2.0 - CS%h0) * R0 - CS%a2 * T0)

  c = 0.096 * CS%a1 * CS%h0 * (CS%Sc**2.0 * R0)**(-1.0/3.0) * R0 &
    - CS%h0 * ((2.0 - CS%h0) * R0 * (R0 - 2.0 * CS%a2 * T0) + CS%a2**2.0 * T0**2.0)

  d = -0.048 * CS%a1 * (CS%Sc**2.0 * R0)**(-1.0/3.0) &
    * R0 * (R0 - 2.0 * CS%a2 * T0)

  call cubsolve(b/a, c/a, d/a, roots)

  ! Select the physically valid root: real, negative lower layer inflow
  mask = 0
  n    = 0
  do i = 1, 3
    if (roots(i,1) < 0.0 .and. roots(i,2) == 0.0) then
      mask(i) = 1
      n = n + 1
    end if
  end do

  if (n == 0) then
    call MOM_error(WARNING, "MOM_EBM estuary_box_model: no valid EBM solution found.")
    ul0 = 0.0
  else if (n == 1) then
    ul0 = sum(roots(1:3,1) * real(mask))
  else
    call MOM_error(WARNING, "MOM_EBM estuary_box_model: multiple valid EBM solutions found.")
    ul0 = 0.0
  end if

  ! Upper layer salinity and volume fluxes at EBM mouth
  uu0 = R0 / (1.0 - CS%h0) - CS%h0 / (1.0 - CS%h0) * ul0
  S_u = (-S_l * ul0 * CS%h0 - S_l * CS%a2 * T0) / (R0 - ul0 * CS%h0 - CS%a2 * T0)
  Q_l = ul0 * CS%h0 * CS%H * CS%W_h * c_wave
  Q_u = uu0 * (1.0 - CS%h0) * CS%H * CS%W_h * c_wave

  ! Verify closure of EBM potential energy budget
  u_l   = ul0 * c_wave
  u_u   = uu0 * c_wave
  u_bar = Q_r / (CS%W_h * CS%H)
  h_l   = CS%H * CS%h0
  rho_u = CS%rho_ref * (1.0 + CS%beta_S * S_u)

  AD = 0.5 * CS%g * rho_l * u_l * h_l**2.0 &
     + 0.5 * CS%g * (rho_u * u_u - rho_r * u_r) * (CS%H**2.0 - h_l**2.0)

  HD = -0.5 * CS%a2 * CS%g * (rho_l - rho_u) * (CS%H**2.0 - h_l**2.0) * u_t / PI

  VD = -0.5 * CS%g * (rho_u - rho_l) &
     * (rho_l + rho_u - 2.0 * rho_r) / (rho_l - rho_r) &
     * 0.024 * CS%a1 * CS%H**2.0 &
     * (c_wave**4.0 / (u_bar * CS%Sc**2.0))**(1.0/3.0)

  LF = -0.25 * CS%g * Q_l / CS%W_h &
     * ( (rho_u**2.0 + 2.0 * rho_l * rho_r - 2.0 * rho_u * rho_r) &
         * (CS%H - h_l) &
       - rho_r**2.0 * CS%H + rho_l**2.0 * h_l) / (rho_l - rho_r)

  ! GMM
  ! TODO: make this a 2D field and add option to save as a diagnostic?
  ERR_EBM = AD - HD - VD - LF

  ! Effective upper layer salinity for MOM6; Q_l is negative
  S_u = -Q_l * S_l / Q_u

end subroutine estuary_box_model

!> Solves the depressed cubic equation x^3 + ax^2 + bx + c = 0 analytically.
!! Returns all three roots; roots(:,1) are the real parts and roots(:,2) the
!! imaginary parts.
subroutine cubsolve(a, b, c, roots)

  real,                  intent(in)  :: a, b, c  !< Coefficients of x^3 + ax^2 + bx + c = 0 [nondim]
  real, dimension(3,2),  intent(out) :: roots     !< Roots [nondim]: column 1 = real, column 2 = imaginary

  ! local variables
  real :: Q     ! Intermediate cubic parameter [nondim]
  real :: R     ! Intermediate cubic parameter [nondim]
  real :: Rsqu  ! R squared [nondim]
  real :: Qcub  ! Q cubed [nondim]
  real :: SQ    ! Square root of Q [nondim]
  real :: theta ! Angle for the three-root case [nondim]
  real :: X     ! Intermediate root quantity [nondim]
  real :: Y     ! Intermediate root quantity [nondim]
  real :: XY    ! Sum X + Y [nondim]

  real, parameter :: PI = 4.0 * atan(1.0) !< Ratio of circumference to diameter [nondim]

  Q    = (a**2.0 - 3.0 * b) / 9.0
  R    = (2.0 * a**3.0 - 9.0 * a * b + 27.0 * c) / 54.0
  Rsqu = R**2.0
  Qcub = Q**3.0

  if (Rsqu < Qcub) then  ! Three distinct real roots
    theta      = acos(R / sqrt(Qcub))
    SQ         = sqrt(Q)
    roots(1,1) = -2.0 * SQ * cos(theta / 3.0) - a / 3.0
    roots(2,1) = -2.0 * SQ * cos((theta + 2.0 * PI) / 3.0) - a / 3.0
    roots(3,1) = -2.0 * SQ * cos((theta - 2.0 * PI) / 3.0) - a / 3.0
    roots(1,2) = 0.0
    roots(2,2) = 0.0
    roots(3,2) = 0.0
    return
  end if

  ! One real root and two conjugate complex roots
  X = -(abs(R) + sqrt(Rsqu - Qcub))**(1.0/3.0)
  if (R < 0.0) X = -X

  if (X == 0.0) then
    Y = 0.0
  else
    Y = Q / X
  end if

  XY         = X + Y
  roots(1,1) = XY - a / 3.0
  roots(1,2) = 0.0
  roots(2,1) = -0.5 * XY - a / 3.0
  roots(3,1) = -0.5 * XY - a / 3.0
  roots(2,2) =  sqrt(3.0) * (X - Y) / 2.0
  roots(3,2) = -sqrt(3.0) * (X - Y) / 2.0

end subroutine cubsolve

end module MOM_EBM
