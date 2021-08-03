!> Provides the ocean stochastic equation of state
module MOM_stoch_eos
! This file is part of MOM6. See LICENSE.md for the license.
use MOM_grid,            only : ocean_grid_type
use MOM_hor_index,       only : hor_index_type
use MOM_file_parser,     only : get_param, param_file_type
use MOM_random,          only : PRNG,random_2d_constructor,random_2d_norm
use MOM_time_manager,    only : time_type
use MOM_io,              only : vardesc, var_desc
use MOM_restart,         only : MOM_restart_CS,is_new_run
use MOM_diag_mediator,   only : register_diag_field,post_data,diag_ctrl,safe_alloc_ptr
use MOM_variables,       only : thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type
use MOM_restart,         only : register_restart_field
use MOM_isopycnal_slopes,only : vert_fill_TS
!use random_numbers_mod, only : getRandomNumbers,initializeRandomNumberStream,randomNumberStream

implicit none
#include <MOM_memory.h>

public MOM_stoch_eos_init
public MOM_stoch_eos_run
public MOM_calc_varT

real,private ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: l2_inv
real,private ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: rgauss
real, parameter,private :: tfac=0.27
real, parameter,private :: amplitude=0.624499 !0.39 is variance
integer        ,private :: seed
type(PRNG)  ::  rn_CS
!integer(8) :: npts
!real :: pi
logical first_time
  
type, public :: MOM_stoch_eos_CS 
  real,public  ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: pattern
                    !< Random pattern for stochastic EOS
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: phi
                    !< temporal correlation stochastic EOS (deugging)
  logical :: use_stoch_eos  !< If true, use the stochastic equation of state (Stanley et al. 2020)
  real :: stanley_coeff !< Coefficient correlating the temperature gradient 
                        !and SGS T variance; if <0, turn off scheme in all codes
  real :: stanley_a !a in exp(aX) in stochastic coefficient 
  real :: kappa_smooth    !< A diffusivity for smoothing T/S in vanished layers [Z2 T-1 ~> m2 s-1]
  integer :: id_stoch_eos  = -1, id_stoch_phi  = -1, id_tvar_sgs = -1
end type MOM_stoch_eos_CS


contains
  subroutine MOM_stoch_eos_init(G,Time,param_file,stoch_eos_CS,restart_CS,diag)
! initialization subroutine called by MOM.F90,
  type(param_file_type), intent(in)    :: param_file  !< structure indicating parameter file to parse
  type(ocean_grid_type), intent(in)    :: G
  type(time_type),       intent(in)    :: Time
  type(MOM_stoch_eos_CS), intent(inout) :: stoch_eos_CS
  type(MOM_restart_CS),  pointer       :: restart_CS
  type(diag_ctrl),       target, intent(inout) :: diag       !< to control diagnostics
  integer :: i,j
  type(vardesc)      :: vd
  seed=0
  ! contants
  !pi=2*acos(0.0)
  call get_param(param_file, "MOM_stoch_eos", "STOCH_EOS", stoch_eos_CS%use_stoch_eos, &
                 "If true, stochastic perturbations are applied "//&
                 "to the EOS in the PGF.", default=.false.)
  call get_param(param_file, "MOM_stoch_eos", "STANLEY_COEFF", stoch_eos_CS%stanley_coeff, &
                 "Coefficient correlating the temperature gradient "//&
                 "and SGS T variance.", default=-1.0)
  call get_param(param_file, "MOM_stoch_eos", "STANLEY_A", stoch_eos_CS%stanley_a, &
                 "Coefficient a which scales chi in stochastic perturbation of the "//&
                 "SGS T variance.", default=1.0)
  call get_param(param_file, "MOM_stoch_eos", "KD_SMOOTH", stoch_eos_CS%kappa_smooth, &
                 "A diapycnal diffusivity that is used to interpolate "//&
                 "more sensible values of T & S into thin layers.", &
                 units="m2 s-1", default=1.0e-6)

  !don't run anything if STANLEY_COEFF < 0
  if (stoch_eos_CS%stanley_coeff >= 0.0) then
  
    ALLOC_(stoch_eos_CS%pattern(G%isd:G%ied,G%jsd:G%jed)) ; stoch_eos_CS%pattern(:,:) = 0.0
    vd = var_desc("stoch_eos_pattern","nondim","Random pattern for stoch EOS",'h','1')
    call register_restart_field(stoch_eos_CS%pattern, vd, .false., restart_CS)
    ALLOC_(stoch_eos_CS%phi(G%isd:G%ied,G%jsd:G%jed)) ; stoch_eos_CS%phi(:,:) = 0.0
    ALLOC_(l2_inv(G%isd:G%ied,G%jsd:G%jed)) 
    ALLOC_(rgauss(G%isd:G%ied,G%jsd:G%jed)) 
    call get_param(param_file, "MOM_stoch_eos", "SEED_STOCH_EOS", seed, &
                 "Specfied seed for random number sequence ", default=0)
    call random_2d_constructor(rn_CS, G%HI, Time, seed)
    call random_2d_norm(rn_CS, G%HI, rgauss)
    ! fill array with approximation of grid area needed for decorrelation
    ! time-scale calculation
    do j=G%jsc,G%jec
       do i=G%isc,G%iec
          l2_inv(i,j)=1.0/(G%dxT(i,j)**2+G%dyT(i,j)**2)
       enddo
    enddo 
    if (is_new_run(restart_CS)) then
       do j=G%jsc,G%jec
          do i=G%isc,G%iec
             stoch_eos_CS%pattern(i,j)=amplitude*rgauss(i,j)
          enddo
       enddo
    endif

    !register diagnostics
    stoch_eos_CS%id_tvar_sgs = register_diag_field('ocean_model', 'tvar_sgs', diag%axesTL, Time, &
      'Parameterized SGS Temperature Variance ', 'None')
    if (stoch_eos_CS%use_stoch_eos) then
      stoch_eos_CS%id_stoch_eos = register_diag_field('ocean_model', 'stoch_eos', diag%axesT1, Time, &
        'random pattern for EOS', 'None')
      stoch_eos_CS%id_stoch_phi = register_diag_field('ocean_model', 'stoch_phi', diag%axesT1, Time, &
        'phi for EOS', 'None')
    endif
  endif
  
  end subroutine MOM_stoch_eos_init

  subroutine MOM_stoch_eos_run(G,u,v,delt,Time,stoch_eos_CS,diag)
  type(ocean_grid_type), intent(in)    :: G
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real,                  intent(in)    :: delt 
  type(time_type),       intent(in)    :: Time
  type(MOM_stoch_eos_CS), intent(inout) :: stoch_eos_CS
  type(diag_ctrl),       target, intent(inout) :: diag       !< to control diagnostics
! locals
  integer                                ::  i,j
  integer :: yr,mo,dy,hr,mn,sc
  real                                   :: phi,ubar,vbar

  call random_2d_constructor(rn_CS, G%HI, Time, seed)
  call random_2d_norm(rn_CS, G%HI, rgauss)
  ! advance AR(1)
  do j=G%jsc,G%jec
     do i=G%isc,G%iec
        ubar=0.5*(u(I,j,1)*G%mask2dCu(I,j)+u(I-1,j,1)*G%mask2dCu(I-1,j))
        vbar=0.5*(v(i,J,1)*G%mask2dCv(i,J)+v(i,J-1,1)*G%mask2dCv(i,J-1))
        phi=exp(-1*delt*tfac*sqrt((ubar**2+vbar**2)*l2_inv(i,j)))
        stoch_eos_CS%pattern(i,j)=phi*stoch_eos_CS%pattern(i,j) + amplitude*sqrt(1-phi**2)*rgauss(i,j)
        stoch_eos_CS%phi(i,j)=phi
     enddo
  enddo
  
  !print*,'stoch_run ubar',minval(u),maxval(u),minval(v),maxval(v)
  !print*,'stoch_run rp',minval(random_pattern),maxval(random_pattern),minval(rgauss),maxval(rgauss)
  !print*,'stoch_run phi',minval(phi_out),maxval(phi_out),minval(rgauss),maxval(rgauss)
  !print*,'indicies',G%isd,G%ied,G%idg_offset,G%jsd,G%jed,G%jdg_offset
  !print*,'rgauss bounds',lbound(rgauss,1),lbound(rgauss,2),ubound(rgauss,1),ubound(rgauss,2)
  !print*,'rp bounds',lbound(random_pattern,1),lbound(random_pattern,2),ubound(random_pattern,1),ubound(random_pattern,2)
  !print*,'subset',minval(rgauss(G%idg_offset:G%idg_offset+G%isd,G%jdg_offset:G%jdg_offset+G%jsd)),&
  !                maxval(rgauss(G%idg_offset:G%idg_offset+G%isd,G%jdg_offset:G%jdg_offset+G%jsd))
  !print*,'PJP should be posting',stoch_eos_CS%id_stoch_eos,stoch_eos_CS%id_stoch_phi
  !print*,'PJP min/max',minval(stoch_eos_CS%phi),maxval(stoch_eos_CS%phi),minval(stoch_eos_CS%pattern),maxval(stoch_eos_CS%pattern)
  !if (stoch_eos_CS%id_stoch_eos > 0) call post_data(stoch_eos_CS%id_stoch_eos, stoch_eos_CS%pattern, diag)!, mask=G%mask2dT)
  !if (stoch_eos_CS%id_stoch_phi > 0) call post_data(stoch_eos_CS%id_stoch_phi, stoch_eos_CS%phi, diag)!, mask=G%mask2dT)
  
  end subroutine MOM_stoch_eos_run


  subroutine MOM_calc_varT(G,GV,h,tv,stoch_eos_CS,dt)
  type(ocean_grid_type), intent(in)    :: G
  type(verticalGrid_type),                   intent(in)  :: GV  !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)  :: h   !< Layer thickness [H ~> m]
  type(thermo_var_ptrs), intent(inout) :: tv   !< Thermodynamics structure
  type(MOM_stoch_eos_CS), intent(inout) :: stoch_eos_CS
  real, intent(in)  :: dt    !< Time increment [T ~> s]
! locals
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: &
    T, &          ! The temperature (or density) [degC], with the values in
                  ! in massless layers filled vertically by diffusion.
    S             ! The filled salinity [ppt], with the values in
                  ! in massless layers filled vertically by diffusion.
  integer                                ::  i,j,k
  real :: hl(5)              ! Copy of local stencil of H [H ~> m]
  real :: dTdi2, dTdj2       ! Differences in T variance [degC2]

  ! This block does a thickness weighted variance calculation and helps control for
  ! extreme gradients along layers which are vanished against topography. It is
  ! still a poor approximation in the interior when coordinates are strongly tilted.
  if (.not. associated(tv%varT)) call safe_alloc_ptr(tv%varT, G%isd, G%ied, G%jsd, G%jed, GV%ke)
  
  call vert_fill_TS(h, tv%T, tv%S, stoch_eos_CS%kappa_smooth*dt, T, S, G, GV, halo_here=1, larger_h_denom=.true.)

  do k=1,G%ke
     do j=G%jsc,G%jec
        do i=G%isc,G%iec
          hl(1) = h(i,j,k) * G%mask2dT(i,j)
          hl(2) = h(i-1,j,k) * G%mask2dCu(I-1,j)
          hl(3) = h(i+1,j,k) * G%mask2dCu(I,j)
          hl(4) = h(i,j-1,k) * G%mask2dCv(i,J-1)
          hl(5) = h(i,j+1,k) * G%mask2dCv(i,J)

          ! SGS variance in i-direction [degC2]
          dTdi2 = ( ( G%mask2dCu(I  ,j) * G%IdxCu(I  ,j) * ( T(i+1,j,k) - T(i,j,k) ) &
                + G%mask2dCu(I-1,j) * G%IdxCu(I-1,j) * ( T(i,j,k) - T(i-1,j,k) ) &
                ) * G%dxT(i,j) * 0.5 )**2
          ! SGS variance in j-direction [degC2]
          dTdj2 = ( ( G%mask2dCv(i,J  ) * G%IdyCv(i,J  ) * ( T(i,j+1,k) - T(i,j,k) ) &
                + G%mask2dCv(i,J-1) * G%IdyCv(i,J-1) * ( T(i,j,k) - T(i,j-1,k) ) &
                ) * G%dyT(i,j) * 0.5 )**2
          tv%varT(i,j,k) = stoch_eos_CS%stanley_coeff * ( dTdi2 + dTdj2 )
          ! Turn off scheme near land
          tv%varT(i,j,k) = tv%varT(i,j,k) * (minval(hl) / (maxval(hl) + GV%H_subroundoff))
        enddo
     enddo
  enddo
  ! if stochastic, perturb
  if (stoch_eos_CS%use_stoch_eos) then
     do k=1,G%ke
        do j=G%jsc,G%jec
           do i=G%isc,G%iec
               tv%varT(i,j,k) = exp (stoch_eos_CS%stanley_a * stoch_eos_CS%pattern(i,j)) * tv%varT(i,j,k)
           enddo
        enddo
     enddo
  endif
  end subroutine MOM_calc_varT

end module MOM_stoch_eos

