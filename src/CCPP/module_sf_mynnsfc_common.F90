!>\file module_sf_mynnedmf_common.F90
!!Define Model-specific constants/parameters.
!!This module will be used at the initialization stage
!!where all model-specific constants are read and saved into
!!memory. This module is then used again in the MYNNSFC_* scheme. All
!!MYNN-specific constants are declared globally in the main
!module (module_sf_mynnsfc) further below:

!>\ingroup gp_mynnsfc
!! This module defines model-specific constants/parameters.
!====================================================================

module module_sf_mynnsfc_common
  use machine,  only : kind_phys

  implicit none
  save

  ! To be specified from dycore
  real(kind_phys):: cp           != 7.*r_d/2. (J/kg/K)
  real(kind_phys):: cpv          != 4.*r_v    (J/kg/K) Spec heat H2O gas
  real(kind_phys):: cice         != 2106.     (J/kg/K) Spec heat H2O ice
  real(kind_phys):: cliq         != 4190.     (J/kg/K) Spec heat H2O liq
  real(kind_phys):: p608         != R_v/R_d-1.
  real(kind_phys):: ep_2         != R_d/R_v
  real(kind_phys):: grav         != accel due to gravity
  real(kind_phys):: karman       != von Karman constant
  real(kind_phys):: t0c          != svpt0 temperature of water at freezing, 273.15 K
  real(kind_phys):: rcp          != r_d/cp
  real(kind_phys):: r_d          != 287.  (J/kg/K) gas const dry air
  real(kind_phys):: r_v          != 461.6 (J/kg/K) gas const water
  real(kind_phys):: xlf          != 0.35E6 (J/kg) fusion at 0 C
  real(kind_phys):: xlv          != 2.50E6 (J/kg) vaporization at 0 C
  real(kind_phys):: xls          != 2.85E6 (J/kg) sublimation
  real(kind_phys):: rvovrd       != r_v/r_d != 1.608

  ! Specified locally
  real(kind_phys),parameter:: tref  = 300.0   ! reference temperature (K)
  ! real(kind_phys),parameter:: p1000mb=100000.0
  ! real(kind_phys),parameter:: svp1  = 0.6112 !(kPa)
  ! real(kind_phys),parameter:: svp2  = 17.67  !(dimensionless)
  ! real(kind_phys),parameter:: svp3  = 29.65  !(K)
  real(kind_phys),parameter:: wfa_max = 800e6   !kg-1
  real(kind_phys),parameter:: wfa_min = 1e6     !kg-1
  real(kind_phys),parameter:: ifa_max = 300e6   !kg-1
  real(kind_phys),parameter:: ifa_min = 0.0     !kg-1
  real(kind_phys),parameter:: wfa_ht  = 2000.   !meters
  real(kind_phys),parameter:: ifa_ht  = 10000.  !meters
   
  ! To be derived in the init routine
  real(kind_phys):: ep_3    !<= 1.-ep_2 != 0.378
  real(kind_phys):: gtr    !<= grav/tref
  real(kind_phys):: rk     !<= cp/r_d
  real(kind_phys):: tv0    !<=  p608*tref
  real(kind_phys):: tv1    !<= (1.+p608)*tref
  real(kind_phys):: xlscp  !<= (xlv+xlf)/cp
  real(kind_phys):: xlvcp  !<= xlv/cp
  real(kind_phys):: g_inv  !<= 1./grav
  real(kind_phys):: ep_1    !Rv/Rd-1

end module module_sf_mynnsfc_common
