!> \file module_sf_mynnsfc_driver.F90
!!  This serves as the interface between the WRF surface driver and the MYNN
!!  surface-layer scheme(s) in module_sfc_mynnsfc_land.F90, module_sfc_mynnsfc_ice.F90, and
!!  module_sfc_mynnsfc_water.F90.

!>\ingroup gsd_mynnsfc
!> The following references best describe the code within
!!  Olson, J. B., T. Smirnova, J. S. Kenyon, D. Turner, J. M. Brown, W. Zheng,
!!  and B. W. Green, 2021: A description of the MYNN Surface-Layer Scheme.
!!  NOAA Tech. Memo. OAR GSL-67, 26 pp., doi:10.25923/f6a8-bc75.
!!
!! Note that the old tech note (Olson et al. 2021, doi:10.25923/f6a8-bc75) is
!! still the best reference for the older version of the scheme that is used in
!! HRRRv4 and RAPv5. An updated tech note is underway... The physics in this
!! universalized and modularized version is still very similar to that in the
!! 2021 tech note, but different model-dependent configuration options have been
!! change to MYNN-specific options, such as:
!!
!!   old option "iz0tlnd" is now "sf_mynn_sfcflux_land":
!!   (default) =0: Zilitinkevich (1995); Czil now set to 0.095
!!             =1: Czil_new (modified according to Chen & Zhang 2008)
!!             =2: Modified Yang et al (2002, 2008) - generalized for all landuse
!!             =3: constant zt = z0/7.4 (original form; Garratt 1992)
!!
!!   old option "isftcflx" is now "sf_mynn_sfcflux_water":
!!             =0: z0, zt, and zq from COARE 3.0 (Fairall et al. 2003)
!!   (default) =1: z0, zt, and zq from COARE 3.5 (Edson et al 2013)
!!             =2: z0 from Davis et al (2008), zt & zq from COARE 3.5
!!             =3: z0 from Davis et al (2008), zt & zq from Garratt (1992)
!!             =4: z0 from Taylor and Yelland (2004), zt and zq from COARE 3.5
!!
!! This way, the new sf_mynn* options exist in all model frameworks (WRF/MPAS/CCPP)
!! and they carry the same meaning, wheras before, they carried a different meaning for
!! different schemes.
!=================================================================================================================
 module module_sf_mynnsfc_driver
!=================================================================================================================
 use module_sf_mynnsfc_common,only: kind_phys,cp,lsm_ruc => ruclsmscheme
 use module_sf_mynnsfc_land,  only: mynnsfc_land, psi_init
 use module_sf_mynnsfc_water, only: mynnsfc_water
 use module_sf_mynnsfc_ice,   only: mynnsfc_ice
 
 implicit none

 !Global configuration options, to be moved to namelist variables:
 integer, parameter :: psi_opt      = 0       !0: mynn, 1: gfs
 logical, parameter :: compute_flux = .true.  !whether or not to compute the surface fluxes
 logical, parameter :: compute_diag = .true.  !whether or not to compute the 2- and 10-m diagnostics
 
 private
 public:: mynnsfc_driver, &
          mynnsfc_init,   &
          mynnsfclay_finalize      

 contains
  
!=================================================================================================================
!>\section arg_table_sf_mynn_init
!!
 subroutine mynnsfc_init(allowed_to_read,errmsg,errflg)

    logical, intent(in)           :: allowed_to_read
    character(len=*), intent(out) :: errmsg
    integer, intent(out)          :: errflg

    !Fill the PSIM and PSIH tables. This code was leveraged from
    !module_sf_sfclayrev.F (from Pedro Jimenez and Jimy Dudhia).
    !This subroutine returns a blended form from Dyer and Hicks (1974)
    !and Grachev et al (2000) for unstable conditions and the form
    !from Cheng and Brutsaert (2005) for stable conditions.

     call psi_init(psi_opt,errmsg,errflg)

 end subroutine mynnsfc_init
!=================================================================================================================
!>\section arg_table_sf_mynn_finalize
!!
 subroutine mynnsfclay_finalize(errmsg,errflg)
!=================================================================================================================

!--- output arguments:                                                                                                                                                                              
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg

  errmsg = ' '
  errflg = 0

 end subroutine mynnsfclay_finalize  
!=================================================================================================================
 subroutine mynnsfc_driver(                                                &
        !3d variables
        u3d      , v3d      , t3d      , qv3d     , p3d      , dz8w      , &
        th3d     , rho3d    ,                                              &
        !GFS-related input
        sigmaf   , vegtype  , shdmax   , ivegsrc  ,                        &  !intent(in)
        z0pert   , ztpert   , redrag   , sfc_z0_type ,                     &  !intent(in)
        !2d variables
        psfcpa   , chs      , chs2     , cqs      , cqs2     , cpm       , &
        znt      , ust      , ustm     , pblh     , mavail   , zol       , &
        mol      , rmol     , regime   , psim     , psih     , xland     , &
        hfx      , qfx      , lh       , tsk      , flhc     , flqc      , &
        qgh      , qsfc     , u10      , v10      , th2      , t2        , &
        q2       , snowh    , gz1oz0   , wspd     , br       , dx        , &
        ch       , qcg      , ck       , cka      , cd       , cda       , &
        stress   , hflx     , qflx     , cm       , fm       , fh        , &
        fm10     , fh2      , tsurf    ,                                   &
        !configuration options
        spp_pbl  , pattern_spp_pbl     ,                                   &
        sf_mynn_sfcflux_water          ,                                   &
        sf_mynn_sfcflux_land           ,                                   &
        isfflx   , restart  , cycling  , initflag , flag_iter, flag_lsm  , &
        !model information
        itimestep,                                                         &
        ids      , ide      , jds      , jde      , kds      , kde       , &
        ims      , ime      , jms      , jme      , kms      , kme       , &
        its      , ite      , jts      , jte      , kts      , kte       , &
        errmsg   , errflg                                                  &
        )
!-------------------------------------------------------------------
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- t3d         3d temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- p3d         3d pressure (pa)
!-- rho3d       3d density (kg/m3)
!-- dz8w        3d dz between full levels (m)
!-- cp          heat capacity at constant pressure for dry air (j/kg/k)
!-- g           acceleration due to gravity (m/s^2)
!-- rovcp       r/cp
!-- r           gas constant for dry air (j/kg/k)
!-- xlv         latent heat of vaporization for water (j/kg)
!-- psfcpa      surface pressure (pa)
!-- znt         roughness length (m)
!-- ust         u* in similarity theory (m/s)
!-- ustm        u* in similarity theory (m/s) w* added to wspd. this is
!               used to couple with tke scheme but not in mynn.
!               (as of now, ustm = ust in this version)
!-- pblh        pbl height from previous time (m)
!-- mavail      surface moisture availability (between 0 and 1)
!-- zol         z/l height over monin-obukhov length
!-- mol         t* (similarity theory) (k)
!-- rmol        reciprocal of m-o length (/m)
!-- regime      flag indicating pbl regime (stable, unstable, etc.)
!-- psim        similarity stability function for momentum
!-- psih        similarity stability function for heat
!-- xland       land mask (1 for land, 2 for water)
!-- hfx         upward heat flux at the surface (w/m^2)
!-- qfx         upward moisture flux at the surface (kg/m^2/s)
!-- lh          net upward latent heat flux at surface (w/m^2)
!-- tsk         surface temperature (k)
!-- flhc        exchange coefficient for heat (w/m^2/k)
!-- flqc        exchange coefficient for moisture (kg/m^2/s)
!-- chs         heat/moisture exchange coefficient for lsm (m/s)
!-- qgh         lowest-level saturated mixing ratio
!-- qsfc        qv (specific humidity) at the surface
!-- qsfcmr      qv (mixing ratio) at the surface
!-- u10         diagnostic 10m u wind
!-- v10         diagnostic 10m v wind
!-- th2         diagnostic 2m theta (k)
!-- t2          diagnostic 2m temperature (k)
!-- q2          diagnostic 2m mixing ratio (kg/kg)
!-- snowh       snow height (m)
!-- gz1oz0      log((z1+znt)/znt) where znt is roughness length 
!-- wspd        wind speed at lowest model level (m/s)
!-- br          bulk richardson number in surface layer
!-- isfflx      isfflx=1 for surface heat and moisture fluxes
!-- dx          horizontal grid size (m)
!-- svp1        constant for saturation vapor pressure (=0.6112 kpa)
!-- svp2        constant for saturation vapor pressure (=17.67 dimensionless)
!-- svp3        constant for saturation vapor pressure (=29.65 k)
!-- svpt0       constant for saturation vapor pressure (=273.15 k)
!-- ep1         constant for virtual temperature (rv/rd - 1) (dimensionless)
!-- ep2         constant for spec. hum. calc (rd/rv = 0.622) (dimensionless)
!-- ep3         constant for spec. hum. calc (1 - rd/rv = 0.378 ) (dimensionless)
!-- karman      von karman constant
!-- ck          enthalpy exchange coeff at 10 meters
!-- cd          momentum exchange coeff at 10 meters
!-- cka         enthalpy exchange coeff at the lowest model level
!-- cda         momentum exchange coeff at the lowest model level
!-- sf_mynn_sfcflux_land chooses the bulk-flux algorithm used over land (described above)
!-- sf_mynn_sfcflux_water chooses the bulk-flux algorithm used over water (described above)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!=================================================================================================================

!--- input configuration arguments:
 integer,intent(in):: ids,ide,jds,jde,kds,kde, &
                      ims,ime,jms,jme,kms,kme, &
                      its,ite,jts,jte,kts,kte
 logical,intent(in),optional:: restart,cycling
 integer,intent(in):: itimestep
 integer,intent(in),optional:: initflag
 integer,intent(in):: isfflx
 integer,intent(in),optional:: sf_mynn_sfcflux_water
 integer,intent(in),optional:: sf_mynn_sfcflux_land
 integer,intent(in),optional:: flag_lsm
 integer,intent(in),optional:: spp_pbl
 integer,intent(in),optional:: ivegsrc
 integer,intent(inout),optional:: sfc_z0_type ! option for calculating surface roughness length over ocean
 logical,intent(inout),optional:: redrag ! reduced drag coeff. flag for high wind over sea (j.han)
 logical,intent(inout),optional:: flag_iter
 
 !Input data needed for GFS-related options
 integer, dimension(ims:ime,jms:jme), optional, intent(in) :: vegtype
 real(kind_phys),dimension(ims:ime,jms:jme),optional,intent(in):: &
      sigmaf,shdmax,z0pert,ztpert

 !threshold for choosing snow/ice points (In WRF, snowh is in meters)
 real,parameter:: snow_thresh = 0.05 !5 cm

 !phycical constants - now coming in through the common module
! real(kind_phys),intent(in):: svp1,svp2,svp3,svpt0
! real(kind_phys),intent(in):: ep1,ep2,karman
! real(kind_phys),intent(in):: cp,g,rovcp,r,xlv
 real(kind_phys), parameter :: &
    zero     = 0.0,   &
    one      = 1.0,   &
    two      = 2.0

 real(kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    dz8w,  &
    qv3d,  &
    p3d,   &
    t3d,   &
    u3d,   &
    v3d,   &
    rho3d, &
    th3d

 real(kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    pattern_spp_pbl

 real(kind_phys),intent(in),dimension(ims:ime,jms:jme),optional:: &
    tsurf
 
 real(kind_phys),intent(in),dimension(ims:ime,jms:jme):: &
    mavail, &
    pblh,   &
    xland,  &
    tsk,    &
    qcg,    &
    psfcpa, &
    snowh,  &
    dx

 !--- output arguments:
 character(len=*), intent(inout) :: errmsg
 integer,          intent(inout) :: errflg

 real(kind_phys),intent(out),dimension(ims:ime,jms:jme):: &
    u10,   &
    v10,   &
    th2,   &
    t2,    &
    q2

 real(kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
    cm,    &
    ck,    &
    cka,   &
    cd,    &
    cda,   &
    ustm,  &
    fm,    &
    fh,    &
    fm10,  &
    fh2,   &
    hflx,  &
    qflx,  &
    stress

!--- inout arguments:
 real(kind_phys),intent(inout),dimension(ims:ime,jms:jme):: &
    regime, &
    hfx,    &
    qfx,    &
    lh,     &
    mol,    &
    rmol,   &
    qsfc,   &
    qgh,    &
    znt,    &
    zol,    &
    ust,    &
    cpm,    &
    chs2,   &
    cqs2,   &
    cqs,    &
    chs,    &
    ch,     &
    flhc,   &
    flqc,   &
    gz1oz0, &
    wspd,   &
    br,     &
    psim,   &
    psih

!--- local variables and arrays:
 integer:: i,j,k,vegtype_1,iter,loc_z0_type
 logical:: loc_redrag,loc_iter
 real(kind_phys),dimension(ims:ime,jms:jme):: qstar,wstar

!intermediate single-point variables will be *_1
 real(kind_phys) :: mavail_1,pblh_1,xland_1,tsk_1,psfcpa_1,           &
                     qcg_1,snowh_1,dx_1
 real(kind_phys) :: u_1,v_1,u_2,v_2,qv_1,p_1,t_1,rho_1,dz8w_1,        &
                     dz8w_2,rstoch_1
 real(kind_phys) :: regime_1,hfx_1,qfx_1,lh_1,mol_1,rmol_1,           &
                     qgh_1,qsfc_1,znt_1,zol_1,ust_1,cpm_1,chs2_1,     &
                     cqs_1,cqs2_1,chs_1,ch_1,flhc_1,flqc_1,gz1oz0_1,  &
                     wspd_1,br_1,psim_1,psih_1
 real(kind_phys) :: u10_1,v10_1,th2_1,t2_1,q2_1
 real(kind_phys) :: cd_1,cda_1,ck_1,cka_1,ustm_1
 !in WRF and MPAS, these variables do not propogate outside of this driver. Use as placeholders for CCPP driver:
 real(kind_phys) :: wstar_1,qstar_1,tsurf_1,cm_1,stress_1,fm_1,fh_1,  &
                     fm10_1,fh2_1,hflx_1,qflx_1
 real(kind_phys) :: sigmaf_1,shdmax_1,z0pert_1,ztpert_1

!-----------------------------------------------------------------------------------------------------------------
 errmsg = ' '
 errflg = 0

 iter = 1
 do j = jts,jte
    do i = its,ite
       u_1      = u3d(i,kts,j)
       v_1      = v3d(i,kts,j)

       !variables that may be used prior to first calculation - sometimes depending on
       !model configuration 
       !if (initflag == 1) then
       if (itimestep== 1) then
          if (ust(i,j) .lt. 1e-4 .or. ust(i,j) .gt. 4.0) then
             ust(i,j)   = max(0.04_kind_phys*sqrt(u_1*u_1 + v_1*v_1),0.001)
          endif
          !qfx(i)   = zero
          !hfx(i)   = zero
          mol_1   = zero
          !qsfc(i)  = qv1d(i)/(1.+qv1d(i))
          qsfc_1  = qsfc(i,j) !should be available (LSM variable)
          qstar_1 = zero
       endif
       ust_1   = ust(i,j)
       stress_1= ust_1**2
       mol_1   = mol(i,j)
       qsfc_1  = qsfc(i,j)
       qstar_1 = qstar(i,j)

       !spp - input only
       if(spp_pbl==1) then
          rstoch_1 = pattern_spp_pbl(i,kts,j)
       else
          rstoch_1 = zero
       endif

       !convert 3D vars to single point - input only
       dz8w_1   = dz8w(i,kts,j)
       qv_1     = qv3d(i,kts,j)
       p_1      = p3d(i,kts,j)
       t_1      = t3d(i,kts,j)
       rho_1    = rho3d(i,kts,j)
       !2nd model level winds - for diags with high-resolution grids:
       dz8w_2   = dz8w(i,kts+1,j)
       u_2      = u3d(i,kts+1,j)
       v_2      = v3d(i,kts+1,j)

       !2d input only variables:
       mavail_1 = mavail(i,j)
       pblh_1   = pblh(i,j)
       xland_1  = xland(i,j)
       tsk_1    = tsk(i,j)
       psfcpa_1 = psfcpa(i,j)
       qcg_1    = qcg(i,j)
       snowh_1  = snowh(i,j)
       dx_1     = dx(i,j)

       !inout arguments:
       regime_1 = regime(i,j)
       hfx_1    = hfx(i,j)
       qfx_1    = qfx(i,j)
       hflx_1   = hfx_1/(cp*rho_1) !ccpp
       qflx_1   = qfx_1/rho_1      !ccpp
       lh_1     = lh(i,j)
       rmol_1   = rmol(i,j)
       qgh_1    = qgh(i,j)
       znt_1    = znt(i,j)
       zol_1    = zol(i,j)
       cpm_1    = cpm(i,j)
       chs2_1   = chs2(i,j)
       cqs2_1   = cqs2(i,j)
       cqs_1    = cqs(i,j)
       chs_1    = chs(i,j)
       ch_1     = ch(i,j)
       flhc_1   = flhc(i,j)
       flqc_1   = flqc(i,j)
       gz1oz0_1 = gz1oz0(i,j)
       wspd_1   = wspd(i,j)
       br_1     = br(i,j)
       psim_1   = psim(i,j)
       psih_1   = psih(i,j)

       !output arguments:
       u10_1    = zero
       v10_1    = zero
       th2_1    = zero
       t2_1     = zero
       q2_1     = zero
       wstar_1  = zero

       !optional configuration inputs:
       loc_redrag = .false.
       if(present(redrag)) loc_redrag = redrag  !CCPP option from GFS
       loc_iter   = .true.
       if(present(flag_iter)) loc_iter = flag_iter !CCPP option. In WRF or MPAS, just set = True 
       loc_z0_type= 0
       if(present(sfc_z0_type)) loc_z0_type = sfc_z0_type !CCPP option from GFS
       !optional array inputs:
       if(present(tsurf)) then
          tsurf_1  = tsurf(i,j)
       else
          tsurf_1  = tsk_1
       endif
       if(present(sigmaf)) then
          sigmaf_1  = sigmaf(i,j)	
       else
          sigmaf_1  = zero
       endif
       if(present(shdmax)) then
          shdmax_1  = shdmax(i,j)	
       else
          shdmax_1  = zero
       endif
       if(present(z0pert)) then
          z0pert_1  = z0pert(i,j)	
       else
          z0pert_1  = zero
       endif
       if(present(ztpert)) then
          ztpert_1  = ztpert(i,j)	
       else
          ztpert_1  = zero
       endif
       if(present(ztpert)) then
          vegtype_1 = vegtype(i,j)
       else
          vegtype_1  = 1
       endif
       
       !optional output arguments:
       if(present(ck) .and. present(cka) .and. present(cd) .and. present(cda)) then
          ck_1  = zero
          cka_1 = zero
          cd_1  = zero
          cda_1 = zero
       endif
       if(present(ustm)) then
          ustm_1 = ustm(i,j)
       else
          ustm_1 = zero
       endif
       if(present(fm) .and. present(fh) .and. present(fm10) .and. present(fh2)) then
          fm_1     = zero  !ccpp
          fh_1     = zero  !ccpp
          fm10_1   = zero  !ccpp
          fh2_1    = zero  !ccpp
       endif

       if (((xland_1-1.5) .lt. zero) .and. (snowh_1 .lt. snow_thresh)) then !5 cm threshold for binary snow/no-snow 
          call mynnsfc_land( &
                 !model info
                 i        = i         , j        = j         , itimestep= itimestep, flag_iter = loc_iter  , &
                 dx       = dx_1      , xland    = xland_1   , iter     = iter     ,                         &
                 !3d input - transformed to single point
                 u_1      = u_1       , v_1      = v_1       , t_1      = t_1      , qv_1      = qv_1      , &
                 p_1      = p_1       , dz8w_1   = dz8w_1    , rho_1    = rho_1    , u_2       = u_2       , &
                 v_2      = v_2       , dz8w_2   = dz8w_2    ,                                               &
                 !GFS input
                 sigmaf   = sigmaf_1  , vegtype  = vegtype_1 , shdmax   = shdmax_1 , ivegsrc   = ivegsrc   , &  !optional intent(in)
                 z0pert   = z0pert_1  , ztpert   = ztpert_1  , redrag   =loc_redrag,sfc_z0_type=loc_z0_type, &  !optional intent(in)
                 !2d variables - transformed to single point
                 chs      = chs_1     , chs2     = chs2_1    , cqs2     = cqs2_1   , cqs       = cqs_1     , &
                 pblh     = pblh_1    , rmol     = rmol_1    , znt      = znt_1    , psfcpa    = psfcpa_1  , &
                 ust      = ust_1     , ustm     = ustm_1    , stress   = stress_1 , &
                 mavail   = mavail_1  , zol      = zol_1     , mol      = mol_1    , tsurf     = tsurf_1   , &
                 psim     = psim_1    , psih     = psih_1    , hfx      = hfx_1    , qfx       = qfx_1     , &
                 tskin    = tsk_1     , u10      = u10_1     , v10      = v10_1    , th2       = th2_1     , &
                 t2       = t2_1      , q2       = q2_1      , flhc     = flhc_1   , flqc      = flqc_1    , &
                 snowh    = snowh_1   , qgh      = qgh_1     , qsfc     = qsfc_1   ,                         &
                 lh       = lh_1      , gz1oz0   = gz1oz0_1  , wspd     = wspd_1   , rb        = br_1      , &
                 cpm      = cpm_1     , ch       = ch_1      , cm       = cm_1     , rstoch_1  = rstoch_1  , &
                 wstar    = wstar_1   , qstar    = qstar_1   ,                                               &
                 ck       = ck_1      , cka      = cka_1     , cd       = cd_1     , cda       = cda_1     , &
                 psix     = fm_1      , psit     = fh_1      , psix10   = fm10_1   , psit2     = fh2_1     , &
                 !configuration options
                 spp_sfc  = spp_pbl   , isfflx   = isfflx    , sf_mynn_sfcflux_land= sf_mynn_sfcflux_land  , &
                 flag_restart= restart, flag_cycle= cycling  , compute_flux        = compute_flux          , &
                 psi_opt  = psi_opt   , lsm      = flag_lsm  , compute_diag        = compute_diag          , &
                 lsm_ruc  = lsm_ruc   ,                                                                      &
                 !error management
                 errmsg   = errmsg    , errflg   = errflg                                                    &
                    )
       endif

       if (((xland_1-1.5) .lt. zero) .and. (snowh_1 .ge. snow_thresh)) then !5 cm threshold for binary snow/no-snow 
          call mynnsfc_ice( &
                 !model info
                 i        = i         , j        = j         , itimestep= itimestep, flag_iter = loc_iter  , &
                 dx       = dx_1      , xland    = xland_1   , iter     = iter     ,                         &
                 !3d input - transformed to single point
                 u_1      = u_1       , v_1      = v_1       , t_1      = t_1      , qv_1      = qv_1      , &
                 p_1      = p_1       , dz8w_1   = dz8w_1    , rho_1    = rho_1    , u_2       = u_2       , &
                 v_2      = v_2       , dz8w_2   = dz8w_2    ,                                               &
                 !2d variables - transformed to single point
                 chs      = chs_1     , chs2     = chs2_1    , cqs2     = cqs2_1   , cqs       = cqs_1     , &
                 pblh     = pblh_1    , rmol     = rmol_1    , znt      = znt_1    , psfcpa    = psfcpa_1  , &
                 ust      = ust_1     , ustm     = ustm_1    , stress   = stress_1 , &
                 mavail   = mavail_1  , zol      = zol_1     , mol      = mol_1    , tsurf     = tsurf_1   , &
                 psim     = psim_1    , psih     = psih_1    , hfx      = hfx_1    , qfx       = qfx_1     , &
                 tskin    = tsk_1     , u10      = u10_1     , v10      = v10_1    , th2       = th2_1     , &
                 t2       = t2_1      , q2       = q2_1      , flhc     = flhc_1   , flqc      = flqc_1    , &
                 snowh    = snowh_1   , qgh      = qgh_1     , qsfc     = qsfc_1   ,                         &
                 lh       = lh_1      , gz1oz0   = gz1oz0_1  , wspd     = wspd_1   , rb        = br_1      , &
                 cpm      = cpm_1     , ch       = ch_1      , cm       = cm_1     , rstoch_1  = rstoch_1  , &
                 wstar    = wstar_1   , qstar    = qstar_1   ,                                               &
                 ck       = ck_1      , cka      = cka_1     , cd       = cd_1     , cda       = cda_1     , &
                 psix     = fm_1      , psit     = fh_1      , psix10   = fm10_1   , psit2     = fh2_1     , &
                 !configuration options
                 spp_sfc  = spp_pbl   , isfflx    = isfflx   ,                                               &
                 flag_restart= restart, flag_cycle= cycling  , compute_flux       = compute_flux           , &
                 psi_opt  = psi_opt   , lsm      = flag_lsm  , compute_diag       = compute_diag           , &
                 lsm_ruc  = lsm_ruc   ,                                                                      &
                 !error management
                 errmsg   = errmsg    , errflg   = errflg                                                    &
                    )
       endif

       if (((xland_1-1.5) .ge. zero)) then
          call mynnsfc_water( &
                 !model info
                 i        = i         , j        = j         , itimestep= itimestep, flag_iter = loc_iter  , &
                 dx       = dx_1      , xland    = xland_1   , iter     = iter     ,                         &
                 !3d input - transformed to single point
                 u_1      = u_1       , v_1      = v_1       , t_1      = t_1      , qv_1      = qv_1      , &
                 p_1      = p_1       , dz8w_1   = dz8w_1    , rho_1    = rho_1    , u_2       = u_2       , &
                 v_2      = v_2       , dz8w_2   = dz8w_2    ,                                               &
                 !GFS input
                 sigmaf   = sigmaf_1  , vegtype  = vegtype_1 , shdmax   = shdmax_1 , ivegsrc   = ivegsrc   , &  !optional intent(in)
                 z0pert   = z0pert_1  , ztpert   = ztpert_1  , redrag   =loc_redrag, sfc_z0_type=loc_z0_type,&  !optional intent(in)
                 !2d variables - transformed to single point
                 chs      = chs_1     , chs2     = chs2_1    , cqs2     = cqs2_1   , cqs       = cqs_1     , &
                 pblh     = pblh_1    , rmol     = rmol_1    , znt      = znt_1    , psfcpa    = psfcpa_1  , &
                 ust      = ust_1     , ustm     = ustm_1    , stress   = stress_1 , &
                 mavail   = mavail_1  , zol      = zol_1     , mol      = mol_1    , tsurf     = tsurf_1   , &
                 psim     = psim_1    , psih     = psih_1    , hfx      = hfx_1    , qfx       = qfx_1     , &
                 tskin    = tsk_1     , u10      = u10_1     , v10      = v10_1    , th2       = th2_1     , &
                 t2       = t2_1      , q2       = q2_1      , flhc     = flhc_1   , flqc      = flqc_1    , &
                 snowh    = snowh_1   , qgh      = qgh_1     , qsfc     = qsfc_1   ,                         &
                 lh       = lh_1      , gz1oz0   = gz1oz0_1  , wspd     = wspd_1   , rb        = br_1      , &
                 cpm      = cpm_1     , ch       = ch_1      , cm       = cm_1     , rstoch_1  = rstoch_1  , &
                 wstar    = wstar_1   , qstar    = qstar_1   ,                                               &
                 ck       = ck_1      , cka      = cka_1     , cd       = cd_1     , cda       = cda_1     , &
                 psix     = fm_1      , psit     = fh_1      , psix10   = fm10_1   , psit2     = fh2_1     , &
                 !configuration options
                 spp_sfc  = spp_pbl   , isfflx   = isfflx    ,sf_mynn_sfcflux_water= sf_mynn_sfcflux_water , &
                 flag_restart= restart, flag_cycle= cycling  , compute_flux        = compute_flux          , &
                 psi_opt  = psi_opt   ,                        compute_diag        = compute_diag          , &
                 lsm_ruc  = lsm_ruc   , lsm      = flag_lsm  ,                                               &
                 !error management
                 errmsg   = errmsg    , errflg   = errflg                                                    &
                    )
       endif
       
       !inout arguments:
       regime(i,j) = regime_1
       hfx(i,j)    = hfx_1
       qfx(i,j)    = qfx_1
       lh(i,j)     = lh_1
       mol(i,j)    = mol_1
       rmol(i,j)   = rmol_1
       qgh(i,j)    = qgh_1
       qsfc(i,j)   = qsfc_1
       znt(i,j)    = znt_1
       zol(i,j)    = zol_1
       ust(i,j)    = ust_1
       cpm(i,j)    = cpm_1
       chs2(i,j)   = chs2_1
       cqs2(i,j)   = cqs2_1
       cqs(i,j)    = cqs_1
       chs(i,j)    = chs_1
       ch(i,j)     = ch_1
       flhc(i,j)   = flhc_1
       flqc(i,j)   = flqc_1
       gz1oz0(i,j) = gz1oz0_1
       wspd(i,j)   = wspd_1
       br(i,j)     = br_1
       psim(i,j)   = psim_1
       psih(i,j)   = psih_1

       !output arguments:
       u10(i,j)    = u10_1
       v10(i,j)    = v10_1
       th2(i,j)    = th2_1
       t2(i,j)     = t2_1
       q2(i,j)     = q2_1
       wstar(i,j)  = wstar_1
       qstar(i,j)  = qstar_1

       !optional output arguments, mostly for ccpp:
       if(present(hflx))hflx(i,j)= hfx_1/(cp*rho_1)
       if(present(qflx))qflx(i,j)= qfx_1/rho_1
       if(present(ck) .and. present(cka) .and. present(cd) .and. present(cda)) then
          ck(i,j)  = ck_1
          cka(i,j) = cka_1
          cd(i,j)  = cd_1
          cda(i,j) = cda_1
       endif
       if(present(fm) .and. present(fh) .and. present(fm10) .and. present(fh2)) then
          fm(i,j)    = fm_1
          fh(i,j)    = fh_1
          fm10(i,j)  = fm10_1
          fh2(i,j)   = fh2_1
       endif
       if(present(ustm))  ustm(i,j)  = ustm_1
       if(present(stress))stress(i,j)= stress_1
       if(present(cm))    cm(i,j)    = cm_1
       
    enddo !i-loop
 enddo !j-loop
 

 end subroutine mynnsfc_driver

!=================================================================================================================
 end module module_sf_mynnsfc_driver
!=================================================================================================================
