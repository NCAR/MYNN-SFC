! ###########################################################################################
! CI test driver for MYNN SFC scheme
! ###########################################################################################
program driver
  use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
  use module_sf_mynnsfc_common
#if (EM_CORE == 1)
  use module_sf_mynnsfc_wrf_tests
#else
  use module_sf_mynnsfc_ccpp_tests
#endif

  implicit none

#if (EM_CORE == 1)
  call wrf_test(case='wat', sf_mynn_sfcflux_water=1, sf_mynn_sfcflux_land=1, saveoutput=.true.)
  call wrf_test(case='wat', sf_mynn_sfcflux_water=0, sf_mynn_sfcflux_land=1, saveoutput=.false.)
  call wrf_test(case='wat', sf_mynn_sfcflux_water=2, sf_mynn_sfcflux_land=1, saveoutput=.false.)
  call wrf_test(case='wat', sf_mynn_sfcflux_water=3, sf_mynn_sfcflux_land=1, saveoutput=.false.)
  call wrf_test(case='wat', sf_mynn_sfcflux_water=4, sf_mynn_sfcflux_land=1, saveoutput=.false.)

  call wrf_test(case='lnd', sf_mynn_sfcflux_water=1, sf_mynn_sfcflux_land=1, saveoutput=.true.)
  call wrf_test(case='lnd', sf_mynn_sfcflux_water=1, sf_mynn_sfcflux_land=0, saveoutput=.false.)
  call wrf_test(case='lnd', sf_mynn_sfcflux_water=1, sf_mynn_sfcflux_land=2, saveoutput=.false.)

  call wrf_test(case='icy', sf_mynn_sfcflux_water=1, sf_mynn_sfcflux_land=0, saveoutput=.true.)
#else
  call ccpp_test(saveoutput=.true.)
#endif
end program driver
