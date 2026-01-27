! ###########################################################################################
! CI test driver for MYNN SFC scheme
! ###########################################################################################
program driver
  use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
  use module_sf_mynnsfc_common,only: kind_phys,cp,lsm_ruc => ruclsmscheme
  use module_sf_mynnsfc_ccpp_tests
  use module_sf_mynnsfc_wrf_tests
  implicit none
  call ccpp_test()
  call wrf_test()
end program driver
