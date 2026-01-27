! Module for MYNN SFC scheme tests
module module_sf_mynnsfc_wrf_tests
    use module_sf_mynnsfc_driver
    ! public
    !=================================================================================================================    
    implicit none
    logical :: cycling,restart
    integer :: initflag, spp_pbl
    real,dimension(1) :: pattern_spp_pbl

    contains

    subroutine init_mynn_sfc_flags()
      write(*,*) '--- calling  init_mynn_sfc_flags()'
       cycling=.false.
       initflag=1
       spp_pbl=1
       restart=.false.
       pattern_spp_pbl=0.0

      ! for future use
    end subroutine init_mynn_sfc_flags
    !=================================================================================================================    
    subroutine init_input_data_for_test()
      integer :: iostat, line_num
      character(len=2000) :: input_line
      integer, parameter :: input_unit = 10
      integer, parameter :: output_unit = 20
      
      write(*,*) '--- opening data files ---'
      ! Open input file
      close(unit=input_unit)
      open(unit=input_unit, file='./data/input_lnd.txt', status='old', action='read', iostat=iostat)

      if (iostat /= 0) then
          print *, 'Error opening input file'
          stop
      end if

      ! Open output file
      open(unit=output_unit, file='./data/wrf_output_lnd.txt', status='replace', action='write', iostat=iostat)
      write(output_unit,'(A5, A5, A10, A10, A10, A10, A10, A10, A10, A10, A10)')                &
            'itimestep', 'iter', 'T2', 'Q2', 'TH2', 'U10', 'V10', 'HFX', 'LH', 'UST_lnd','PBLH'
      if (iostat /= 0) then
          print *, 'Error opening output file'
          close(output_unit)
          stop
      end if

    end subroutine init_input_data_for_test

    !===============================================================================
    ! Subroutine to process each line
    !===============================================================================
    subroutine process_line(line, out_unit, line_number, flag_iter,  &                         
           U1D,V1D,T1D,QV1D,P1D,dz8w1d,                              &
           U1D2,V1D2,dz2w1d,                                         &
           PSFCPA,PBLH,MAVAIL,XLAND,DX,                              &
           ISFFLX,isftcflx,iz0tlnd,psi_opt,                          &
           compute_flux,compute_diag,                                &
           sigmaf,vegtype,shdmax,ivegsrc,                            & 
           z0pert,ztpert,                                            &
           redrag,sfc_z0_type,                                       &
           itimestep,iter,flag_restart,lsm,lsm_ruc,                  &
           wet,          dry,          icy,                          &
           tskin_wat,    tskin_lnd,    tskin_ice,                    &
           tsurf_wat,    tsurf_lnd,    tsurf_ice,                    &
           qsfc_wat,     qsfc_lnd,     qsfc_ice,                     &
           snowh_wat,    snowh_lnd,    snowh_ice,                    &
           ZNT_wat,      ZNT_lnd,      ZNT_ice,                      &
           UST_wat,      UST_lnd,      UST_ice,                      &
           cm_wat,       cm_lnd,       cm_ice,                       &
           ch_wat,       ch_lnd,       ch_ice,                       &
           rb_wat,       rb_lnd,       rb_ice,                       &
           stress_wat,   stress_lnd,   stress_ice,                   &
           psix_wat,     psix_lnd,     psix_ice,                     &
           psit_wat,     psit_lnd,     psit_ice,                     &
           psix10_wat,   psix10_lnd,   psix10_ice,                   &
           psit2_wat,    psit2_lnd,    psit2_ice,                    &
           HFLX_wat,     HFLX_lnd,     HFLX_ice,                     &
           QFLX_wat,     QFLX_lnd,     QFLX_ice,                     &
           ch,CHS,CHS2,CQS2,CPM,                                     &
           ZNT,USTM,ZOL,MOL,RMOL,                                    &
           PSIM,PSIH,                                                &
           HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                           &
           QGH,QSFC,                                                 &
           U10,V10,TH2,T2,Q2,                                        &
           GZ1OZ0,WSPD,wstar,qstar,                                  &
           spp_sfc,  rstoch1D )         
        
        implicit none
        
        character(len=2000), intent(in) :: line
        integer, intent(in) :: out_unit
        integer, intent(in) :: line_number
        
        ! Variables to store parsed values
        logical, intent(out) :: flag_iter, compute_flux, compute_diag, redrag, flag_restart
        logical, intent(out) :: wet, dry, icy
        real, intent(out) :: U1D, V1D, T1D, QV1D, P1D, dz8w1d, U1D2, V1D2, dz2w1d, &
                PSFCPA, PBLH, MAVAIL, XLAND, DX, sigmaf, shdmax, z0pert,           &
                ztpert
        real, intent(out) :: tskin_wat, tskin_lnd, tskin_ice, tsurf_wat, tsurf_lnd, tsurf_ice, &
                qsfc_wat, qsfc_lnd, qsfc_ice, snowh_wat, snowh_lnd, snowh_ice,                 &
                ZNT_wat, ZNT_lnd, ZNT_ice, UST_wat, UST_lnd, UST_ice,                          &
                cm_wat, cm_lnd, cm_ice, ch_wat, ch_lnd, ch_ice,                                &
                rb_wat, rb_lnd, rb_ice, stress_wat, stress_lnd, stress_ice,                    &                   
                psix_wat,     psix_lnd,     psix_ice,                                          &
                psit_wat,     psit_lnd,     psit_ice,                                          &
                psix10_wat,   psix10_lnd,   psix10_ice,                                        &
                psit2_wat,    psit2_lnd,    psit2_ice,                                         &
                HFLX_wat,     HFLX_lnd,     HFLX_ice,                                          &
                QFLX_wat,     QFLX_lnd,     QFLX_ice,                                          &
                ch,CHS,CHS2,CQS2,CPM,                                                          &
                ZNT,USTM,ZOL,MOL,RMOL,                                                         &
                PSIM,PSIH,                                                                     &
                HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                                                &
                QGH,QSFC,                                                                      &
                U10,V10,TH2,T2,Q2,                                                             &
                GZ1OZ0,WSPD,wstar,qstar,                                                       &
                rstoch1D  
        integer :: read_stat           
        integer, intent(out) :: ISFFLX, isftcflx, iz0tlnd, psi_opt, vegtype,                   &
                   ivegsrc, sfc_z0_type, itimestep, iter, lsm, lsm_ruc, spp_sfc
        
        ! Read values from the line with specified format
        read (line, '(L5,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,'          // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,' // &
                     'I5,I5,I5,I5,'                                     // &
                     'L5,L5,'                                           // &
                     'E15.3,I5,E15.3,I5,'                               // &
                     'E15.3,E15.3,'                                     // &
                     'L5,I5,'                                           // &
                     'I5,I5,L5,I5,I5,'                                  // &
                     'L5,L5,L5,'                                        // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // & !stress_* values
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,'                   // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,'                   // &
                     'E15.3,E15.3,'                                     // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,'       // &
                     'E15.3,E15.3,'                                     // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,'                   // &
                     'E15.3,E15.3,E15.3,E15.3,'                         // &
                     'I5,E15.3)', iostat=read_stat)                        &
                      flag_iter,                                           &
                      U1D,   V1D,   T1D,   QV1D,   P1D,   dz8w1d,          &
                      U1D2,   V1D2,   dz2w1d,                              &
                      PSFCPA,   PBLH,   MAVAIL,   XLAND,   DX,             &
                      ISFFLX,isftcflx,iz0tlnd,psi_opt,                     &
                      compute_flux,compute_diag,                           &
                      sigmaf,   vegtype,   shdmax,   ivegsrc,              &
                      z0pert,   ztpert,                                    &
                      redrag,sfc_z0_type,                                  &
                      itimestep,iter,flag_restart,lsm,lsm_ruc,             &
                            wet,             dry,             icy,         &
                      tskin_wat,       tskin_lnd,       tskin_ice,         &
                      tsurf_wat,       tsurf_lnd,       tsurf_ice,         &
                       qsfc_wat,        qsfc_lnd,        qsfc_ice,         &
                      snowh_wat,       snowh_lnd,       snowh_ice,         &
                        ZNT_wat,         ZNT_lnd,         ZNT_ice,         &
                        UST_wat,         UST_lnd,         UST_ice,         &
                         cm_wat,          cm_lnd,          cm_ice,         &
                         ch_wat,          ch_lnd,          ch_ice,         &
                         rb_wat,          rb_lnd,          rb_ice,         &
                     stress_wat,      stress_lnd,      stress_ice,         &
                     psix_wat,     psix_lnd,     psix_ice,                 & 
                     psit_wat,     psit_lnd,     psit_ice,                 & 
                     psix10_wat,   psix10_lnd,   psix10_ice,               &
                     psit2_wat,    psit2_lnd,    psit2_ice,                &
                     HFLX_wat,     HFLX_lnd,     HFLX_ice,                 &
                     QFLX_wat,     QFLX_lnd,     QFLX_ice,                 &
                     ch,CHS,CHS2,CQS2,CPM,                                 &
                     ZNT,USTM,ZOL,MOL,RMOL,                                &
                     PSIM,PSIH,                                            &
                     HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                       &
                     QGH,QSFC,                                             &
                     U10,V10,TH2,T2,Q2,                                    &
                     GZ1OZ0,WSPD,wstar,qstar,                              &
                     spp_sfc,rstoch1D                                              
              
      ! Debug: print the line being read
      write(0,*) "Line length:", len_trim(line)
      write(0,*) "U1D=",U1D, "flag_iter=",flag_iter, "V1D=", V1D," ZNT", ZNT_lnd, "WSPD=",WSPD, "UST=",UST_lnd, &
                "HFX=",HFX, "LH=",LH
        
    end subroutine process_line

    !=================================================================================================================
 
    subroutine wrf_test

      !use module_sf_mynnsfc_land, only : mynnsfc_land
      ! use module_sf_mynnsfc_driver, only: mynnsfc_init

      integer :: iostat, line_num
      integer, parameter :: i=1, j=1

      character(len=2000) :: input_line
      integer, parameter :: input_unit = 10
      integer, parameter :: output_unit = 20

      logical :: flag_iter 
      logical :: compute_flux, compute_diag, redrag, flag_restart
      logical, dimension(1)  :: wet, dry, icy
      real, dimension(1) ::  U1D, V1D, T1D, QV1D, P1D, dz8w1d, U1D2, V1D2, dz2w1d, &
              PSFCPA, PBLH, MAVAIL, XLAND, DX, sigmaf, shdmax, z0pert,             &
              ztpert
      integer :: spp_sfc
      real, dimension(1) :: tskin_wat, tskin_lnd, tskin_ice, tsurf_wat, tsurf_lnd, tsurf_ice, &
              qsfc_wat, qsfc_lnd, qsfc_ice, snowh_wat, snowh_lnd, snowh_ice,                  &
              ZNT_wat, ZNT_lnd, ZNT_ice, UST_wat, UST_lnd, UST_ice,                           &
              cm_wat, cm_lnd, cm_ice, ch_wat, ch_lnd, ch_ice,                                 &
              rb_wat, rb_lnd, rb_ice, stress_wat, stress_lnd, stress_ice,                     &
                   psix_wat,     psix_lnd,     psix_ice,                                      &
                   psit_wat,     psit_lnd,     psit_ice,                                      &
                 psix10_wat,   psix10_lnd,   psix10_ice,                                      &
                  psit2_wat,    psit2_lnd,    psit2_ice,                                      &
                   HFLX_wat,     HFLX_lnd,     HFLX_ice,                                      &
                   QFLX_wat,     QFLX_lnd,     QFLX_ice,                                      &
                 ch,CHS,CHS2,CQS2,CPM,                                                        &
                 ZNT,USTM,ZOL,MOL,RMOL,                                                       &
                 PSIM,PSIH,                                                                   &
                 HFLX_IN,HFX_IN,QFLX_IN,QFX_IN,LH_IN,FLHC,FLQC,                               &
                 QGH,QSFC,                                                                    &
                 U10,V10,TH2,T2,Q2,                                                           &
                 GZ1OZ0,WSPD,wstar,qstar,                                                     &
                 rstoch1D

      ! Free evolving variables           
      real, dimension(1) :: HFLX, HFX, QFLX, QFX, LH, stress, ust
      integer, dimension(1) :: vegtype
                  
      integer :: read_stat,ISFFLX,isftcflx,iz0tlnd,psi_opt,ivegsrc,sfc_z0_type,itimestep,iter,lsm,lsm_ruc

      ! Variables for error handling
      character(len=512) :: errmsg
      integer :: errflg
      ! Initialize 3D variables
      real, dimension(64) ::  u3d, v3d, t3d, qv3d, p3d, dz8w, th3d, rho3d  ! 64 is random number
      u3d=0.0
      v3d=0.0
      t3d=0.0
      qv3d=0.0
      p3d=0.0
      dz8w=0.0 
      th3d=0.0
      rho3d=0.0


      ! Initialize error variables
      errmsg = ''
      errflg = 0

      write(*,*) '--- entering wrf_test subroutine ---'    
      ! Initialize input data for tests
      call init_input_data_for_test()

      ! Read header
      read(input_unit, '(A)', iostat=iostat) input_line

      ! Process each line
      line_num = 0
      do
          read(input_unit, '(A)', iostat=iostat) input_line
          write(0,*) input_line
          
          ! Check for end of file or error
          if (iostat < 0) exit  ! End of file
          if (iostat > 0) then
              print *, 'Error reading line', line_num + 1
              exit
          end if
          
          line_num = line_num + 1
          
          ! Call subroutine to process the line
          call process_line(input_line, output_unit, line_num, flag_iter,      &
           U1D(1),  V1D(1),  T1D(1),  QV1D(1),  P1D(1),  dz8w1d(1),            &
           U1D2(1),  V1D2(1),  dz2w1d(1),                                      &
           PSFCPA(1),  PBLH(1),  MAVAIL(1),  XLAND(1),  DX(1),                 &
           ISFFLX,  isftcflx,  iz0tlnd,  psi_opt,                              &
           compute_flux,  compute_diag,                                        &
           sigmaf(1),  vegtype(1),  shdmax(1),  ivegsrc,                       &
           z0pert(1),  ztpert(1),                                              &
           redrag,  sfc_z0_type,                                               &
           itimestep,  iter,  flag_restart,  lsm,  lsm_ruc,                    &
                  wet(1),            dry(1),            icy(1),                &
            tskin_wat(1),      tskin_lnd(1),      tskin_ice(1),                &
            tsurf_wat(1),      tsurf_lnd(1),      tsurf_ice(1),                &
             qsfc_wat(1),       qsfc_lnd(1),       qsfc_ice(1),                &
            snowh_wat(1),      snowh_lnd(1),      snowh_ice(1),                &
              ZNT_wat(1),        ZNT_lnd(1),        ZNT_ice(1),                &
              UST_wat(1),        UST_lnd(1),        UST_ice(1),                &
               cm_wat(1),         cm_lnd(1),         cm_ice(1),                &
               ch_wat(1),         ch_lnd(1),         ch_ice(1),                &
               rb_wat(1),         rb_lnd(1),         rb_ice(1),                &
           stress_wat(1),     stress_lnd(1),     stress_ice(1),                &
             psix_wat(1),       psix_lnd(1),       psix_ice(1),                &
             psit_wat(1),       psit_lnd(1),       psit_ice(1),                &
             psix10_wat(1),     psix10_lnd(1),     psix10_ice(1),              &
             psit2_wat(1),      psit2_lnd(1),      psit2_ice(1),               &
             HFLX_wat(1),       HFLX_lnd(1),       HFLX_ice(1),                &
             QFLX_wat(1),       QFLX_lnd(1),       QFLX_ice(1),                &
             ch(1),  CHS(1),  CHS2(1),  CQS2(1),  CPM(1),                      &
             ZNT(1),  USTM(1),  ZOL(1),  MOL(1),  RMOL(1),                     &
             PSIM(1),  PSIH(1),                                                &
             HFLX_IN(1),  HFX_IN(1),  QFLX_IN(1),  QFX_IN(1),                  &
             LH_IN(1),  FLHC(1),  FLQC(1),                                     &
             QGH(1),  QSFC(1),                                                 &
             U10(1),  V10(1),  TH2(1),  T2(1),  Q2(1),                         &
             GZ1OZ0(1),  WSPD(1),  wstar(1),  qstar(1),                        &
             spp_sfc,  rstoch1D(1))

          ! Construct 3D arrays
          u3d(1)=U1D(1)
          u3d(2)=U1D2(1)
          v3d(1)=V1D(1)
          v3d(2)=V1D2(1)
          t3d(1)=T1D(1)
          qv3d(1)=QV1D(1)
          p3d(1)=P1D(1)
          dz8w(1)=dz8w1d(1)
          dz8w(2)=dz2w1d(1)

          ! write(*,*) 'u3d ',u3d

         ! Initialize MYNN SFC
          write(*,*) '--- calling  mynnsfc_land() ---'

          call mynnsfc_init(allowed_to_read=.true.,errmsg=errmsg,errflg=errflg)
          !write(*,*) 'psih_stab_in_wrftests=',psih_stab

          call mynnsfc_driver(u3d=u3d , v3d =v3d , t3d=t3d , qv3d =qv3d, p3d =p3d, dz8w=dz8w,        &
              th3d=th3d, rho3d=rho3d,                                                                &
              !GFS-related input
              sigmaf=sigmaf, vegtype=vegtype, shdmax=shdmax, ivegsrc=ivegsrc,                        &  !intent(in)
              z0pert=z0pert, ztpert=ztpert, redrag=redrag, sfc_z0_type=sfc_z0_type,                  &  !intent(in)
              !2d variables
              psfcpa=psfcpa , chs=CHS, chs2=CHS2, cqs=CHS2 , cqs2=CHS2, cpm=CPM,                     &
              znt=ZNT_lnd, ust=ust, ustm=USTM, pblh=pblh, mavail=mavail, zol=ZOL ,                   &
              mol=MOL, rmol=RMOL, psim=PSIM , psih=PSIH, xland=XLAND,                                &
              hfx=hfx, qfx=QFX , lh=LH, tsk=tskin_lnd, flhc=FLHC, flqc=FLQC,                         &
              qsfc=QSFC, u10=U10, v10=V10, th2 =TH2, t2=T2 ,                                         &
              q2=Q2, snowh=snowh_lnd, gz1oz0=GZ1OZ0, wspd=WSPD, br=rb_lnd, dx=DX,                    &
              ch=ch_lnd, ck=CQS2, cka=CQS2, cd=CQS2, cda=CQS2 ,                                      &
              stress=stress , hflx=hflx, qflx=qflx, cm=cm_lnd, fm=psix_lnd, fh=psit_lnd,             &
              fm10=psix10_lnd, fh2=psit2_lnd,                                                        &
              tsurf=tsurf_lnd    ,                                                                   &
              !configuration options
              spp_pbl=spp_pbl, pattern_spp_pbl=pattern_spp_pbl,                                                 &
              sf_mynn_sfcflux_water=iz0tlnd          ,                                               &
              sf_mynn_sfcflux_land=iz0tlnd           ,                                               &
              isfflx=ISFFLX, restart=restart, cycling=cycling, initflag=1,                           &
              flag_iter=flag_iter, flag_lsm=lsm,                                                     &
              !model information
              itimestep=itimestep,                                                                   &
              ids=1     , ide=1      , jds=1    , jde=1     , kds=1     , kde=1     ,                &
              ims=1     , ime=1      , jms=1    , jme=1     , kms=1     , kme=1     ,                &
              its=1     , ite=1      , jts=1    , jte=1     , kts=1     , kte=1     ,                &
              errmsg=errmsg , errflg=errflg                                                          &
              )


         write(0,*) "T2=",t2,'chs=',ch,'ust=',UST_lnd,'hfx=',hfx
         write(0,*) "Read status:", read_stat
         open(output_unit, file = './data/wrf_output_lnd.txt')
         write(output_unit,'(I5, I5, F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)')  &
              itimestep,iter,t2,q2,th2,u10,v10,hfx,lh,ust,pblh

      end do
    end subroutine wrf_test

end module module_sf_mynnsfc_wrf_tests           
