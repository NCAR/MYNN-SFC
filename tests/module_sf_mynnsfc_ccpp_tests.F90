! Module for MYNN SFC scheme tests in CCPP
module module_sf_mynnsfc_ccpp_tests
  implicit none

  contains
    !=================================================================================================================    

    subroutine init_mynn_sfc_flags_for_test_all_true()
      write(*,*) '--- calling  init_mynn_sfc_flags_for_test_all_true()'
      ! for future use
    end subroutine init_mynn_sfc_flags_for_test_all_true
    !=================================================================================================================    
    subroutine init_input_data_for_test(saveoutput)
      logical, intent(in) :: saveoutput
      integer :: iostat, line_num
      character(len=2000) :: input_line
      integer, parameter :: input_unit = 10
      integer, parameter :: output_unit = 20
      
      write(*,*) '--- opening data files ---'
      ! Open input file
      open(unit=input_unit, file='./data/ccpp_input_lnd.txt', status='old', action='read', iostat=iostat)

      if (iostat /= 0) then
          print *, 'Error opening input file'
          stop
      end if

      ! Open output file
      if (saveoutput) then
          open(unit=output_unit, file='./data/ccpp_output_lnd.txt', status='replace', action='write', iostat=iostat)
          write(output_unit,'(A5, A5, A10, A10, A10, A10, A10, A10, A10, A10, A10)')                &
                'itimestep', 'iter', 'T2', 'Q2', 'TH2', 'U10', 'V10', 'HFX', 'LH', 'UST_lnd','PBLH'
          if (iostat /= 0) then
              print *, 'Error opening output file'
              close(output_unit)
              stop
          end if
      end if

    end subroutine init_input_data_for_test


    !=================================================================================================================
    subroutine ccpp_test(saveoutput)

      use module_sf_mynnsfc_driver, only : mynnsfc_driver
      ! use module_sf_mynnsfc, only : SFCLAY1D_mynn

      logical, intent(in) :: saveoutput
      integer :: iostat, line_num
      integer, parameter :: ids=1, ide=1, jds=1, jde=1, kds=1,kde=2
      integer, parameter :: ims=1, ime=1, jms=1, jme=1, kms=1,kme=2
      integer, parameter :: its=1, ite=1, jts=1, jte=1, kts=1,kte=2

      character(len=2000) :: input_line
      integer, parameter :: input_unit = 10
      integer, parameter :: output_unit = 20
      integer, parameter :: n = 1  ! Number of points

      logical :: flag_iter 
      logical :: compute_flux, compute_diag, redrag, flag_restart

      logical, dimension(ims:ime,jms:jme):: wet, dry, icy
      real, dimension(ims:ime,kms:kme,jms:jme) :: u3d, v3d, t3d, qv3d, p3d, dz8w, th3d, rho3d
      real, dimension(ims:ime,jms:jme) :: sigmaf, shdmax, z0pert, ztpert
      integer, dimension(ims:ime,jms:jme) :: vegtype
      real, dimension(ims:ime,jms:jme) :: psfcpa, chs, chs2, cqs, cqs2, cpm, znt, ust, ustm,  &
              pblh, mavail, zol, mol, rmol, psim, psih, xland, qgh, hfx, qfx, lh, tsk, flhc, &
              flqc, qsfc, u10, v10, th2, t2, lakemask, q2, snowh, gz1oz0, wspd, br, dx, ch, ck,   &
              cka, cd, cda, stress
      real, dimension(ims:ime,jms:jme) :: tskin_wat, tskin_lnd, tskin_ice, tsurf_wat, tsurf_lnd, tsurf_ice, &
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
                 HFLX,QFLX,      wstar,qstar,  rstoch1D 

      integer :: ivegsrc, read_stat,ISFFLX,isftcflx,iz0tlnd,psi_opt,itimestep, iter,          &
            lsm, lsm_ruc,spp_sfc, sfc_z0_type
      character(len=2000) :: line
      integer :: out_unit
      integer :: line_number, i

      ! Variables for error handling
      character(len=512) :: errmsg
      integer :: errflg
      ! Initialize error variables
      errmsg = ''
      errflg = 0

      write(*,*) '--- entering ccpp_test subroutine'    
      ! Initialize input data for tests
      call init_input_data_for_test(saveoutput)

      ! Read header
      read(input_unit, '(A)', iostat=iostat) input_line

      ! Process each line
      line_num = 0
      i=0
      do
          read(input_unit, '(A)', iostat=iostat) input_line
          !write(0,*) input_line
          
          ! Check for end of file or error
          if (iostat < 0) exit  ! End of file
          if (iostat > 0) then
              print *, 'Error reading line', line_num + 1
              exit
          end if
          
          line_num = line_num + 1
          i=i+1

          read (input_line, '(L5,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,'          // &
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
                        u3d(1,1,1), v3d(1,1,1), t3d(1,1,1),qv3d(1,1,1), p3d(1,1,1), dz8w(1,1,1),&
                        u3d(1,2,1), v3d(1,2,1), dz8w(1,2,1),                              &
                        psfcpa(1,1),pblh(1,1),  mavail(1,1), xland(1,1), dx(1,1),             &
                        isfflx, isftcflx,iz0tlnd,psi_opt,                           &
                        compute_flux,compute_diag,                                  &
                        sigmaf(1,1),  vegtype(1,1),   shdmax(1,1),   ivegsrc,          &
                        z0pert(1,1),   ztpert(1,1),                                     &
                        redrag,sfc_z0_type,                                   &
                        itimestep,iter,flag_restart,lsm,lsm_ruc,                    &
                              wet(1,1),             dry(1,1),             icy(1,1),       &
                        tskin_wat(1,1),       tskin_lnd(1,1),       tskin_ice(1,1),       &
                        tsurf_wat(1,1),       tsurf_lnd(1,1),       tsurf_ice(1,1),         &
                         qsfc_wat(1,1),        qsfc_lnd(1,1),        qsfc_ice(1,1),         &
                        snowh_wat(1,1),       snowh_lnd(1,1),       snowh_ice(1,1),         &
                          ZNT_wat(1,1),         ZNT_lnd(1,1),         ZNT_ice(1,1),         &
                          UST_wat(1,1),         UST_lnd(1,1),         UST_ice(1,1),         &
                           cm_wat(1,1),          cm_lnd(1,1),          cm_ice(1,1),         &
                           ch_wat(1,1),          ch_lnd(1,1),          ch_ice(1,1),         &
                           rb_wat(1,1),          rb_lnd(1,1),          rb_ice(1,1),         &
                       stress_wat(1,1),      stress_lnd(1,1),      stress_ice(1,1),         &
                       psix_wat(1,1),     psix_lnd(1,1),     psix_ice(1,1),                 & 
                       psit_wat(1,1),     psit_lnd(1,1),     psit_ice(1,1),                 & 
                       psix10_wat(1,1),   psix10_lnd(1,1),   psix10_ice(1,1),               &
                       psit2_wat(1,1),    psit2_lnd(1,1),    psit2_ice(1,1),                &
                       HFLX_wat(1,1),     HFLX_lnd(1,1),     HFLX_ice(1,1),                 &
                       QFLX_wat(1,1),     QFLX_lnd(1,1),     QFLX_ice(1,1),                 &
                       ch(1,1),CHS(1,1),CHS2(1,1),CQS2(1,1),CPM(1,1),                                 &
                       ZNT(1,1),USTM(1,1),ZOL(1,1),MOL(1,1),RMOL(1,1),                                &
                       PSIM(1,1),PSIH(1,1),                                            &
                       HFLX(1,1),HFX(1,1),QFLX(1,1),QFX(1,1),LH(1,1),FLHC(1,1),FLQC(1,1),                       &
                       QGH(1,1),QSFC(1,1),                                             &
                       U10(1,1),V10(1,1),TH2(1,1),T2(1,1),Q2(1,1),                                    &
                       GZ1OZ0(1,1),WSPD(1,1),wstar(1,1),qstar(1,1),                              &
                       spp_sfc,rstoch1D(1,1)             

       ! calculate th3d and rho3d, which is not available in the ccpp input
        th3d(1,1,1)=t3d(1,1,1)*(1.0e5/p3d(1,1,1))**(287.04/1004.5)
        rho3d(1,1,1)=p3d(1,1,1)/(287.04*t3d(1,1,1)*(1.0+0.608*qv3d(1,1,1)))
        write(0,*) "t3d=",t3d(1,1,1),"th3d=",th3d(1,1,1), "lsm=", lsm
       ! call mynnsfc_driver
        write(*,*) '--- calling  mynnsfc_driver()'
        if (dry(1,1)) then
            CALL mynnsfc_driver(                                                  &
                 u3d=u3d,v3d=v3d,t3d=t3d,qv3d=qv3d,p3d=p3d,dz8w=dz8w,             &
                 th3d=th3d,rho3d=rho3d,                                           &
                 PSFCPA=PSFCPA,PBLH=pblh,MAVAIL=mavail,XLAND=xland,DX=dx,         &
                 ISFFLX=isfflx,sf_mynn_sfcflux_water=isftcflx,flag_lsm=lsm,       &
                 sf_mynn_sfcflux_land=iz0tlnd,                                    &
                 sigmaf=sigmaf,vegtype=vegtype,shdmax=shdmax,ivegsrc=ivegsrc,     &  !intent(in)
                 z0pert=z0pert,ztpert=ztpert,                                     &  !intent(in)
                 redrag=redrag,sfc_z0_type=sfc_z0_type,                           &  !intent(in)
                 itimestep=itimestep,flag_iter=flag_iter,                         &
                 restart=flag_restart,dry=dry(1,1),wet=wet(1,1),icy=icy(1,1),     &  !intent(in)
                 tsk=tskin_lnd, tsurf=tsurf_lnd, qsfc=qsfc_lnd, snowh=snowh_lnd,  &
                 znt=znt_lnd, ust=ust_lnd,cm=cm_lnd,ch=ch_lnd, br=rb_lnd,         &  !intent(inout)
                 stress=stress_lnd,fm=psix_lnd, fh=psit_lnd,fm10=psix10_lnd,      &  !intent(inout)
                 fh2=psit2_lnd, hflx=hflx_lnd, qflx=qflx_lnd, CHS=chs,CHS2=chs2,  &
                 CQS2=cqs2,CQS=cqs,CPM=cpm,USTM=ustm,ZOL=zol,MOL=mol,RMOL=rmol,   &
                 psim=psim,psih=psih,HFX=hfx,QFX=qfx,LH=lh,FLHC=flhc,FLQC=flqc,   &
                 QGH=qgh,U10=u10,V10=v10,TH2=th2,T2=t2,Q2=q2,                     &
                 GZ1OZ0=GZ1OZ0,WSPD=wspd,                                         &
                 spp_pbl=spp_sfc,pattern_spp_pbl=rstoch1D(1,1),                   &
                 ids=1,ide=1, jds=1,jde=1, kds=1,kde=kde,                         &
                 ims=1,ime=1, jms=1,jme=1, kms=1,kme=kme,                         &
                 its=1,ite=1, jts=1,jte=1, kts=1,kte=kte,                         &
                 errmsg=errmsg, errflg=errflg)
        end if

        if (wet(1,1)) then
            CALL mynnsfc_driver(                                                  &
                 u3d=u3d,v3d=v3d,t3d=t3d,qv3d=qv3d,p3d=p3d,dz8w=dz8w,             &
                 th3d=th3d,rho3d=rho3d,                                           &
                 PSFCPA=PSFCPA,PBLH=pblh,MAVAIL=mavail,XLAND=xland,DX=dx,         &
                 ISFFLX=isfflx,sf_mynn_sfcflux_water=isftcflx,flag_lsm=lsm,       &
                 sf_mynn_sfcflux_land=iz0tlnd,                                    &
                 sigmaf=sigmaf,vegtype=vegtype,shdmax=shdmax,ivegsrc=ivegsrc,     &  !intent(in)
                 z0pert=z0pert,ztpert=ztpert,                                     &  !intent(in)
                 redrag=redrag,sfc_z0_type=sfc_z0_type,                           &  !intent(in)
                 itimestep=itimestep,flag_iter=flag_iter,                         &
                 restart=flag_restart,dry=dry(1,1),wet=wet(1,1),icy=icy(1,1),     &  !intent(in)
                 tsk=tskin_wat, tsurf=tsurf_wat, qsfc=qsfc_wat, snowh=snowh_wat,  &
                 znt=znt_wat, ust=ust_wat,cm=cm_wat,ch=ch_wat, br=rb_wat,         &  !intent(inout)
                 stress=stress_wat,fm=psix_wat, fh=psit_wat,fm10=psix10_wat,      &  !intent(inout)
                 fh2=psit2_wat, hflx=hflx_wat, qflx=qflx_wat, CHS=chs,CHS2=chs2,  &
                 CQS2=cqs2,CQS=cqs,CPM=cpm,USTM=ustm,ZOL=zol,MOL=mol,RMOL=rmol,   &
                 psim=psim,psih=psih,HFX=hfx,QFX=qfx,LH=lh,FLHC=flhc,FLQC=flqc,   &
                 QGH=qgh,U10=u10,V10=v10,TH2=th2,T2=t2,Q2=q2,                     &
                 GZ1OZ0=GZ1OZ0,WSPD=wspd,                                         &
                 spp_pbl=spp_sfc,pattern_spp_pbl=rstoch1D(1,1),                   &
                 ids=1,ide=1, jds=1,jde=1, kds=1,kde=kde,                         &
                 ims=1,ime=1, jms=1,jme=1, kms=1,kme=kme,                         &
                 its=1,ite=1, jts=1,jte=1, kts=1,kte=kte,                         &
                 errmsg=errmsg, errflg=errflg)
        end if

        if (icy(1,1)) then
            CALL mynnsfc_driver(                                                  &
                 u3d=u3d,v3d=v3d,t3d=t3d,qv3d=qv3d,p3d=p3d,dz8w=dz8w,             &
                 th3d=th3d,rho3d=rho3d,                                           &
                 PSFCPA=PSFCPA,PBLH=pblh,MAVAIL=mavail,XLAND=xland,DX=dx,         &
                 ISFFLX=isfflx,sf_mynn_sfcflux_water=isftcflx,flag_lsm=lsm,       &
                 sf_mynn_sfcflux_land=iz0tlnd,                                    &
                 sigmaf=sigmaf,vegtype=vegtype,shdmax=shdmax,ivegsrc=ivegsrc,     &  !intent(in)
                 z0pert=z0pert,ztpert=ztpert,                                     &  !intent(in)
                 redrag=redrag,sfc_z0_type=sfc_z0_type,                           &  !intent(in)
                 itimestep=itimestep,flag_iter=flag_iter,                         &
                 restart=flag_restart,dry=dry(1,1),wet=wet(1,1),icy=icy(1,1),     &  !intent(in)
                 tsk=tskin_ice, tsurf=tsurf_ice, qsfc=qsfc_ice, snowh=snowh_ice,  &
                 znt=znt_ice, ust=ust_ice,cm=cm_ice,ch=ch_ice, br=rb_ice,         &  !intent(inout)
                 stress=stress_ice,fm=psix_ice, fh=psit_ice,fm10=psix10_ice,      &  !intent(inout)
                 fh2=psit2_ice, hflx=hflx_ice, qflx=qflx_ice, CHS=chs,CHS2=chs2,  &
                 CQS2=cqs2,CQS=cqs,CPM=cpm,USTM=ustm,ZOL=zol,MOL=mol,RMOL=rmol,   &
                 psim=psim,psih=psih,HFX=hfx,QFX=qfx,LH=lh,FLHC=flhc,FLQC=flqc,   &
                 QGH=qgh,U10=u10,V10=v10,TH2=th2,T2=t2,Q2=q2,                     &
                 GZ1OZ0=GZ1OZ0,WSPD=wspd,                                         &
                 spp_pbl=spp_sfc,pattern_spp_pbl=rstoch1D(1,1),                   &
                 ids=1,ide=1, jds=1,jde=1, kds=1,kde=kde,                         &
                 ims=1,ime=1, jms=1,jme=1, kms=1,kme=kme,                         &
                 its=1,ite=1, jts=1,jte=1, kts=1,kte=kte,                         &
                 errmsg=errmsg, errflg=errflg)
        end if

         write(0,*) "T2=",T2(1,1)
         write(0,*) "Read status:", read_stat
         if (saveoutput) then
             open(output_unit, file = './data/ccpp_output_lnd.txt')
             write(output_unit,'(I5, I5, F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F15.2,F10.2)') itimestep, &
                  iter,T2(1,1),Q2(1,1),TH2(1,1),U10(1,1),V10(1,1),HFX(1,1),LH(1,1),UST_lnd(1,1),PBLH(1,1)
         end if

        end do

    end subroutine ccpp_test
    !=================================================================================================================

end module module_sf_mynnsfc_ccpp_tests           
