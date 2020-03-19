module GEMSTOOL_WriteResults_m

contains

subroutine GEMSTOOL_WriteResults


   implicit none

!  Section 1. OPEN FILES of RESULTS, write Headers
!  ***********************************************

!  revision 09/03/12; Nstokes = 3/4 possible, V3 or V4 files.

!  REGULAR OUTPUT Filenames, either TOA-up or BOA-dn or both. Vector or scalar

  write(Char1,'(I1)')nstokes
  if ( nstokes .eq. 1 ) then
    filename1(1)  = 'results/regS'//CHAR1//'_TOAUP.out'
    filename1(2)  = 'results/regS'//CHAR1//'_BOADN.out'
  else
    filename1(1)  = 'results/regV'//CHAR1//'_TOAUP.out'
    filename1(2)  = 'results/regV'//CHAR1//'_BOADN.out'
  endif

!  FLUX OUTPUT Filenames, either TOA-up or BOA-dn or both. Vector or scalar

  if ( do_additional_mvout .or. do_mvout_only ) then
     if ( nstokes .eq. 1 ) then
        filename1_mf(1)  = 'results/regS'//CHAR1//'_TOAUP.mfout'
        filename1_mf(2)  = 'results/regS'//CHAR1//'_BOADN.mfout'
     else
        filename1_mf(1)  = 'results/regV'//CHAR1//'_TOAUP.mfout'
        filename1_mf(2)  = 'results/regV'//CHAR1//'_BOADN.mfout'
     endif
  endif

!  Cloudy/Clear

   cloud_annex = '_Clear_M0'
   if ( cloud_method .eq. 1 ) cloud_annex = '_Clear_M1'
   if ( cloud_method .eq. 2 ) cloud_annex = '_Clear_M2'
   if ( do_clouds )  cloud_annex(2:6) = 'Cloud'

   do i = 1, nd
      filename1(i)  = filename1(i) (1:LEN(TRIM(filename1(i))))//cloud_annex
      if ( do_additional_mvout .or. do_mvout_only ) then
         filename1_mf(i)  = filename1_mf(i) (1:LEN(TRIM(filename1_mf(i))))//cloud_annex
      endif
   enddo

!  open files, using the direction mask. Regular Output

  DO i = 1, nd
     dm = dirmask(i)
     open(98+I,file=filename1(dm), status = 'replace')
     if ( nstokes .eq. 1 ) then
        write(98+I, '(a,a/)') '       Wav       SZA      VZA      RAA      ',' I'
     else if ( nstokes .eq. 3 ) then
        write(98+I, '(a,a/)') '       Wav       SZA      VZA      RAA      ',&
                                     ' I              Q             U            DoLP '
     else if ( nstokes .eq. 4 ) then
        write(98+I, '(a,a/)') '       Wav       SZA      VZA      RAA      ',&
                                     ' I              Q             U             V            DoLP            DoCP'
     endif
  enddo

!  open files, using the direction mask. Flux Output
!    Downwelling has additional Direct-Flux entries

  if ( do_additional_mvout .or. do_mvout_only ) then
     DO i = 1, nd
        dm = dirmask(i)
        open(78+I,file=filename1_mf(dm), status = 'replace')
        if ( dm.eq.2 ) then
           if ( nstokes .eq. 1 ) then
              write(78+I, '(a,a,a/)') '       Wav       SZA',       &
                                    '     ADirect I     Actinic I', &
                                    '     FDirect I     Regflux I'
           else if ( nstokes .eq. 3 ) then
              write(78+I, '(a,a,a/)') '       Wav       SZA',                                   &
                                    '     ADirect I     Actinic I     Actinic Q     Actinic U', &
                                    '     FDirect I     Regflux I     Regflux Q     Regflux U'
           else if ( nstokes .eq. 4 ) then
              write(78+I, '(a,a,a/)') '       Wav       SZA',                                                 &
                                    '     ADirect I     Actinic I     Actinic Q     Actinic U     Actinic V', &
                                    '     FDirect I     Regflux I     Regflux Q     Regflux U     Regflux V'
           endif
        else
           if ( nstokes .eq. 1 ) then
              write(78+I, '(a,a,a/)') '       Wav       SZA     Actinic I     Regflux I'
           else if ( nstokes .eq. 3 ) then
              write(78+I, '(a,a,a/)') '       Wav       SZA',                     &
                                    '     Actinic I     Actinic Q     Actinic U', &
                                    '     Regflux I     Regflux Q     Regflux U'
           else if ( nstokes .eq. 4 ) then
              write(78+I, '(a,a,a/)') '       Wav       SZA',                                   &
                                    '     Actinic I     Actinic Q     Actinic U     Actinic V', &
                                    '     Regflux I     Regflux Q     Regflux U     Regflux V'
           endif
        endif
     enddo
  endif

!  contour output

  if ( do_contour_output ) then

     if (do_ss_external.and.FO_do_regular_ps)unit1        = 73
     if (do_ss_external.and..not.FO_do_regular_ps)unit1   = 74
     if (.not.do_ss_external.and.VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR)unit1    = 75
     if (.not.do_ss_external.and.VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING)unit1 = 76

     if (do_ss_external.and.FO_do_regular_ps)unit2        = 93
     if (do_ss_external.and..not.FO_do_regular_ps)unit2   = 94
     if (.not.do_ss_external.and.VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR)unit2    = 95
     if (.not.do_ss_external.and.VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING)unit2 = 96
     write(cf1,'(I2)')unit1
     write(cf2,'(I2)')unit2

     open(unit1,file='results/SingleScatterRadiances_'//cf1//'.dat_FD',status='replace')
     open(unit2,file='results/MultipleScatterRadiances_'//cf2//'.dat_FD',status='replace')

  endif


   return
end subroutine GEMSTOOL_WriteResults


end module GEMSTOOL_WriteResults_m
