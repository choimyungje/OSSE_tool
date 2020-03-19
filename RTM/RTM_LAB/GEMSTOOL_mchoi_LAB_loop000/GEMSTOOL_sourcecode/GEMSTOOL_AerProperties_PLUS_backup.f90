Module GEMSTOOL_AerProperties_Plus_m

!  Rob Fix, 10/20/16.
!  New facility for reading and Writing User-defined optical properties

!  Use Type structures (GEMSTOOL)

   use GEMSTOOL_Input_Types_m

!  Auxiliary files from the GEMSTOOL_Read_Configfiles_m

   use GEMSTOOL_Read_Configfiles_m, only : Get_FileUnit, &
                                           Read_aerosol_oprop_file , &
                                           Read_aerosol_oprop_file_plus, &
                                           Read_aerosol_oprop_file_plus_short, &
                                           Write_aerosol_oprop_file, &
                                           Write_aerosol_oprop_file_plus

!  Use Loading routines
!  --------------------

   use GEMSTOOL_aerload_routines_m

!  Use Mie code
!  ------------

!  parameters, sourcecode (single-mode masters)

   use RTSMie_parameters_m
   use RTSMie_sourcecode_m
   use RTSMie_sourcecode_plus_m

!  Use Tmatrix code
!  ----------------

!  parameters, single-mode masters

   use tmat_parameters
   use tmat_master_m
   use tmat_master_PLUS_m

!  All routines are public

public  :: GEMSTOOL_AERWFS_BOOKKEEP, GEMSTOOL_AER_PROPERTIES_PLUS 
private :: fpk

contains

subroutine GEMSTOOL_AERWFS_BOOKKEEP &
      ( maxaerwfs, GEMSTOOL_INPUTS,                        & ! Inputs & InOut
        n_AerBulk_Total_wfs, n_AerBulk_Group_wfs,          & ! OUTPUT, aerosol Jacobian Bookkeeping
        AerBulk_Mask_wfs, AerBulk_pars,                    & ! OUTPUT, aerosol Jacobian Bookkeeping
        do_LinearRef, do_LinearPSD, Do_LinearEps,          & ! OUTPUT, aerosol Jacobian Bookkeeping
        NLIN, NPSD, Parslist,                              & ! OUTPUT, aerosol Jacobian Bookkeeping
        FAIL2, MTU_Bookkeep )                                ! Exception handling

!  =============================================================================
!                        AEROSOL JACOBIAN CREATION
!  =============================================================================

!  Aerosol Loading:
!    Loading Jacobians
!       Cases 1-3 : dloading_Dtau      = derivative of profile w.r.t the total aerosol optical depth at wavelength w0
!       Case 2    : dloading_Dpars(1)  = derivative of profile w.r.t the relaxation parameter (Exponential)
!       Case 3    : dloading_Dpars(1)  = derivative of profile w.r.t the GDF Peak height
!       Case 3    : dloading_Dpars(2)  = derivative of profile w.r.t the GDF Half Width 

!  Generation of Aerosol optical properties
!     1. Call to the Linearized Mie/Tmatrix program
!     2. Convert Mie/Tmatrix output (Microsopic) to IOP output (macroscopic) 

!  Here is a breakdown of the aerosol weighting functions
!  ------------------------------------------------------

!  There are 5 groups, according to the input control

!  Group 1.  do_AerBulk_LoadPars_Jacobians --> w.r.t. Profile loading parameters
!           # 1  Total optical depth of aerosol (aerau_input_w0) at reference wavelength w0
!           # 2  EITHER Relaxation parameter in EXP loading; OR peakheight parameter in GDF loading
!           # 3  halfwidth parameter in GDF Loading

!  Group 2.  do_AerBulk_RefIndex_Jacobians --> w.r.t. Refractive Index components
!           # 1  Real part of the refractive index for Mode-1
!           # 2  Imag part of the refractive index for Mode-1
!           # 3  Real part of the refractive index for Mode-2 (if present)
!           # 4  Imag part of the refractive index for Mode-2 (if present)

!  Group 3.  do_AerBulk_ShapeFac_Jacobians --> w.r.t. Shape Factor (T-matrix only)
!           # 1  Shape Factor for Mode-1
!           # 2  Shape Factor for Mode-2 (if present)

!  Group 4.  do_AerBulk_SizeDist_Jacobians --> w.r.t. Size Distribution Parameters
!           # 1  LN Mode radius  of Particle Size Distribution for Mode-1
!           # 2  LN Std. dev.    of Particle Size Distribution for Mode-1
!           # 3  3rd parameter   of Particle Size Distribution for Mode-1 (if present)
!           # 4  LN Mode radius  of Particle Size Distribution for Mode-2 (if present)
!           # 5  LN Std. dev.    of Particle Size Distribution for Mode-2 (if present)
!           # 6  3rd parameter   of Particle Size Distribution for Mode-2 (if present)

!  Group 5.  do_AerBulk_BmodFrac_Jacobian --> w.r.t.Bimodal fraction
!           # 1  Bimodal fraction for PSDs                                (if present)

!  Total possible number of Jacobians = 16 in all 5 Groups

!  ORDER OF APPEARANCE for A BIMODAL TMATRIX with both PSDs having 3 parameters, with GDF LOADING
!  ----------------------------------------------------------------------------------------------

!   WF #   GROUP  PSD MODE #   Jacobian parameter

!     1      2       1         Real part of the refractive index for Mode-1
!     2      2       1         Imag part of the refractive index for Mode-1
!     3      3       1         Shape Factor for Mode-1
!     4      4       1         Parameter 1 of Particle Size Distribution for Mode-1
!     5      4       1         Parameter 2 of Particle Size Distribution for Mode-1
!     6      4       1         Parameter 3 of Particle Size Distribution for Mode-1
!     7      2       2         Real part of the refractive index for Mode-2
!     8      2       2         Imag part of the refractive index for Mode-2
!     9      3       2         Shape Factor for Mode-2
!     10     4       2         Parameter 1 of Particle Size Distribution for Mode-2
!     11     4       2         Parameter 2 of Particle Size Distribution for Mode-2
!     12     4       2         Parameter 3 of Particle Size Distribution for Mode-2
!     13     5       -         Bimodal fraction for PSDs
!     14     1       -         Total optical depth of aerosol (aerau_input_w0) at reference wavelength w0
!     15     1       -         Peak Height parameter in GDF loading
!     16     1       -         Half Width  parameter in GDF Loading

!  ORDER OF APPEARANCE for A BIMODAL MIE Aerosol with both PSDs having 2 parameters, with GDF LOADING
!  --------------------------------------------------------------------------------------------------

!   WF #   GROUP  PSD MODE #   Jacobian parameter

!     1      2       1         Real part of the refractive index for Mode-1
!     2      2       1         Imag part of the refractive index for Mode-1
!     3      4       1         Parameter 1 of Particle Size Distribution for Mode-1 e.g. Log-Normal Mode Radius
!     4      4       1         Parameter 2 of Particle Size Distribution for Mode-1 e.g. Log-Normal Standard Deviation
!     5      2       2         Real part of the refractive index for Mode-2
!     6      2       2         Imag part of the refractive index for Mode-2
!     7      4       2         Parameter 1 of Particle Size Distribution for Mode-2 e.g. Log-Normal Mode Radius
!     8      4       2         Parameter 2 of Particle Size Distribution for Mode-2 e.g. Log-Normal Standard Deviation
!     9      5       -         Bimodal fraction for PSDs
!     10     1       -         Total optical depth of aerosol (aerau_input_w0) at reference wavelength w0
!     11     1       -         Peak Height parameter in GDF loading
!     12     1       -         Half Width  parameter in GDF Loading

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)
  
!  Input variables
!  ---------------

!  External Dimensioning

   integer, INTENT (IN) :: maxaerwfs

!  Type Structure inputs. May set some flags here

   TYPE(GEMSTOOL_Config_Inputs), intent(INOUT) :: GEMSTOOL_INPUTS

!  output
!  ======

!  General output for  Bookkeeping
!    Total number of aerosol weighting functions = n_AerBulk_Total_wfs
!    Group number of aerosol weighting functions = n_AerBulk_Group_wfs(5)
!    Total collection of aerosol parameters      = AerBulk_pars
!    Mode lists of aerosol parameters            = Parslist

   INTEGER     , INTENT (OUT)   :: n_AerBulk_Total_wfs, n_AerBulk_Group_wfs(5), AerBulk_Mask_wfs(maxaerwfs)
   REAL(fpk)   , INTENT (OUT)   :: AerBulk_pars(maxaerwfs)
!mick mod 1/24/2017 - moved Do_LinearRef, Do_LinearPSD, Do_LinearEps, nlin(2), & npsd(2) from internal variables
!                     to output variables that they may be defined here in "GEMSTOOL_AERWFS_BOOKKEEP" whether they 
!                     be defined by Mie or Tmatrix parameters or User file data 
   LOGICAL     , INTENT (OUT)   :: Do_LinearRef, Do_LinearPSD, Do_LinearEps
   INTEGER     , INTENT (OUT)   :: nlin(2), npsd(2)
   REAL(fpk)   , INTENT (OUT)   :: Parslist(10,2)

!  Exception handling

   LOGICAL      , intent(out) :: FAIL2
   CHARACTER*(*), intent(out) :: MTU_Bookkeep(3)

!  Local variables
!  ===============

!  Linearization w.r.t refractive index components
!  Linearization w.r.t shape factor (T-matrix only)
!  Linearization for PSD parameters
!  Linearization for Loading parameters
!  Linearization for Bimodal Fraction
!  Overall Linearization

   !logical    :: Do_LinearRef
   !logical    :: Do_LinearEps
   !logical    :: Do_LinearPSD
   logical    :: do_LinearLoad
   logical    :: do_LinearBmf
   logical    :: do_aer_Jacobians

!  Local mask

   integer :: distpars_mask(9)
   data distpars_mask / 2,2,3,2,3,3,3,3,0 /

!  Help Variables

   logical      :: do_bimodal, RefWvl_dummy
   integer      :: m, m1, m2, n1, qm1, qm2, qmt, kd, n_Aer_wfs_1, n_Aer_wfs_2 
   integer      :: nm, nmodes, istatus
   !integer      :: nlin(2), npsd(2)
   real(fpk)    :: Wavelength_dummy
   character*90 :: AeroInFile_Plus
   integer      :: InFileUnit(2)

!mick mod 1/24/2017 - added aerosol file preamble capacity
!  Rob mod 3/7/17 - expanded preamble to include ref wavelength and PSD information (11 lines in all)
!  Aerosol output file description header
!   integer, parameter :: nPLines = 6
   integer, parameter :: nPLines = 11
   logical            :: DoAeroFilePreamble(1)

!  Initialize Aerosol Jacobian Bookkeeping

   n_AerBulk_Total_Wfs = 0
   n_AerBulk_Group_Wfs = 0
   AerBulk_Mask_Wfs    = 0
   AerBulk_pars        = 0.0d0
   Parslist            = 0.0d0
!mick fix 1/24/2017 - initialize
   fail2 = .false.

!  First wavelength is always the reference

   Wavelength_dummy    = 0.001d0 * GEMSTOOL_INPUTS%AerLoad%reference_w0
   RefWvl_dummy        = .true.

!@@@@@@@@@@@@@@@@@@@@@@@@ START JACOBIAN BOOKKEEPING @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Proxy bimodal flag, number of modes

   do_bimodal = GEMSTOOL_INPUTS%MieTmatUser%do_bimodal
   nmodes = 1 ; if ( do_bimodal ) nmodes = 2
   do_LinearBmf  = do_bimodal .and. GEMSTOOL_INPUTS%LinControl%do_AerBulk_BmodFrac_Jacobian

!  Proxies for the Loading

   do_LinearLoad = GEMSTOOL_INPUTS%LinControl%do_AerBulk_LoadPars_Jacobians  !!CHECK!!Y.Jung

!  Set variables using User-defined aerosols, as alternative
!  =========================================================

   if ( GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then

!  For each mode, Open file and read headers using the "Short" routine

      do nm = 1, nmodes
         call Get_FileUnit(InFileUnit(nm))
         AeroInFile_Plus='UserAeroFiles/in/' // trim(GEMSTOOL_INPUTS%MieTmatUser%User_Aerofile_names(nm))//'_PLUS'
write(*,*) AeroInFile_Plus
!pause
         open(unit=InFileUnit(nm), file=trim(AeroInFile_Plus),iostat=istatus,status='old')
         if ( istatus > 0 ) then
            MTU_Bookkeep(1) = 'Failure : Aerosol optical property input file, short version'
            MTU_Bookkeep(2) = 'File    : ' // trim(AeroInFile_Plus)
            MTU_Bookkeep(3) = 'Problem : Error Opening file (not found)'
            fail2 = .true. ; close(unit=InFileUnit(nm)) ; return
         endif

!mick mod 1/24/2017 - added aerosol file preamble capacity
         DoAeroFilePreamble(1) = .true.
         Call Read_aerosol_oprop_file_plus_short ( &
             InFileUnit(nm), Wavelength_dummy, RefWvl_dummy,                 & !Input
             1, 1, nPLines, DoAeroFilePreamble,                              & !I/O
             do_LinearRef, do_LinearPSD, NLIN(nm), NPSD(nm), ParsList(1,nm), & !Output
             fail2, MTU_Bookkeep )                                             !Output
         if ( fail2 ) then
            MTU_Bookkeep(2) = 'File   : ' // trim(AeroInFile_Plus)
            close(unit=InFileUnit(nm)) ; return
         endif
         close(unit=InFileUnit(nm))
      enddo
!mick mod 1/24/2017 - defined default value of "do_LinearEps" for User case
      do_LinearEps = .false.

!  Set overall Jacobian flag, return if not true

      do_aer_Jacobians = do_LinearBmf .or. do_LinearLoad .or. do_LinearRef .or. do_LinearPSD 
      if ( .not.do_aer_Jacobians ) return

!  Set up aerosol parameters and total number of aerosol WFs, First aerosol mode

      m = 0 ; nm = 1
      if ( do_LinearRef ) then
         m1 = m + 1; m2 = m + NLIN(nm)
         AerBulk_pars(m1:m2) = ParsList(1:NLIN(nm),nm)
         if ( nlin(nm).eq.3 ) then
            n1 =  NLIN(nm) - 1 ; do_LinearEps = .true.
            AerBulk_Mask_wfs(1:N1)     = 2 ; n_AerBulk_Group_wfs(2) = n_AerBulk_Group_wfs(2) + N1 ; m = m + N1
            AerBulk_Mask_wfs(NLIN(nm)) = 3 ; n_AerBulk_Group_wfs(3) = n_AerBulk_Group_wfs(3) + 1  ; m = m + 1
         else
            n1 =  NLIN(nm)
            AerBulk_Mask_wfs(1:N1)     = 2 ; n_AerBulk_Group_wfs(2) = n_AerBulk_Group_wfs(2) + N1 ; m = m + N1
         endif
      endif
      if ( do_LinearPSD ) then
         n1 = npsd(nm)
         m1 = m + 1; m2 = m + n1
         AerBulk_pars    (m1:m2) = ParsList(NLIN(nm)+1:NLIN(nm)+N1,nm)
         AerBulk_Mask_wfs(m1:m2) = 4 ; n_AerBulk_Group_wfs(4) = n_AerBulk_Group_wfs(4) + N1; m = m + N1
      endif
      n_Aer_wfs_1 = m

!  Set up aerosol parameters and total number of aerosol WFs, Second aerosol mode

      if ( nmodes.eq.2 ) then
         nm = 2
         if ( do_LinearRef ) then
            m1 = m + 1; m2 = m + NLIN(nm)
            AerBulk_pars(m1:m2) = ParsList(1:NLIN(nm),nm)
            if ( nlin(nm).eq.3 ) then
               n1 =  NLIN(nm) - 1 ; do_LinearEps = .true.
               AerBulk_Mask_wfs(m+1:m+N1)   = 2 ; n_AerBulk_Group_wfs(2) = n_AerBulk_Group_wfs(2) + N1 ; m = m + N1
               AerBulk_Mask_wfs(m+NLIN(nm)) = 3 ; n_AerBulk_Group_wfs(3) = n_AerBulk_Group_wfs(3) + 1  ; m = m + 1
            else
               n1 =  NLIN(nm)
               AerBulk_Mask_wfs(m+1:m+N1)      = 2 ; n_AerBulk_Group_wfs(2) = n_AerBulk_Group_wfs(2) + N1 ; m = m + N1
            endif
         endif
         if ( do_LinearPSD ) then
            n1 = npsd(nm)
            m1 = m + 1; m2 = m + n1
            AerBulk_pars    (m1:m2) = ParsList(NLIN(nm)+1:NLIN(nm)+N1,nm)
            AerBulk_Mask_wfs(m1:m2) = 4 ; n_AerBulk_Group_wfs(4) = n_AerBulk_Group_wfs(4) + N1; m = m + N1
         endif
         n_Aer_wfs_2 = m - n_Aer_wfs_1
      endif

!  Set up bookkeeping for Bimodal fraction

      if ( nmodes.eq.2 ) then
         do_LinearBmf  = .true.
         AerBulk_pars(m+1)     = GEMSTOOL_INPUTS%MieTmatUser%bimodal_fraction
         AerBulk_Mask_wfs(m+1) = 5 ; n_AerBulk_Group_wfs(5) = n_AerBulk_Group_wfs(5) + 1
      endif 

!  done

   endif

!  Set variables from Configuration file inputs (Mie or Tmat, not the USER-defined aerosols)
!  ============================================

   if ( .not.GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then

!  Set local Mie/Tmatrix Linearization variables, for the 5 input GROUPS

      do_LinearBmf  = do_bimodal .and. GEMSTOOL_INPUTS%LinControl%do_AerBulk_BmodFrac_Jacobian
      do_LinearLoad = GEMSTOOL_INPUTS%LinControl%do_AerBulk_LoadPars_Jacobians  !!CHECK!!Y.Jung

      do_LinearRef  = GEMSTOOL_INPUTS%LinControl%do_AerBulk_RefIndex_Jacobians
      do_LinearPSD  = GEMSTOOL_INPUTS%LinControl%do_AerBulk_SizeDist_Jacobians

      do_LinearEps  = GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols .and. &
                      GEMSTOOL_INPUTS%LinControl%do_AerBulk_ShapeFac_Jacobians

!  Overall Jacobian flag, return if not true

      do_aer_Jacobians = do_LinearBmf .or. do_LinearLoad .or. do_LinearRef .or. do_LinearPSD .or. do_LinearEps
      if ( .not. do_aer_Jacobians ) return

!  Set up aerosol parameters and total number of aerosol WFs, First aerosol mode
!mick mod 1/24/2017 - added defining of nlin(1) here

      nlin(1) = 2 ; if (do_LinearEps) nlin(1) = nlin(1) + 1
      npsd(1) = distpars_mask(GEMSTOOL_INPUTS%MieTmatUser%PSDIndex(1))
      m = 0
      if ( do_LinearRef ) then
         AerBulk_pars(1) = GEMSTOOL_INPUTS%MieTmatUser%nreal(1) ; AerBulk_pars(2) = GEMSTOOL_INPUTS%MieTmatUser%nimag(1)
         AerBulk_Mask_wfs(1:2) = 2 ;  n_AerBulk_Group_wfs(2) = n_AerBulk_Group_wfs(2) + 2 ; m = m + 2
      endif
      if ( do_LinearEps ) then
         AerBulk_pars(3) = GEMSTOOL_INPUTS%MieTmatUser%Tmat_eps(1)
         AerBulk_Mask_wfs(3) = 3 ; n_AerBulk_Group_wfs(3) = n_AerBulk_Group_wfs(3) + 1 ; m = m + 1
      endif
      if ( do_LinearPSD ) then
         do kd = 1, npsd(1)
            m = m + 1 ; AerBulk_pars(m) = GEMSTOOL_INPUTS%MieTmatUser%PSDpars(kd,1)
            AerBulk_Mask_wfs(m) = 4 ; n_AerBulk_Group_wfs(4) = n_AerBulk_Group_wfs(4) + 1 
         end do
      endif
      n_Aer_wfs_1 = m ; Parslist(1:n_Aer_wfs_1,1) =  AerBulk_pars(1:m)

!  Set up aerosol parameters and total number of aerosol WFs, Second aerosol mode
!mick mod 1/24/2017 - added defining of nlin(2) here

      if ( do_bimodal ) then
         nlin(2) = 2 ; if (do_LinearEps) nlin(2) = nlin(2) + 1
         npsd(2) = distpars_mask(GEMSTOOL_INPUTS%MieTmatUser%PSDIndex(2))
         m = n_Aer_wfs_1
         if ( do_LinearRef ) then
            AerBulk_pars(m+1) = GEMSTOOL_INPUTS%MieTmatUser%nreal(2) ; AerBulk_pars(m+2) = GEMSTOOL_INPUTS%MieTmatUser%nimag(2)
            AerBulk_Mask_wfs(m+1:m+2) = 2 ;  n_AerBulk_Group_wfs(2) = n_AerBulk_Group_wfs(2) + 2  ; m = m + 2  
         endif
         if ( do_LinearEps ) then
            AerBulk_pars(m+1) = GEMSTOOL_INPUTS%MieTmatUser%Tmat_eps(2)
            AerBulk_Mask_wfs(m+1) = 3 ; n_AerBulk_Group_wfs(3) = n_AerBulk_Group_wfs(3) + 1  ; m = m + 1
         endif
         if ( do_LinearPSD ) then
            do kd = 1, npsd(2)
               m = m + 1 ; AerBulk_pars(m) = GEMSTOOL_INPUTS%MieTmatUser%PSDpars(kd,2)
               AerBulk_Mask_wfs(m) = 4 ; n_AerBulk_Group_wfs(4) = n_AerBulk_Group_wfs(4) + 1 
            end do
         endif
         n_Aer_wfs_2 = m - n_Aer_wfs_1
         Parslist(1:n_Aer_wfs_2,2) =  AerBulk_pars(n_Aer_wfs_1+1:n_Aer_wfs_1+n_Aer_wfs_2)
         if ( do_LinearBmf ) then
            AerBulk_pars(m+1)     = GEMSTOOL_INPUTS%MieTmatUser%bimodal_fraction
            AerBulk_Mask_wfs(m+1) = 5 ; n_AerBulk_Group_wfs(5) = n_AerBulk_Group_wfs(5) + 1
         endif 
      endif

!  Done

   endif

!  Total number of WFs

    qm1 = n_Aer_wfs_1
    if ( do_bimodal )                   qm1 = qm1 + n_Aer_wfs_2
    if ( do_bimodal.and.do_LinearBmf )  qm1 = qm1 + 1

!  Add WFs for the loading
!    Total loading + relaxation for the Exponential profile
!    Total loading + peak height + half width for the GDF profile
!    Total loading only for the Uniform profile.

   if ( do_LinearLoad ) then
      if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 2) then
         qmt = qm1 + 2
      else if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 3 ) then
         qmt = qm1 + 3
      else
         qmt = qm1 + 1
      endif
   else
      qmt = qm1
   endif

!  Add 1 or 2-3 more parameters for loading 

   qm2 = qm1 + 1
   AerBulk_pars(qm2) = GEMSTOOL_INPUTS%AerLoad%aertau_input_w0
   AerBulk_Mask_wfs(qm2) = 1 ; n_AerBulk_Group_wfs(1) = 1 
   if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 2 ) then
      AerBulk_pars(qmt) = GEMSTOOL_INPUTS%AerLoad%exploading_relaxation
      AerBulk_Mask_wfs(qmt) = 1   ; n_AerBulk_Group_wfs(1) = n_AerBulk_Group_wfs(1) + 1 
   else if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 3 ) then
      AerBulk_pars(qm2+1) = GEMSTOOL_INPUTS%AerLoad%gdfloading_peakheight
      AerBulk_pars(qmt)   = GEMSTOOL_INPUTS%AerLoad%gdfloading_halfwidth
      AerBulk_Mask_wfs(qm2+1) = 1  ; AerBulk_Mask_wfs(qmt) = 1 ; n_AerBulk_Group_wfs(1) = n_AerBulk_Group_wfs(1) + 2
   endif

!  total number of weighing functions

   n_AerBulk_Total_Wfs = qmt

!  debug
!   write(*,*)'PARSLIST',parslist(1:5,1)
!   write(*,*)'PARSLIST',parslist(1:5,2)
!   write(*,*)qmt,do_bimodal, do_LinearBmf
!   pause

!@@@@@@@@@@@@@@@@@@@@@@@@ END JACOBIAN BOOKKEEPING @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   return
end subroutine GEMSTOOL_AERWFS_BOOKKEEP

!

subroutine GEMSTOOL_AER_PROPERTIES_PLUS &
      ( maxlayers, maxwav, maxaermoms, maxaerwfs, interpolate_aerosols, do_wavnums, & ! Dimensions
        nlayers, nmuller, nwav, wav, height_grid, GEMSTOOL_INPUTS, momsize_cutoff,  & ! Inputs
        n_AerBulk_Total_wfs, n_AerBulk_Group_wfs, AerBulk_pars,                     & ! Inputs, aerosol Jacobian Bookkeeping
        do_LinearRef, do_LinearPSD, Do_LinearEps, NLIN, NPSD, ParsList,             & ! Inputs, aerosol Jacobian Bookkeeping
        aerlayerflags, Loading, Dloading_Dtau, Dloading_Dpars,                      & ! OUTPUT, aerosol Loading
        n_scatmoms, aertau_unscaled, aod_scaling, aerosol_distchars,                & ! OUTPUT, aerosol optical properties
        aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms,                           & ! OUTPUT, aerosol optical properties
        L_aertau_unscaled, L_aod_scaling, L_aerosol_deltau,                         & ! OUTPUT, aerosol Linearized OPs
        L_aerosol_ssalbs, L_aerosol_scatmoms,                                       & ! OUTPUT, aerosol Linearized OPs
        fail1, fail2, Message_Loading, Messages_Optical )                             ! Exception handling

!  =============================================================================
!                        AEROSOL JACOBIAN CREATION
!  =============================================================================

!  aerosol Loading:
!    Loading = optical depth profile
!    Loading Jacobians
!       Cases 1-3 : dloading_Dtau      = derivative of profile w.r.t the total aerosol optical depth at wavelength w0
!       Case 2    : dloading_Dpars(1)  = derivative of profile w.r.t the relaxation parameter (Exponential)
!       Case 3    : dloading_Dpars(1)  = derivative of profile w.r.t the GDF Peak height
!       Case 3    : dloading_Dpars(2)  = derivative of profile w.r.t the GDF Half Width 

!  Generation of Aerosol optical properties
!     1. Call to the Linearized Mie/Tmatrix program
!     2. Convert Mie/Tmatrix output (Microsopic) to IOP output (macroscopic) 

!  Here is a breakdown of the aerosol weighting functions
!  ------------------------------------------------------

!  There are 5 groups, according to the input control

!  Group 1.  do_AerBulk_LoadPars_Jacobians --> w.r.t. Profile loading parameters
!           # 1  Total optical depth of aerosol (aerau_input_w0) at reference wavelength w0
!           # 2  EITHER Relaxation parameter in EXP loading; OR peakheight parameter in GDF loading
!           # 3  halfwidth parameter in GDF Loading

!  Group 2.  do_AerBulk_RefIndex_Jacobians --> w.r.t. Refractive Index components
!           # 1  Real part of the refractive index for Mode-1
!           # 2  Imag part of the refractive index for Mode-1
!           # 3  Real part of the refractive index for Mode-2 (if present)
!           # 4  Imag part of the refractive index for Mode-2 (if present)

!  Group 3.  do_AerBulk_ShapeFac_Jacobians --> w.r.t. Shape Factor (T-matrix only)
!           # 1  Shape Factor for Mode-1
!           # 2  Shape Factor for Mode-2 (if present)

!  Group 4.  do_AerBulk_SizeDist_Jacobians --> w.r.t. Size Distribution Parameters
!           # 1  LN Mode radius  of Particle Size Distribution for Mode-1
!           # 2  LN Std. dev.    of Particle Size Distribution for Mode-1
!           # 3  3rd parameter   of Particle Size Distribution for Mode-1 (if present)
!           # 4  LN Mode radius  of Particle Size Distribution for Mode-2 (if present)
!           # 5  LN Std. dev.    of Particle Size Distribution for Mode-2 (if present)
!           # 6  3rd parameter   of Particle Size Distribution for Mode-2 (if present)

!  Group 5.  do_AerBulk_BmodFrac_Jacobian --> w.r.t.Bimodal fraction
!           # 1  Bimodal fraction for PSDs                                (if present)

!  Total possible number of Jacobians = 16 in all 4 Groups

!  ORDER OF APPEARANCE for A BIMODAL TMATRIX with both PSDs having 3 parameters, with GDF LOADING
!  ----------------------------------------------------------------------------------------------

!   WF #   GROUP  PSD MODE #   Jacobian parameter

!     1      2       1         Real part of the refractive index for Mode-1
!     2      2       1         Imag part of the refractive index for Mode-1
!     3      3       1         Shape Factor for Mode-1
!     4      4       1         Parameter 1 of Particle Size Distribution for Mode-1
!     5      4       1         Parameter 2 of Particle Size Distribution for Mode-1
!     6      4       1         Parameter 3 of Particle Size Distribution for Mode-1
!     7      2       2         Real part of the refractive index for Mode-2
!     8      2       2         Imag part of the refractive index for Mode-2
!     9      3       2         Shape Factor for Mode-2
!     10     4       2         Parameter 1 of Particle Size Distribution for Mode-2
!     11     4       2         Parameter 2 of Particle Size Distribution for Mode-2
!     12     4       2         Parameter 3 of Particle Size Distribution for Mode-2
!     13     5       -         Bimodal fraction for PSDs
!     14     1       -         Total optical depth of aerosol (aerau_input_w0) at reference wavelength w0
!     15     1       -         Peak Height parameter in GDF loading
!     16     1       -         Half Width  parameter in GDF Loading

!  ORDER OF APPEARANCE for A BIMODAL MIE Aerosol with both PSDs having 2 parameters, with GDF LOADING
!  --------------------------------------------------------------------------------------------------

!   WF #   GROUP  PSD MODE #   Jacobian parameter

!     1      2       1         Real part of the refractive index for Mode-1
!     2      2       1         Imag part of the refractive index for Mode-1
!     3      4       1         Parameter 1 of Particle Size Distribution for Mode-1 e.g. Log-Normal Mode Radius
!     4      4       1         Parameter 2 of Particle Size Distribution for Mode-1 e.g. Log-Normal Standard Deviation
!     5      2       2         Real part of the refractive index for Mode-2
!     6      2       2         Imag part of the refractive index for Mode-2
!     7      4       2         Parameter 1 of Particle Size Distribution for Mode-2 e.g. Log-Normal Mode Radius
!     8      4       2         Parameter 2 of Particle Size Distribution for Mode-2 e.g. Log-Normal Standard Deviation
!     9      5       -         Bimodal fraction for PSDs
!     10     1       -         Total optical depth of aerosol (aerau_input_w0) at reference wavelength w0
!     11     1       -         Peak Height parameter in GDF loading
!     12     1       -         Half Width  parameter in GDF Loading

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)
  
!  Input variables
!  ---------------

!  External Dimensioning

   integer, INTENT (IN) :: maxlayers, maxwav
   integer, INTENT (IN) :: maxaermoms, maxaerwfs

!  Flag for interpolation of aerosols

   logical, INTENT(in)  :: interpolate_aerosols

!  General Bookkeeping input
!    Total number of aerosol weighting functions   = n_AerBulk_Total_wfs
!    Group number of aerosol weighting functions   = n_AerBulk_Group_wfs(5)
!    Total collection of aerosol parameters        = AerBulk_pars
!    Mode list of aerosol parameters               = ParsList

   INTEGER         , INTENT (in)   :: n_AerBulk_Total_wfs, n_AerBulk_Group_wfs(5)
   REAL    (fpk)   , INTENT (in)   :: AerBulk_pars(maxaerwfs)
!mick mod 1/24/2017 - added Do_LinearRef, Do_LinearPSD, Do_LinearEps, nlin(2), & npsd(2) to 
!                     input variables that they all may be defined in "GEMSTOOL_AERWFS_BOOKKEEP" 
!                     whether they be defined by Mie or Tmatrix parameters or User file data 
   LOGICAL         , INTENT (in)   :: Do_LinearRef, Do_LinearEps
   LOGICAL         , INTENT (inout):: Do_LinearPSD
   INTEGER         , INTENT (in)   :: nlin(2), npsd(2)
   REAL    (fpk)   , INTENT (in)   :: ParsList(10,2)

!  Flag for using wavnumber output

   logical, INTENT(in)  :: do_wavnums

!  Numbers

   integer, INTENT (IN) ::  nlayers
   integer, INTENT (IN) ::  nwav

!  Heights and wavelengths

   REAL    (fpk)   , INTENT (IN)   :: wav(maxwav)
   REAL    (fpk)   , INTENT (IN)   :: height_grid(0:maxlayers)

!    nmuller = 1 (scalar code), = 6 (vector code)

   integer, INTENT (IN) ::  nmuller

!  Aerosol moment size cutoff. DEFAULT = 0.001

   REAL    (fpk),    INTENT (IN)   :: momsize_cutoff

!  Type Structure inputs

   TYPE(GEMSTOOL_Config_Inputs) :: GEMSTOOL_INPUTS

!  Mie/Tmatrix PSD inputs
!  ======================

!  PSD inputs (distribution index, PSD parameters)

!      psd_Index      - Index for particle size distribution of spheres
!      psd_pars       - Parameters characterizing PSD (up to 3 allowed)

!    PSD_index = 1 : TWO-PARAMETER GAMMA with alpha and b given
!    PSD_index = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given
!    PSD_index = 3 : BIMODAL GAMMA with equal mode weights
!    PSD_index = 4 : LOG-NORMAL with rg and sigma given
!    PSD_index = 5 : LOG-NORMAL with reff and veff given
!    PSD_index = 6 : POWER LAW
!    PSD_index = 7 : MODIFIED GAMMA with alpha, rc and gamma given
!    PSD_index = 8 : MODIFIED GAMMA with alpha, b and gamma given

!  FixR1R2 : If  set, Use Internal routine to calculate R1 and R2 (outputs)
!            If Not set, Use Input R1 and R2 for PSD limits.
!  R1, R2         - Minimum and Maximum radii (Microns)
!  N_REAL, N_IMAG - real and imaginary parts, refractive index (N-i.GE.0)

!  Mie-specific inputs
!  ===================

!  Limiting particle size value. Set to 10000.0 default.
!   If you exceed this, program will tell you to increase dimensioning.

!  R1R2_cutoff particle size for setting R1 and R2 internally

!  PSD quadrature control
!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.

!  Tmatrix-specific inputs
!  =======================

!  Logical flag for using equal-surface-area sepcification

!      NKMAX.LE.988 is such that NKMAX+2 is the                        
!           number of Gaussian quadrature points used in               
!           integrating over the size distribution for particles
!           MKMAX should be set to -1 for Monodisperse

!      NDGS - parameter controlling the number of division points      
!             in computing integrals over the particle surface.        
!             For compact particles, the recommended value is 2.       
!             For highly aspherical particles larger values (3, 4,...) 
!             may be necessary to obtain convergence.                  
!             The code does not check convergence over this parameter. 
!             Therefore, control comparisons of results obtained with  
!             different NDGS-values are recommended.

!      EPS (Shape_factor) and NP (Spheroid type) - specify the shape of the particles.                
!             For spheroids NP=-1 and EPS is the ratio of the          
!                 horizontal to rotational axes.  EPS is larger than   
!                 1 for oblate spheroids and smaller than 1 for       
!                 prolate spheroids.                                   
!             For cylinders NP=-2 and EPS is the ratio of the          
!                 diameter to the length.                              
!             For Chebyshev particles NP must be positive and 
!                 is the degree of the Chebyshev polynomial, while     
!                 EPS is the deformation parameter                     

!      Accuracy       - accuracy of the computations

!  output
!  ======

!  Loading output

    real(fpk)   , dimension ( maxlayers ),  INTENT (OUT)      :: Loading
    real(fpk)   , dimension ( maxlayers ),  INTENT (OUT)      :: DLoading_Dtau
    real(fpk)   , dimension ( maxlayers, 2 ),  INTENT (OUT)   :: Dloading_Dpars

!  aerosol layer flags

    LOGICAL,      DIMENSION ( maxlayers ), intent(out)  :: AERLAYERFLAGS 

!  AOP output: Number of exapnsion coefficients to be used
!    -- Now wavelength dependent, 10/20/16

    integer, INTENT (OUT)   ::  n_scatmoms(maxwav)

!  Unscaled profiles of optical depth

    real(fpk)   , INTENT (OUT)  :: aertau_unscaled(maxlayers)

!  AOD scaling

    real(fpk)   , INTENT (OUT) :: aod_scaling(maxwav)

!  Aod output, final

    REAL(fpk)   , DIMENSION( maxlayers, maxwav )      , INTENT (OUT)  :: AEROSOL_DELTAU

!  AOPS

    REAL(fpk)   , DIMENSION( maxwav )                 , INTENT (OUT)  :: AEROSOL_SSALBS 
    REAL(fpk)   , DIMENSION( 6, 0:maxaermoms, maxwav ), INTENT (OUT)  :: AEROSOL_SCATMOMS

!  Aerosol distribution characterisstics.
!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF
   real(fpk),     DIMENSION( 5, 2 )           , INTENT (OUT)  :: AEROSOL_DISTCHARS

!  Linearized AOD Scaling and profile

    real(fpk)   , INTENT (OUT)  :: L_aod_scaling     ( maxwav,    maxaerwfs )
    real(fpk)   , INTENT (OUT)  :: L_aertau_unscaled ( maxlayers, maxaerwfs )

!  Linearized AOPs

    REAL(fpk)   , DIMENSION ( maxlayers, maxwav, maxaerwfs )      , INTENT (OUT) :: L_AEROSOL_DELTAU
    REAL(fpk)   , DIMENSION ( maxwav, maxaerwfs )                 , INTENT (OUT) :: L_AEROSOL_SSALBS 
    REAL(fpk)   , DIMENSION ( 6, 0:maxaermoms, maxwav, maxaerwfs ), INTENT (OUT) :: L_AEROSOL_SCATMOMS

!  Exception handling

   logical,        INTENT (OUT)           :: fail1, fail2
   character*(*),  INTENT (OUT)           :: message_Loading(3)
   character*(*),  INTENT (OUT)           :: Messages_Optical(5)

!  LOCAL VARIABLES
!  @@@@@@@@@@@@@@@

!  Mie/Tmatrix LOCAL INPUT variables (Not part of module input)
!  ============================================================

!      Do_Expcoeffs      - Boolean flag for computing Expansion Coefficients
!      Do_Fmatrix        - Boolean flag for computing F-matrix at equal-angles

   logical    :: Do_Expcoeffs
   logical    :: Do_Fmatrix

!      Do_Monodisperse   - Boolean flag for Doing a Monodisperse calculation
!                          If set, the PSD stuff will be turned off internally

   LOGICAL    :: do_Monodisperse

!  F-matrix Angular control input (NOT REQUIRED HERE)
!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles. (NPNA)
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER    :: n_Fangles
   REAL(fpk)  :: Fangles(max_Mie_angles)

!  Monoradius     - Monodisperse radius size (Microns)

   REAL(fpk)  :: Monoradius

!  NPNA - number of equidistant scattering angles (from 0
!             to 180 deg) for which the scattering matrix is
!             calculated.

!  (Tmatrix only). Style flag.
!    * This is checked and re-set (if required) for Monodisperse case

   logical    :: Do_psd_OldStyle

!  linearization w.r.t refractive index components

   logical    :: Do_LinearRefdummy

!  Linearization for PSD parameters

   logical    :: Do_LinearPSDdummy

!  linearization w.r.t shape factor (T-matrix only)

   !logical    :: Do_LinearEps

!  Lienarization for Loading parameters

   logical    :: do_LinearLoad

!  Lienarization for Bimodal Fraction

   logical    :: do_LinearBmf

!  Aerosol codes OUTPUT variables
!  ==============================

!  10/19/16. Maximum 2 modes for GEMSTOOL.
!  same settings as in the Type-structures file for "GEMSTOOL_MieTmatUser"

   integer, parameter :: maxAerModes = 2

!  10/19/16.Local dimension for the User aerosols.

   integer, parameter :: max_user_angles = max(max_Mie_angles,NPL1)

!  Bulk optical parameters
!     linearizations w.r.t. PSD parameters
!     linearizations w.r.t. Refractive index parameters

   real(fpk)    :: MTU_bulk      (3,  maxAerModes)
   real(fpk)    :: LPSD_MTU_bulk (3,3,MaxAerModes)
   real(fpk)    :: LRFE_MTU_bulk (3,3,MaxAerModes)     ! 3 in second dimension, could be used by Tmat

!jlee added, extinction coefficients at reference wavelength for each mode and regime

   real(fpk)    :: MTU_bulk_ext_ref(maxAerModes) 
   real(fpk)    :: LRFE_MTU_bulk_ext_ref(3,maxAerModes), LPSD_MTU_bulk_ext_ref(3,maxAerModes)

!  Number of Expansion coefficients and Asymmetry parameter

   integer      :: MTU_ncoeffs(maxAerModes)

!  Asymmetry parameter
!     linearizations w.r.t. PSD parameters
!     linearizations w.r.t. Refractive index parameters

   real(fpk)    :: MTU_asymm     (  maxAerModes)
   real(fpk)    :: LPSD_MTU_asymm(3,maxAerModes)
   real(fpk)    :: LRFE_MTU_asymm(3,maxAerModes)     ! 3 in second dimension, could be used by Tmat

!  Distribution parameters ( 1 = Normalization.  2 = Cross-section. 3 = Volume.  4 = REFF.  5 = VEFF )
!     linearizations w.r.t. PSD parameters

   real(fpk)    :: MTU_dist      (5,  maxAerModes)
   real(fpk)    :: LPSD_MTU_dist (5,3,maxAerModes)

!  Exception handling

   character*90 :: MTU_messages(3)

!  Aerosol data from Mie subroutine: Expansion coefficients and F-matrix (optional)
!     linearizations w.r.t. PSD parameters
!     linearizations w.r.t. Refractive index parameters

   real(fpk)    :: Mie_expcoeffs      (6,0:max_Mie_angles,  maxAerModes)
   real(fpk)    :: LPSD_Mie_expcoeffs (6,0:max_Mie_angles,3,maxAerModes)
   real(fpk)    :: LRFE_Mie_expcoeffs (6,0:max_Mie_angles,2,maxAerModes)

   real(fpk)    :: Mie_Fmatrix      (4,max_Mie_angles,  maxAerModes)
   real(fpk)    :: LPSD_Mie_Fmatrix (4,max_Mie_angles,3,maxAerModes)
   real(fpk)    :: LRFE_Mie_Fmatrix (4,max_Mie_angles,2,maxAerModes)

!  Aerosol data from Tmatrix subroutine: Expansion coefficients and F-matrix (optional)
!     linearizations w.r.t. PSD parameters
!     linearizations w.r.t. Refractive index parameters

   real(fpk)    :: Tmat_expcoeffs      (NPL1,6,  maxAerModes)
   real(fpk)    :: LPSD_Tmat_expcoeffs (NPL1,6,3,maxAerModes)
   real(fpk)    :: LRFE_Tmat_expcoeffs (NPL1,6,3,maxAerModes)

   real(fpk)    :: Tmat_Fmatrix      (MAXNPA,6,  maxAerModes)
   real(fpk)    :: LPSD_Tmat_Fmatrix (MAXNPA,6,3,maxAerModes)
   real(fpk)    :: LRFE_Tmat_Fmatrix (MAXNPA,6,3,maxAerModes)

!  Other aerosol data-related variables
!  ====================================

!  Aerosol data from user input file: Expansion coefficients

   real(fpk)    :: User_expcoeffs      (6,0:max_user_angles,  maxAerModes)
   real(fpk)    :: LPSD_User_expcoeffs (6,0:max_user_angles,3,maxAerModes)
   real(fpk)    :: LRFE_User_expcoeffs (6,0:max_user_angles,3,maxAerModes)

!  Universal aerosol file output

   real(fpk)    :: MTU_expcoeffs      (6,0:max_user_angles,  maxAerModes)
   real(fpk)    :: LPSD_MTU_expcoeffs (6,0:max_user_angles,3,maxAerModes)
   real(fpk)    :: LRFE_MTU_expcoeffs (6,0:max_user_angles,3,maxAerModes)

!  Aerosol I/O data file names

   character*90 :: AeroInFile,      AeroOutFile
   character*90 :: AeroInFile_Plus, AeroOutFile_Plus

!  Local variables
!  ===============

!  fixed Proxies

   integer      :: Tmat_Sphtype, Tmat_ndgs, Tmat_nkmax, nblocks, nweights
   real(fpk)    :: Tmat_accuracy, R1, R2, R1R2_cut, xlimit
   logical      :: FixR1R2, Do_EqSaSphere

!  regime-dependent Local aerosol properties

   integer      :: PSDindex
   real(fpk)    :: PSDpars(3)
   real(fpk)    :: nreal, nimag
   real(fpk)    :: Tmateps

!  AOP output: Reference values
!    * These are the values at reference wavelength w0
!    * They are set for wavelength w0, then used again

   real(fpk)    :: extinction_ref, extinction
   real(fpk)    :: L_extinction_ref(maxaerwfs), L_extinction(maxaerwfs)

!  Help Variables
!  --------------

!  revised 10/19/16 to include many new variables

!  Local mask

   integer :: distpars_mask(9)
   data distpars_mask / 2,2,3,2,3,3,3,3,0 /

   character*10 :: ctype
   character*2  :: cnm
   character*4  :: cwav
   character*9  :: AerChar, CharSrc(3)

   logical      :: RefWvl, do_bimodal
   logical      :: MakeAeroFile(maxAerModes)
!   logical      :: ReadAeroFile(maxAerModes)

   integer      :: NLINTMAT, q, qm, qm1, qm2, q1, q2, qmt, nmodes
   integer      :: NPSDdummy(2), NLINdummy(2)
   integer      :: k, n, nm, l, istatus, point_index, w, wc, m, wavmask(maxwav)

   integer      :: OpropSrc(maxAerModes)
   integer      :: InFileUnit(maxAerModes),OutFileUnit(maxAerModes)

   real(fpk)    :: Csca_total,Csca_frac(maxAerModes),sca_wgt(maxAerModes)
   real(fpk)    :: wavelength, Micron_wavelength
   real(fpk)    :: frac(maxAerModes)
   real(fpk)    :: ParsListDummy(10,2)

   real(fpk)    :: LRFE_Csca_frac(3,maxAerModes),LPSD_Csca_frac(3,maxAerModes), deriv(3), div(2)
   real(fpk)    :: L_Csca_total(maxaerwfs),L_sca_wgt(maxaerwfs,maxAerModes)

   integer      :: n_scatmoms_w, n_aerwavs, local_nwav, wa, wa1, wa2, wastart, wdum
   real(fpk)    :: fa1, fa2, lam, lamstart, lamfinish
   logical      :: trawl
   logical      :: do_aer_Jacobians

!  Interpolation arrays
!     UV-Vis : 56  is enough wavelengths for 10 nm intervals over a 250-800  nm range
!     SWIr   : 136 is enough wavelengths for 10 nm intervals over a 750-2100 nm range
!    NOTE    : CAN MAKE THESE ARRAYS ALLOCATABLE ----------THINK ABOUT IT !!!!!!!!!
!    NOTE    : 16 is the maximum number of aerosol weighting functions

   real(fpk)    :: aerwav(136)
   real(fpk)    :: local_aodscaling(136),      local_aerssalbs(136),      local_scatmoms(6,0:maxaermoms,136)
   real(fpk)    :: local_L_aodscaling(136,16), local_L_aerssalbs(136,16), local_L_scatmoms(6,0:maxaermoms,136,16)
   integer      :: local_nscatmoms(136)

!mick mod 1/24/2017 - added aerosol file preamble capacity
!  Rob mod 3/7/17 - expanded preamble to include ref wavelength and PSD information (11 lines in all)
!  Rob fix, 3/7/17 - AeroFilePreamble needs to have a mode index
!  Aerosol output file description header
!   integer, parameter :: nPLines = 6
   integer, parameter :: nPLines = 11
   real(fpk)          :: aerwavlo,aerwavhi,wavlo,wav2,wavhi,wavres
   logical            :: DoAeroFilePreamble(maxAerModes)
!   character(len=80)  :: AeroFilePreamble(nPLines)
   character(len=80)  :: AeroFilePreamble(nPLines,maxAerModes)

!  Local control parameters
!  ------------------------

   logical, parameter :: do_iopchecker    = .false.
!   logical, parameter :: do_iopchecker    = .true.

!  New Local control parameters 10/19/16
!  -------------------------------------

   logical, parameter :: jlee_method    = .false.
!   logical, parameter :: jlee_method    = .true.

!  verbose and debug parameters

   logical :: Verbose  = .true.
   logical :: Debug    = .true.

! Tmat_Verbose flag replaces Tmat_Progress

   logical :: Tmat_Verbose  = .true.

!  Initialize output
!  =================

!  Initialize exception handling

   fail1 = .false.
   fail2 = .false.
   Message_Loading  = ' '
   Messages_optical = ' '

!  Initialize Aerosol Loading

   Loading        = 0.0d0
   Dloading_Dtau  = 0.0d0
   Dloading_Dpars = 0.0d0
   aerlayerflags = .false.

!  Initialize optical properties

   extinction_ref    = 0.0d0
   aod_scaling       = 0.0d0
   aertau_unscaled   = 0.0d0
   aerosol_deltau    = 0.0d0
   aerosol_ssalbs    = 0.0d0
   aerosol_scatmoms  = 0.0d0
   aerosol_distchars = 0.0d0

   L_aod_scaling      = 0.0d0
   L_aertau_unscaled  = 0.0d0
   L_aerosol_deltau   = 0.0d0
   L_aerosol_ssalbs   = 0.0d0
   L_aerosol_scatmoms = 0.0d0

!  Initialize nscatmoms

   n_scatmoms = 50000

!  Set local Mie/Tmatrix variables

   do_monodisperse = .false.             ! Always false
   Do_Fmatrix      = .false.             ! Always false
   DO_Expcoeffs    = .false.             ! this will be set later on
   monoradius      =  0.0d0
   do_psd_oldstyle = .false.

!  Preliminary bookkeeping, necessary before we start
!  --------------------------------------------------

!  Number of modes (1 or 2). Unimodal or Bimodal only.

   nmodes = 1 ; frac(1) = 1.0_fpk
   do_bimodal = GEMSTOOL_INPUTS%MieTmatUser%do_bimodal
   if ( do_bimodal ) then
      nmodes = 2
      frac(1) =  GEMSTOOL_INPUTS%MieTmatUser%bimodal_fraction
      frac(2) =  1.0_fpk - frac(1)
   endif

!  Bimodal and Loading Jacobian flags

   do_LinearBmf  = do_bimodal .and. GEMSTOOL_INPUTS%LinControl%do_AerBulk_BmodFrac_Jacobian
   do_LinearLoad =                  GEMSTOOL_INPUTS%LinControl%do_AerBulk_LoadPars_Jacobians  !! CHECK!! Y.Jung

!  Proxies: local Mie/Tmatrix Linearization variables, for the 5 input GROUPS
!    -- These will emerge from the file-reads, when User-defined aerosols are in force
!mick fix 1/24/2017 - variables now all defined in subroutine "GEMSTOOL_AERWFS_BOOKKEEP"
!                     and then passed in (this code block turned off)

   !if ( GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then
   !   !Use control inputs from user aerosol file
   !else
   !   !Use control inputs GEMSTOOL config files (mainly "GEMSTOOL_LinControl_AerBulk.cfg")
   !   do_LinearRef  = GEMSTOOL_INPUTS%LinControl%do_AerBulk_RefIndex_Jacobians
   !   do_LinearPSD  = GEMSTOOL_INPUTS%LinControl%do_AerBulk_SizeDist_Jacobians
   !   do_LinearEps  = GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols .and. &
   !                   GEMSTOOL_INPUTS%LinControl%do_AerBulk_ShapeFac_Jacobians
   !   nlin(1) = 2 ; if (do_LinearEps) nlin(1) = nlin(1) + 1
   !   npsd(1) = distpars_mask(GEMSTOOL_INPUTS%MieTmatUser%PSDIndex(1))
   !   if ( do_bimodal ) then
   !      nlin(2) = 2 ; if (do_LinearEps) nlin(2) = nlin(2) + 1
   !      npsd(2) = distpars_mask(GEMSTOOL_INPUTS%MieTmatUser%PSDIndex(2))
   !   endif
   !endif

!  Overall Jacobian flag

   do_aer_Jacobians = do_LinearBmf .or. do_LinearLoad .or. do_LinearRef .or. do_LinearPSD

!  Loading offsets

   if ( do_aer_Jacobians ) then
      qm1 = SUM(n_AerBulk_Group_wfs(2:5))
      qm2 = qm1 + 1
      qmt = qm1 + n_AerBulk_Group_wfs(1)
   endif

! Check this
   if ( qmt .ne. n_AerBulk_Total_wfs ) then
      write(*,*) qmt,n_AerBulk_Total_wfs ; stop 'bookekeeping incorrect'
   endif
   qmt = n_AerBulk_Total_wfs 

!  ============================
!  Now form the aerosol Loading
!  ============================

!    Loading = optical depth profile

!    Derivatives : 
!       Cases 1-3 : dloading_Dtau      = derivative of profile w.r.t the total aerosol optical depth at wavelength w0
!       Case 2    : dloading_Dpars(1)  = derivative of profile w.r.t the relaxation parameter (Exponential)
!       Case 3    : dloading_Dpars(1)  = derivative of profile w.r.t the GDF Peak height
!       Case 3    : dloading_Dpars(2)  = derivative of profile w.r.t the GDF Half Width 

!  Case 1: Uniform layer of aerosol
!  ********************************

!   write(*,*)GEMSTOOL_INPUTS%AerLoad%loading_case ; pause 'GRONK'

   if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 1 ) then

      CALL profiles_uniform  &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,  & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_upperboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_lowerboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%aertau_input_w0,            & ! Inputs
            Loading, Dloading_Dtau,                             & ! output
            fail1, message_Loading(1), message_Loading(2) )       ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'Uniform aerosol Loading failed'

!  Case 2: Exponential decay profile
!  *********************************

   else if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 2 ) then

      CALL profiles_expone &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,  & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_upperboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_lowerboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%exploading_relaxation,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%aertau_input_w0,            & ! Inputs
            Loading, Dloading_Dtau, Dloading_Dpars(:,1),        & ! output
            fail1, message_Loading(1), message_Loading(2) )       ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'Exponential aerosol Loading failed'

!  Case 3: GDF (quasi-Gaussian) profile
!  ************************************

   else if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 3 ) then

      CALL profiles_gdfone &
          ( maxlayers, nlayers, height_grid, do_aer_Jacobians,  & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_upperboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%gdfloading_peakheight,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%loading_lowerboundary,      & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%gdfloading_halfwidth,       & ! Inputs
            GEMSTOOL_INPUTS%AerLoad%aertau_input_w0,            & ! Inputs
            Loading, Dloading_Dtau, Dloading_Dpars(:,1), Dloading_Dpars(:,2),     & ! output
            fail1, message_Loading(1), message_Loading(2) )                         ! Exception Handling

      if ( fail1 ) message_Loading(3) = 'GDF aerosol Loading failed'

   endif

!  Return if failure at this stage

   if ( fail1 ) return

!  Assign Unscaled loadings

   do n = 1, nlayers
      aerlayerflags(n)   = ( Loading(n) .ne. 0.0d0  )
      aertau_unscaled(n) =   loading(n)
   enddo

!  Assign Linearizations of Unscaled Loadings

   if ( do_aer_Jacobians ) then
      do n = 1, nlayers
         aertau_unscaled(n) = loading(n)
!         do q = 1, qm1
!            L_aertau_unscaled(n,q) = loading(n)
!         enddo
         L_aertau_unscaled(n,qm1+1) = Dloading_Dtau(n)
         if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 2 ) then
            L_aertau_unscaled(n,qm1+2) = Dloading_Dpars(n,1)
         else if ( GEMSTOOL_INPUTS%AerLoad%loading_case .eq. 3 ) then
            L_aertau_unscaled(n,qm1+2) = Dloading_Dpars(n,1)
            L_aertau_unscaled(n,qm1+3) = Dloading_Dpars(n,2)
         endif
      enddo
   endif

!  debug

!   do n = 1, nlayers
!      if (aerlayerflags(n) ) write(*,*)n,aertau_unscaled(n),(L_aertau_unscaled(n,qm1+q),q=1,3)
!   enddo
!   pause 'aerosol loading'

!  Interpolation setup
!  ===================

   if ( interpolate_aerosols ) then
!mick mod 1/24/2017 - provide user feedback in interpolated aerosol case
!                     (here and below)
      !write(*,*)
      write(*,'(a)') '    Interpolating aerosols ...'
      trawl = .true. ; wa = 0
      if ( .not.do_wavnums ) then
         !Wavelength set
         lamstart = 200.0d0 ; lam = lamstart        !  Smallest wavelength is 200 nm (UVN application)
         do while (trawl)
            lam = lam + 20.0d0
            if ( lam .ge. wav(1) ) then
               wa = wa + 1 ; aerwav(wa) = lam - 20.0d0
               if ( lam .gt. wav(nwav) ) then
                  wa = wa + 1 ; aerwav(wa) = lam ; trawl = .false.
               endif
            endif
         enddo
         n_aerwavs = wa ; local_nwav = n_aerwavs
      else
         !Wavenumber set
         lamfinish = 2100.0d0 ; lam = lamfinish     !  Largest wavelength is 2100 nm (NSW application)
         do while (trawl)
            lam = lam - 20.0d0
            if ( lam .le. 1.0d+07/wav(1) ) then
               wa = wa + 1 ; aerwav(wa) = lam + 20.0d0
               if ( lam .lt. 1.0d+07/wav(nwav) ) then
                  wa = wa + 1 ; aerwav(wa) = lam ; trawl = .false.
               endif
            endif
         enddo
         n_aerwavs = wa ; local_nwav = n_aerwavs
      endif
      write(*,'(a,i2)')   '    Number of aerosol optical property wavelengths: ',n_aerwavs
      if ( .not. do_wavnums ) then
         write(*,'(a,f8.3)') '    Lo wavelength (nm): ',aerwav(1)
         write(*,'(a,f8.3)') '    Hi wavelength (nm): ',aerwav(n_aerwavs)
      else
         write(*,'(a,f8.3)') '    Lo wavelength (nm): ',aerwav(n_aerwavs)
         write(*,'(a,f8.3)') '    Hi wavelength (nm): ',aerwav(1)
      endif
      write(*,'(a,f8.3)') '    Resolution    (nm): ',20.0d0
   else
      local_nwav = nwav
   endif

!  Check to see if reference wavelength is one of the set. Initialize Mask.

   point_index = 0
   if ( .not. do_wavnums ) then
      !Wavelength setup
      if ( interpolate_aerosols ) then
         do w = 1, n_aerwavs
            wavmask(w) = w ; if ( aerwav(w) .eq. GEMSTOOL_INPUTS%AerLoad%reference_w0 ) point_index = w
         enddo
      else
         do w = 1, nwav
            wavmask(w) = w ; if ( wav(w)    .eq. GEMSTOOL_INPUTS%AerLoad%reference_w0 ) point_index = w
         enddo
      endif
   else
      !Wavenumber setup
      if ( interpolate_aerosols ) then
         do w = 1, n_aerwavs
            wavmask(w) = w ; if ( aerwav(w) .eq. GEMSTOOL_INPUTS%AerLoad%reference_w0 ) point_index = w
         enddo
      else
         do w = 1, nwav
            wavmask(w) = w
         enddo
      endif
   endif

!  Mask to use IF reference wavelength IS one of list of wavelengths. [UVN only]

   if ( point_index .ne. 0 ) then
     wavmask(1) = point_index
     wc = 1
     do w = 1, point_index - 1
       wc = wc + 1 ; wavmask(wc) = w
     enddo
     do w = point_index + 1, local_nwav
       wc = wc + 1 ; wavmask(wc) = w
     enddo
   endif

!  ==================
!  OPTICAL PROPERTIES
!  ==================

!  Make aerosol file flags. These will be set if you are doing Mie or Tmatrix

    MakeAeroFile = .true.
    if ( GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then
       MakeAeroFile(1:nmodes) = .false.
    endif
!mick mod 1/24/2017 - added aerosol file preamble capacity
    DoAeroFilePreamble = .false.

!  Optical sources

   OpropSrc = 0 ; CharSrc = ' '
   if ( GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols ) then
      OpropSrc(1) = 1 ; if ( GEMSTOOL_INPUTS%MieTmatUser%do_bimodal) OpropSrc(2) = 1 ; CharSrc(1) = '(MieScat)'
   else if ( GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols ) then
      OpropSrc(1) = 2 ; if ( GEMSTOOL_INPUTS%MieTmatUser%do_bimodal) OpropSrc(2) = 2 ; CharSrc(2) = '(Tmatrix)'
   else
      OpropSrc(1) = 3 ; if ( GEMSTOOL_INPUTS%MieTmatUser%do_bimodal) OpropSrc(2) = 3 ; CharSrc(3) = '(UserDef)'
   endif
 
!  First, open files for READING user-defined aerosol optical property input if requested
!  (note: these files closed @ ~line 2650)
!mick mod 1/24/2017 - added filename "Debug" output

   if ( GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then
      do nm = 1, nmodes
         call Get_FileUnit(InFileUnit(nm))
         AeroInFile='UserAeroFiles/in/' // trim(GEMSTOOL_INPUTS%MieTmatUser%User_Aerofile_names(nm)) 
         if ( do_Aer_Jacobians ) then
            AeroInFile_plus = Trim(AeroInFile)//'_PLUS'
            open(unit=InFileUnit(nm), file=trim(AeroInFile_plus),iostat=istatus,status='old')
            if ( Debug ) then
              if (nm .eq. 1) write(*,*)
              write(*,*) '   Using ' // trim(AeroInFile_plus)
            endif
         else
            open(unit=InFileUnit(nm), file=trim(AeroInFile),iostat=istatus,status='old')
            if ( Debug ) then
              if (nm .eq. 1) write(*,*)
              write(*,*) '   Using ' // trim(AeroInFile)
            endif
         endif
         if ( istatus > 0 ) then
            Messages_Optical(1) = 'Failure: Aerosol optical property input file'
            Messages_Optical(2) = 'File   : ' // trim(AeroInFile)
            Messages_Optical(3) = 'Problem: Error opening file'
            fail2 = .true. ; close(unit=InFileUnit(nm)) ; return
         endif

!mick mod 1/24/2017 - added aerosol file preamble capacity
         DoAeroFilePreamble(nm) = .true.
      enddo
   endif

!  Second, open files for WRITING Mie or T-matrix aerosol optical property output
!   - This is independent of the above read-user code

   if ( .not. GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then
     do nm = 1, nmodes
       if ( MakeAeroFile(nm) ) then
         call Get_FileUnit(OutFileUnit(nm)) ; write(cnm,'(i2.2)') nm
         if ( OpropSrc(nm) .eq. 1 ) then
            AeroOutFile='UserAeroFiles/out/' // 'TestProps_mode' // cnm // '_Mie_prop.dat'
         elseif ( OpropSrc(nm) .eq. 2 ) then
            AeroOutFile='UserAeroFiles/out/' // 'TestProps_mode' // cnm // '_Tmat_prop.dat'
         else
            AeroOutFile='UserAeroFiles/out/' // 'TestProps_mode' // cnm // '_User_prop.dat'
         endif
         if ( do_Aer_Jacobians ) then
            AeroOutFile_plus = Trim(AeroOutFile)//'_PLUS'
            open(unit=OutFileUnit(nm), file=trim(AeroOutFile_plus), status='replace')
         else
            open(unit=OutFileUnit(nm), file=trim(AeroOutFile), status='replace')
         endif

!mick mod 1/24/2017 - added aerosol file preamble capacity
!  Define aerosol optical property file preamble

         if ( .not. do_wavnums ) then
            !already in nm
            wavlo = wav(1)
            wav2  = wav(2)
            wavhi = wav(nwav)

            aerwavlo = aerwav(1)
            aerwavhi = aerwav(n_aerwavs)
         else
            !convert from wavenumbers in cm^-1 to wavelengths in nm
            wavlo = 1.0d+07/wav(nwav)   
            wav2  = 1.0d+07/wav(nwav-1)
            wavhi = 1.0d+07/wav(1)

            aerwavlo = aerwav(n_aerwavs)
            aerwavhi = aerwav(1)
         endif
         wavres = dabs(wav2-wavlo)

         DoAeroFilePreamble(nm) = .true.

!  Original 6-lines Preamble.
!         write(AeroFilePreamble(1),'(6x,a)')    'Inputs used to create aerosol file (all wvls in nm):'
!         write(AeroFilePreamble(2),'(6x,a,l1)') 'Aerosols interpolated: ',interpolate_aerosols
!         if ( interpolate_aerosols ) then
!            write(AeroFilePreamble(3),'(6x,2(a,f8.3))') 'Band Lo Wvl  : ',wavlo ,'   Aero Lo Wvl  : ',aerwavlo
!            write(AeroFilePreamble(4),'(6x,2(a,f8.3))') 'Band Hi Wvl  : ',wavhi ,'   Aero Hi Wvl  : ',aerwavhi
!            write(AeroFilePreamble(5),'(6x,2(a,f8.3))') 'Band Res     : ',wavres,'   Aero Res     : ',20.0d0
!            write(AeroFilePreamble(6),'(6x,2(a,i8  ))') 'Band Num Wvls: ',nwav  ,'   Aero Num Wvls: ',n_aerwavs
!         else 
!            write(AeroFilePreamble(3),'(6x,a,f8.3)')    'Band Lo Wvl  : ',wavlo
!            write(AeroFilePreamble(4),'(6x,a,f8.3)')    'Band Hi Wvl  : ',wavhi
!            write(AeroFilePreamble(5),'(6x,a,f8.3)')    'Band Res     : ',wavres
!            write(AeroFilePreamble(6),'(6x,a,i8  )')    'Band Num Wvls: ',nwav
!         endif

!  New 11-line preamble

         write(AeroFilePreamble(1,nm),'(6x,a)')    'Inputs used to create aerosol file (all wvls in nm):'
         write(AeroFilePreamble(2,nm),'(6x,a,l1)') 'Aerosols interpolated: ',interpolate_aerosols
         write(AeroFilePreamble(3,nm),'(6x,a,f8.3)') 'Aerosol RefWavelength: ',GEMSTOOL_INPUTS%AerLoad%reference_w0
         if ( interpolate_aerosols ) then
            write(AeroFilePreamble(4,nm),'(6x,2(a,f8.3))') 'Band Lo Wvl  : ',wavlo ,'   Aero Lo Wvl  : ',aerwavlo
            write(AeroFilePreamble(5,nm),'(6x,2(a,f8.3))') 'Band Hi Wvl  : ',wavhi ,'   Aero Hi Wvl  : ',aerwavhi
            write(AeroFilePreamble(6,nm),'(6x,2(a,f8.3))') 'Band Res     : ',wavres,'   Aero Res     : ',20.0d0
            write(AeroFilePreamble(7,nm),'(6x,2(a,i8  ))') 'Band Num Wvls: ',nwav  ,'   Aero Num Wvls: ',n_aerwavs
         else 
            write(AeroFilePreamble(4,nm),'(6x,a,f8.3)')    'Band Lo Wvl  : ',wavlo
            write(AeroFilePreamble(5,nm),'(6x,a,f8.3)')    'Band Hi Wvl  : ',wavhi
            write(AeroFilePreamble(6,nm),'(6x,a,f8.3)')    'Band Res     : ',wavres
            write(AeroFilePreamble(7,nm),'(6x,a,i8  )')    'Band Num Wvls: ',nwav
         endif
         PSDIndex = GEMSTOOL_INPUTS%MieTmatUser%PSDIndex(nm)
         PSDpars  = GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1:3,nm)
         write(AeroFilePreamble(8,nm), '(6x,a,i9)')      'PSD Index      : ',PSDIndex
         write(AeroFilePreamble(9,nm), '(6x,a,f9.5)')    'PSD parameter 1: ',PSDpars(1)
         write(AeroFilePreamble(10,nm),'(6x,a,f9.5)')    'PSD parameter 2: ',PSDpars(2)
         write(AeroFilePreamble(11,nm),'(6x,a,f9.5)')    'PSD parameter 3: ',PSDpars(3)

       endif
     enddo
   endif

!  Prepare the reference-wavelength Mie/Tmatrix/User Inputs
!  ========================================================

   if ( point_index .eq. 0 ) then

!  Only require extinction coefficient if flagged
!  Set the local Mie program inputs (bulk properties only)

      Do_Expcoeffs = .false.
      RefWvl       = .true.

!  Dummy index for writing

      wdum = 0

!  Reference wavelength

      wavelength        = GEMSTOOL_INPUTS%AerLoad%reference_w0
      Micron_wavelength = wavelength/1000.0d0

!  Progress. 2/2/15 Format changed from f8.3 to f12.5

      if ( Verbose .or. Debug ) then
         if ( Debug ) write(*,*)
         if ( do_aer_jacobians ) then
            write(*,'(4x,a,5x,f12.5)') 'Aerosol Linearized calculation - Doing reference wavelength   : ',wavelength
         else 
            write(*,'(4x,a,5x,f12.5)') 'Aerosol Regular calculation - Doing reference wavelength   : ',wavelength
         endif
      endif

!  Reference extinction

      extinction_ref   = 0.0_fpk
      L_extinction_ref = 0.0_fpk

!  Jacobian counting

      qm = 0

!  Mode loop

      do nm = 1, nmodes

!  Aerosol type

         AerChar = CharSrc(OPropSrc(nm))

!  Mie/Tmat Proxies, mode dependent
!mick fix 1/24/2017 - NLIN defined in subroutine "GEMSTOOL_AERWFS_BOOKKEEP" now

         if ( OPropSrc(nm) .eq. 1 .or.  OPropSrc(nm) .eq. 2 ) then
            PSDIndex = GEMSTOOL_INPUTS%MieTmatUser%PSDIndex(nm)
            PSDpars  = GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1:3,nm)
            R1       = GEMSTOOL_INPUTS%MieTmatUser%R1(nm)
            R2       = GEMSTOOL_INPUTS%MieTmatUser%R2(nm)
            FixR1R2  = GEMSTOOL_INPUTS%MieTmatUser%FixR1R2(nm)
            nreal    = GEMSTOOL_INPUTS%MieTmatUser%nreal(nm)
            nimag    = GEMSTOOL_INPUTS%MieTmatUser%nimag(nm)
         endif

         if ( OPropSrc(nm) .eq. 1 ) then
            !NLIN(nm)  = 2  ! Number of weighting functions
            nblocks   = GEMSTOOL_INPUTS%MieTmatUser%nblocks(nm)
            nweights  = GEMSTOOL_INPUTS%MieTmatUser%nweights(nm)
            xlimit    = GEMSTOOL_INPUTS%MieTmatUser%xparticle_limit
            R1R2_cut  = GEMSTOOL_INPUTS%MieTmatUser%R1R2_cutoff(nm)
         else if ( OPropSrc(nm) .eq. 2 ) then
            !NLIN(nm)  = 2  ! Number of weighting functions
            !if ( do_LinearEps) NLIN(nm)  = 3  ! Number of weighting functions Tmatrix 
            Do_EqSaSphere = GEMSTOOL_INPUTS%MieTmatUser%Do_EqSaSphere
            Tmat_Sphtype  = GEMSTOOL_INPUTS%MieTmatUser%Tmat_Sphtype
            Tmateps       = GEMSTOOL_INPUTS%MieTmatUser%Tmat_eps(nm)
            Tmat_nkmax    = GEMSTOOL_INPUTS%MieTmatUser%Tmat_nkmax(nm)
            Tmat_ndgs     = GEMSTOOL_INPUTS%MieTmatUser%Tmat_ndgs(nm)
            Tmat_accuracy = GEMSTOOL_INPUTS%MieTmatUser%Tmat_accuracy  
         endif
         
!  Progress

         if ( Verbose .or. Debug ) then
            write(*,'(5x,a,1x,i3,1x,A9)') '-- Doing mode # ',nm, charSrc(OpropSrc(nm))
         endif

!  SINGLE CALL. [Tmat_Verbose flag introduced 6/6/14]

         if ( do_aer_Jacobians ) then
           if ( OpropSrc(nm) .eq. 1 ) then  !Mie
             call RTSMie_main_plus  &
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_LinearRef, Do_LinearPSD,         & ! Input top level flags
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, nblocks, nweights,             & ! Inputs (PSD stuff)
                xlimit, R1R2_cut, n_Fangles, Fangles, Micron_wavelength, nreal, nimag,         & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm), Mie_expcoeffs(:,:,nm),         & ! Outputs regular
                Mie_Fmatrix(:,:,nm), MTU_dist(:,nm),                                           & ! Outputs regular
                LPSD_MTU_bulk(:,:,nm),   LPSD_MTU_asymm(:,nm),   LPSD_Mie_expcoeffs(:,:,:,nm), & ! Outputs linearized
                LPSD_Mie_Fmatrix(:,:,:,nm), LPSD_MTU_dist(:,:,nm),                             & ! Outputs linearized
                LRFE_MTU_bulk(:,1:2,nm), LRFE_MTU_asymm(1:2,nm), LRFE_Mie_expcoeffs(:,:,:,nm), & ! Outputs linearized
                LRFE_Mie_Fmatrix(:,:,:,nm),                                                    & ! Outputs linearized
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
           elseif ( OpropSrc(nm) .eq. 2 ) then !Tmat
!mick fix 1/24/2017 - modified some array extents in the call
             !call tmat_master_plus ( Tmat_Verbose, &
             !   Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_EqSaSphere,                      & ! Input top level flags
             !   Do_LinearRef, Do_LinearEps, Do_LinearPSD, Do_psd_OldStyle,                     & ! Inputs (Flags)
             !   PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, Tmat_Sphtype, Tmat_nkmax,      & ! Inputs (PSD stuff)
             !   n_Fangles, Tmat_ndgs, Tmateps, Tmat_accuracy, Micron_wavelength, nreal, nimag, & ! Inputs (optical + angles)
             !   MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm),                                & ! Outputs Regular
             !   Tmat_expcoeffs(:,:,nm), Tmat_Fmatrix(:,:,nm),                                  & ! Outputs Regular
             !   LPSD_MTU_bulk(:,:,nm),         LPSD_MTU_asymm(:,nm),                           & ! Outputs linearized
             !   LPSD_Tmat_expcoeffs(:,:,:,nm), LPSD_Tmat_Fmatrix(:,:,:,nm),                    & ! Outputs linearized
             !   LRFE_MTU_bulk(:,1:3,nm),       LRFE_MTU_asymm(1:2,nm),                         & ! Outputs linearized
             !   LRFE_Tmat_expcoeffs(:,:,:,nm), LRFE_Tmat_Fmatrix(:,:,:,nm),                    & ! Outputs linearized
             !   NLINTMAT, MTU_dist(:,nm), LPSD_MTU_dist(:,:,nm),                               & ! Outputs (PSD)
             !   fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
             call tmat_master_plus ( Tmat_Verbose, &
                Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_EqSaSphere,                      & ! Input top level flags
                Do_LinearRef, Do_LinearEps, Do_LinearPSD, Do_psd_OldStyle,                     & ! Inputs (Flags)
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, Tmat_Sphtype, Tmat_nkmax,      & ! Inputs (PSD stuff)
                n_Fangles, Tmat_ndgs, Tmateps, Tmat_accuracy, Micron_wavelength, nreal, nimag, & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm),                                & ! Outputs Regular
                Tmat_expcoeffs(:,:,nm), Tmat_Fmatrix(:,:,nm),                                  & ! Outputs Regular
                LPSD_MTU_bulk(:,:,nm),         LPSD_MTU_asymm(:,nm),                           & ! Outputs linearized
                LPSD_Tmat_expcoeffs(:,:,:,nm), LPSD_Tmat_Fmatrix(:,:,:,nm),                    & ! Outputs linearized
                LRFE_MTU_bulk(:,:,nm),         LRFE_MTU_asymm(:,nm),                           & ! Outputs linearized
                LRFE_Tmat_expcoeffs(:,:,:,nm), LRFE_Tmat_Fmatrix(:,:,:,nm),                    & ! Outputs linearized
                NLINTMAT, MTU_dist(:,nm), LPSD_MTU_dist(:,:,nm),                               & ! Outputs (PSD)
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
           else  !User
             call Read_aerosol_oprop_file_plus  &
                ( max_user_angles, InFileUnit(nm), Micron_wavelength, RefWvl,                                        & ! Inputs
                  maxAerModes, nm, nPLines, DoAeroFilePreamble,                                                      & ! I/O
                  do_LinearRefdummy, do_LinearPSDdummy, NLINdummy(NM), NPSDdummy(NM), ParsListDummy(1,nm),           & ! Outputs
                  MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), User_expcoeffs(1,0,nm),            & ! Outputs
                  LRFE_MTU_bulk(1,1,nm), LRFE_MTU_asymm(1,nm), LRFE_User_expcoeffs(1,0,1,nm),                        & ! Outputs
                  LPSD_MTU_bulk(1,1,nm), LPSD_MTU_asymm(1,nm), LPSD_User_expcoeffs(1,0,1,nm), LPSD_MTU_dist(1,1,nm), & ! Outputs
                  fail2, MTU_messages )                                                                     ! Exception handling
!mick fix 1/24/2017 - correct error handling of file in use
             !if ( fail2 ) MTU_messages(2) = 'File   : ' // trim(GEMSTOOL_INPUTS%MieTmatUser%User_Aerofile_names(nm))
             if ( fail2 ) MTU_messages(2) = 'File   : ' // trim(AeroInFile_plus)
           endif
         else
           if ( OpropSrc(nm) .eq. 1 ) then  !Mie
             call RTSMie_main  &
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,                                     & ! Input top level flags
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, nblocks, nweights,             & ! Inputs (PSD stuff)
                xlimit, R1R2_cut, n_Fangles, Fangles, Micron_wavelength, nreal, nimag,         & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm), Mie_expcoeffs(:,:,nm),         & ! Outputs
                Mie_Fmatrix(:,:,nm), MTU_dist(:,nm),                                           & ! Outputs
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
           elseif ( OpropSrc(nm) .eq. 2 ) then !Tmat
             call tmat_master ( Tmat_Verbose, &
                Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_EqSaSphere, Do_psd_OldStyle,     & ! Input top level flags
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, Tmat_Sphtype, Tmat_nkmax,      & ! Inputs (PSD stuff)
                n_Fangles, Tmat_ndgs, Tmateps, Tmat_accuracy, Micron_wavelength, nreal, nimag, & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm), Tmat_expcoeffs(:,:,nm),        & ! Outputs
                Tmat_Fmatrix(:,:,nm), MTU_dist(:,nm),                                          & ! Outputs
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
           else  !User
             call Read_aerosol_oprop_file  &
              ( max_user_angles, InFileUnit(nm), Micron_wavelength, RefWvl,                             & ! Inputs
                maxAerModes, nm, nPLines, DoAeroFilePreamble,                                           & ! I/O
                MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), User_expcoeffs(1,0,nm), & ! Outputs
                fail2, MTU_messages )                                                                     ! Exception handling
             if ( fail2 ) MTU_messages(2) = 'File   : ' // trim(GEMSTOOL_INPUTS%MieTmatUser%User_Aerofile_names(nm))
           endif
         endif

!  Exception handling on everything

         if ( Fail2 ) then
            ctype = 'Regular   ' ; if ( do_aer_Jacobians ) ctype = 'Linearized'
            do m = 1, 3
               Messages_Optical(m) = adjustl(trim(MTU_messages(m)))
            enddo
            if ( OpropSrc(nm) .eq. 1 .or. OpropSrc(nm) .eq. 2 ) then
               Messages_Optical(4) = 'First call to '//ctype//' '//AerChar//' program, mode # '//cnm//', reference wavelength'
            else
               Messages_Optical(4) = 'First read from '//ctype//' '//AerChar//' file, mode # '//cnm//', reference wavelength'
            endif
            Messages_Optical(5) = 'Failure for '//ctype//', mode # '//cnm//', reference wavelength calculation'
            return
         endif

!  Coefficients Must be copied. Note special care over Tmat output.
!   This is the reference wavelength calculation, do not need coefficients

         MTU_expcoeffs(1,0,nm)        = 0.0_fpk 
         LRFE_MTU_expcoeffs(1,0,1,nm) = 0.0_fpk 
         LPSD_MTU_expcoeffs(1,0,1,nm) = 0.0_fpk

!  Write aerosol optical properties to file (optional)
 
         if ( MakeAeroFile(nm) ) then
            if ( do_aer_Jacobians ) then
               call Write_aerosol_oprop_file_plus  &
                 ( max_user_angles, OutFileUnit(nm), wdum, Micron_wavelength, RefWvl,                               & ! Inputs
                   maxAerModes, nm, nPLines, DoAeroFilePreamble, AeroFilePreamble,                                  & ! I/O
                   do_LinearRef, do_LinearPSD, NLIN(nm), NPSD(nm), ParsList(1,nm),                                  & ! Inputs
                   MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), MTU_expcoeffs(1,0,nm),           & ! Inputs
                   LRFE_MTU_bulk(1,1,nm), LRFE_MTU_asymm(1,nm), LRFE_MTU_expcoeffs(1,0,1,nm),                       & ! Inputs
                   LPSD_MTU_bulk(1,1,nm), LPSD_MTU_asymm(1,nm), LPSD_MTU_expcoeffs(1,0,1,nm),LPSD_MTU_dist(1,1,nm)  ) ! Inputs
            else
               call Write_aerosol_oprop_file  &
                 ( max_user_angles, OutFileUnit(nm), w, Micron_wavelength, RefWvl,                       & ! Inputs
                   maxAerModes, nm, nPLines, DoAeroFilePreamble, AeroFilePreamble,                       & ! I/O
                   MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), MTU_expcoeffs(1,0,nm) ) ! Inputs
            endif
         endif

!  Set the reference quantities. NOTE the linearization Bookkeeping.

         !Composite extinction @ REF wvl produced by mode extinctions from
         !  the Tmat or Mie code weighted by input mode fractions.
         !  --> Case: ref wvl is SEPARATE FROM the current spectral band 

         MTU_bulk_ext_ref(nm) = MTU_bulk(1,nm) ! jlee added, extinction coefficients at reference wavelength
         extinction_ref = extinction_ref + frac(nm)*MTU_bulk(1,nm) 
         if ( do_aer_Jacobians ) then
            if ( Do_LinearRef) then
               q1 = qm+1 ; q2 = qm + nlin(nm) ; qm = qm + nlin(nm)
               L_extinction_ref(qm1:qm2) = L_extinction_ref(qm1:qm2) + frac(nm) * LRFE_MTU_bulk(1,1:nlin(nm),nm)
            endif
            if ( Do_LinearPSD) then
               q1 = qm+1 ; q2 = qm + npsd(nm) ; qm = qm + npsd(nm)
               L_extinction_ref(qm1:qm2) = L_extinction_ref(qm1:qm2) + frac(nm) * LPSD_MTU_bulk(1,1:npsd(nm),1)
            endif
            if ( nm.eq.2 .and. do_LinearBmf ) then
               qm = qm + 1 ; L_extinction_ref(qm+1) = MTU_bulk(1,1) - MTU_bulk(1,2)
            endif
         endif

!  Set the Distributions

         AEROSOL_DISTCHARS(1:5,nm) = MTU_dist(1:5,nm)

!  End mode loop

      enddo

!  Progress

      if ( Verbose .or. Debug ) then
         write(*,'(4x,a,5x,f12.5)') &
            'Aerosol Linearized calculation - Finished reference wavelength: ',wavelength
      endif

!  End reference-wavelength calculation

   endif

!  Prepare General (all-wavelength, all-wavenumber) Mie/Tmat/User inputs
!  =====================================================================

!  Set the local inputs (general)

   Do_Expcoeffs = .true.
   RefWvl       = .false.

!  Start wavelength/wavenumber loop.  First wavelength will be the reference, if in list.
!  (note: this loop ends @ ~line 2465)

   do wc = 1, local_nwav

!  Wavelengths [nm]

      w = wavmask(wc)

!mick note:
!      aerwav(w) - comes from interpolation setup code above
!      wav(w)    - comes from spectral config file (either "GEMSTOOL_Lambdas.cfg" or "GEMSTOOL_Wavenums.cfg")

      if ( interpolate_aerosols ) then
         wavelength = aerwav(w)
      else
        if ( do_wavnums      ) wavelength = 1.0d+07/wav(w) !convert from cm^-1 to nm
        if ( .not.do_wavnums ) wavelength = wav(w)         !in nm
      endif
      Micron_wavelength = wavelength/1000.0d0 !in um

!  Progress.  2/2/15 Format changed from f8.3 to f12.5

      if ( Verbose.or.Debug ) then
         if ( Debug ) write(*,*)
         if ( do_aer_jacobians ) then
            write(*,'(4x,a,1x,i3,1x,f12.5)') 'Aerosol Linearized calculation - Doing point, wavelength : ', wc, wavelength
         else 
            write(*,'(4x,a,1x,i3,1x,f12.5)') 'Aerosol Regular calculation - Doing point, wavelength    : ', wc, wavelength
         endif
      endif

!  Jacobian counting

      qm = 0

!  Mode loop

      do nm = 1, nmodes

!  Aerosol type

         AerChar = CharSrc(OPropSrc(nm)) ; write(cwav,'(I4)')wc ; write(cnm,'(I2)') nm

!  Mie/Tmat Proxies, mode dependent
!mick fix 1/24/2017 - NLIN defined in subroutine "GEMSTOOL_AERWFS_BOOKKEEP" now

         if ( OPropSrc(nm) .eq. 1 .or.  OPropSrc(nm) .eq. 2 ) then
            PSDIndex = GEMSTOOL_INPUTS%MieTmatUser%PSDIndex(nm)
            PSDpars  = GEMSTOOL_INPUTS%MieTmatUser%PSDpars(1:3,nm)
            R1       = GEMSTOOL_INPUTS%MieTmatUser%R1(nm)
            R2       = GEMSTOOL_INPUTS%MieTmatUser%R2(nm)
            FixR1R2  = GEMSTOOL_INPUTS%MieTmatUser%FixR1R2(nm)
            nreal    = GEMSTOOL_INPUTS%MieTmatUser%nreal(nm)
            nimag    = GEMSTOOL_INPUTS%MieTmatUser%nimag(nm)
         endif
         if ( OPropSrc(nm) .eq. 1 ) then
            !NLIN(nm)  = 2
            nblocks   = GEMSTOOL_INPUTS%MieTmatUser%nblocks(nm)
            nweights  = GEMSTOOL_INPUTS%MieTmatUser%nweights(nm)
            xlimit    = GEMSTOOL_INPUTS%MieTmatUser%xparticle_limit
            R1R2_cut  = GEMSTOOL_INPUTS%MieTmatUser%R1R2_cutoff(nm)
         else if ( OPropSrc(nm) .eq. 2 ) then
            !NLIN(nm)  = 2  ! Number of weighting functions
            !if ( do_LinearEps) NLIN(nm)  = 3  ! Number of weighting functions Tmatrix
!mick fix 1/24/2017 - define "Do_EqSaSphere"
            Do_EqSaSphere = GEMSTOOL_INPUTS%MieTmatUser%Do_EqSaSphere
            Tmat_Sphtype  = GEMSTOOL_INPUTS%MieTmatUser%Tmat_Sphtype
            Tmateps       = GEMSTOOL_INPUTS%MieTmatUser%Tmat_eps(nm)
            Tmat_nkmax    = GEMSTOOL_INPUTS%MieTmatUser%Tmat_nkmax(nm)
            Tmat_ndgs     = GEMSTOOL_INPUTS%MieTmatUser%Tmat_ndgs(nm)
            Tmat_accuracy = GEMSTOOL_INPUTS%MieTmatUser%Tmat_accuracy  
         endif

!  Progress

         if ( Verbose .or. Debug ) then
            write(*,'(5x,a,1x,i3,1x,A9)') '-- Doing mode # ',nm, charSrc(OpropSrc(nm))
         endif

!  SINGLE CALL. [Tmat_Verbose flag introduced 6/6/14]

         if ( do_aer_Jacobians ) then
           if ( OpropSrc(nm) .eq. 1 ) then  !Mie
             call RTSMie_main_plus  &
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_LinearRef, Do_LinearPSD,         & ! Input top level flags
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, nblocks, nweights,             & ! Inputs (PSD stuff)
                xlimit, R1R2_cut, n_Fangles, Fangles, Micron_wavelength, nreal, nimag,         & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm), Mie_expcoeffs(:,:,nm),         & ! Outputs regular
                Mie_Fmatrix(:,:,nm), MTU_dist(:,nm),                                           & ! Outputs regular
                LPSD_MTU_bulk(:,:,nm),   LPSD_MTU_asymm(:,nm),   LPSD_Mie_expcoeffs(:,:,:,nm), & ! Outputs linearized
                LPSD_Mie_Fmatrix(:,:,:,nm), LPSD_MTU_dist(:,:,nm),                             & ! Outputs linearized
                LRFE_MTU_bulk(:,1:2,nm), LRFE_MTU_asymm(1:2,nm), LRFE_Mie_expcoeffs(:,:,:,nm), & ! Outputs linearized
                LRFE_Mie_Fmatrix(:,:,:,nm),                                                    & ! Outputs linearized
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
           elseif ( OpropSrc(nm) .eq. 2 ) then !Tmat
!mick fix 1/24/2017 - modified some array extents in the call
             !call tmat_master_plus ( Tmat_Verbose, &
             !   Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_EqSaSphere,                      & ! Input top level flags
             !   Do_LinearRef, Do_LinearEps, Do_LinearPSD, Do_psd_OldStyle,                     & ! Inputs (Flags)
             !   PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, Tmat_Sphtype, Tmat_nkmax,      & ! Inputs (PSD stuff)
             !   n_Fangles, Tmat_ndgs, Tmateps, Tmat_accuracy, Micron_wavelength, nreal, nimag, & ! Inputs (optical + angles)
             !   MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm),                                & ! Outputs Regular
             !   Tmat_expcoeffs(:,:,nm), Tmat_Fmatrix(:,:,nm),                                  & ! Outputs Regular
             !   LPSD_MTU_bulk(:,:,nm),         LPSD_MTU_asymm(:,nm),                           & ! Outputs linearized
             !   LPSD_Tmat_expcoeffs(:,:,:,nm), LPSD_Tmat_Fmatrix(:,:,:,nm),                    & ! Outputs linearized
             !   LRFE_MTU_bulk(:,1:3,nm),       LRFE_MTU_asymm(1:2,nm),                         & ! Outputs linearized
             !   LRFE_Tmat_expcoeffs(:,:,:,nm), LRFE_Tmat_Fmatrix(:,:,:,nm),                    & ! Outputs linearized
             !   NLINTMAT, MTU_dist(:,nm), LPSD_MTU_dist(:,:,nm),                               & ! Outputs (PSD)
             !   fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
             call tmat_master_plus ( Tmat_Verbose, &
                Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_EqSaSphere,                      & ! Input top level flags
                Do_LinearRef, Do_LinearEps, Do_LinearPSD, Do_psd_OldStyle,                     & ! Inputs (Flags)
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, Tmat_Sphtype, Tmat_nkmax,      & ! Inputs (PSD stuff)
                n_Fangles, Tmat_ndgs, Tmateps, Tmat_accuracy, Micron_wavelength, nreal, nimag, & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm),                                & ! Outputs Regular
                Tmat_expcoeffs(:,:,nm), Tmat_Fmatrix(:,:,nm),                                  & ! Outputs Regular
                LPSD_MTU_bulk(:,:,nm),         LPSD_MTU_asymm(:,nm),                           & ! Outputs linearized
                LPSD_Tmat_expcoeffs(:,:,:,nm), LPSD_Tmat_Fmatrix(:,:,:,nm),                    & ! Outputs linearized
                LRFE_MTU_bulk(:,:,nm),         LRFE_MTU_asymm(:,nm),                           & ! Outputs linearized
                LRFE_Tmat_expcoeffs(:,:,:,nm), LRFE_Tmat_Fmatrix(:,:,:,nm),                    & ! Outputs linearized
                NLINTMAT, MTU_dist(:,nm), LPSD_MTU_dist(:,:,nm),                               & ! Outputs (PSD)
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
           else  !User
             call Read_aerosol_oprop_file_plus  &
                ( max_user_angles, InFileUnit(nm), Micron_wavelength, RefWvl,                                        & ! Inputs
                  maxAerModes, nm, nPLines, DoAeroFilePreamble,                                                      & ! I/O
                  do_LinearRefdummy, do_LinearPSDdummy, NLINdummy(NM), NPSDdummy(NM), ParsListDummy(1,nm),           & ! Outputs
                  MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), User_expcoeffs(1,0,nm),            & ! Outputs
                  LRFE_MTU_bulk(1,1,nm), LRFE_MTU_asymm(1,nm), LRFE_User_expcoeffs(1,0,1,nm),                        & ! Outputs
                  LPSD_MTU_bulk(1,1,nm), LPSD_MTU_asymm(1,nm), LPSD_User_expcoeffs(1,0,1,nm), LPSD_MTU_dist(1,1,nm), & ! Outputs
                  fail2, MTU_messages )                                                                     ! Exception handling
!mick fix 1/24/2017 - correct error handling of file in use
             !if ( fail2 ) MTU_messages(2) = 'File   : ' // trim(GEMSTOOL_INPUTS%MieTmatUser%User_Aerofile_names(nm))
             if ( fail2 ) MTU_messages(2) = 'File   : ' // trim(AeroInFile_plus)
           endif
         else
           if ( OpropSrc(nm) .eq. 1 ) then  !Mie
             call RTSMie_main  &
              ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,                                     & ! Input top level flags
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, nblocks, nweights,             & ! Inputs (PSD stuff)
                xlimit, R1R2_cut, n_Fangles, Fangles, Micron_wavelength, nreal, nimag,         & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm), Mie_expcoeffs(:,:,nm),         & ! Outputs
                Mie_Fmatrix(:,:,nm), MTU_dist(:,nm),                                           & ! Outputs
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
           elseif ( OpropSrc(nm) .eq. 2 ) then !Tmat
             call tmat_master ( Tmat_Verbose, &
                Do_Expcoeffs, Do_Fmatrix, do_Monodisperse, Do_EqSaSphere, Do_psd_OldStyle,     & ! Input top level flags
                PSDIndex, PSDpars, MonoRadius, R1, R2, FixR1R2, Tmat_Sphtype, Tmat_nkmax,      & ! Inputs (PSD stuff)
                n_Fangles, Tmat_ndgs, Tmateps, Tmat_accuracy, Micron_wavelength, nreal, nimag, & ! Inputs (optical + angles)
                MTU_bulk(:,nm), MTU_asymm(nm), MTU_ncoeffs(nm), Tmat_expcoeffs(:,:,nm),        & ! Outputs
                Tmat_Fmatrix(:,:,nm), MTU_dist(:,nm),                                          & ! Outputs
                fail2, istatus, MTU_messages(1), MTU_messages(2), MTU_messages(3) )              ! Exception handling
           else  !User
             call Read_aerosol_oprop_file  &
              ( max_user_angles, InFileUnit(nm), Micron_wavelength, RefWvl,                             & ! Inputs
                maxAerModes, nm, nPLines, DoAeroFilePreamble,                                           & ! I/O
                MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), User_expcoeffs(1,0,nm), & ! Outputs
                fail2, MTU_messages )                                                                     ! Exception handling
             if ( fail2 ) MTU_messages(2) = 'File   : ' // trim(GEMSTOOL_INPUTS%MieTmatUser%User_Aerofile_names(nm))
           endif
         endif

!  Exception handling on everything

         if ( Fail2 ) then
            ctype = 'Regular   ' ; if ( do_aer_Jacobians ) ctype = 'Linearized'
            do m = 1, 3
               Messages_Optical(m) = adjustl(trim(MTU_messages(m)))
            enddo
            if ( OpropSrc(nm) .eq. 1 .or. OpropSrc(nm) .eq. 2 ) then
               Messages_Optical(4) = 'Second call to '//ctype//' '//AerChar//' program, mode # '//cnm//', wavelength # '//cwav
            else
               Messages_Optical(4) = 'Second read from '//ctype//' '//AerChar//' file, mode # '//cnm//', wavelength # '//cwav
            endif
            Messages_Optical(5) = 'Failure for '//ctype//', mode # '//cnm
            return
         endif

!   Coefficients Must be copied first. Necessary here....

         if ( do_aer_Jacobians ) then
            if ( OpropSrc(nm) .eq. 1 ) then
               do k = 1, 6
                  MTU_expcoeffs(k,0:MTU_ncoeffs(nm),nm) = Mie_expcoeffs(k,0:MTU_ncoeffs(nm),nm)
                  if ( Do_LinearPSD ) then
                     do q = 1, npsd(nm)
                        LPSD_MTU_expcoeffs(k,0:MTU_ncoeffs(nm),q,nm) = LPSD_Mie_expcoeffs(k,0:MTU_ncoeffs(nm),q,nm)
                     enddo
                  endif
                  if ( Do_LinearRef ) then
                     do q = 1, nlin(nm)
                        LRFE_MTU_expcoeffs(k,0:MTU_ncoeffs(nm),q,nm) = LRFE_Mie_expcoeffs(k,0:MTU_ncoeffs(nm),q,nm)
                     enddo
                  endif
               enddo
            else if ( OpropSrc(nm) .eq. 2 ) then
               do k = 1, 6
                  MTU_expcoeffs(k,0:MTU_ncoeffs(nm)-1,nm) = Tmat_expcoeffs(1:MTU_ncoeffs(nm),k,nm)
                  if ( Do_LinearPSD ) then
                     do q = 1, npsd(nm)
                        LPSD_MTU_expcoeffs(k,0:MTU_ncoeffs(nm)-1,q,nm) = LPSD_Tmat_expcoeffs(1:MTU_ncoeffs(nm),k,q,nm)
                     enddo
                  endif
                  if ( Do_LinearRef ) then
                     do q = 1, nlin(nm)
                        LRFE_MTU_expcoeffs(k,0:MTU_ncoeffs(nm)-1,q,nm) = LRFE_Tmat_expcoeffs(1:MTU_ncoeffs(nm),k,q,nm)
                     enddo
                  endif
               enddo
               MTU_ncoeffs(nm) = MTU_ncoeffs(nm)-1
            else if ( OpropSrc(nm) .eq. 3 ) then
               do k = 1, 6
                  MTU_expcoeffs(k,0:MTU_ncoeffs(nm),nm) = User_expcoeffs(k,0:MTU_ncoeffs(nm),nm)
                  if ( Do_LinearPSD ) then
                     do q = 1, npsd(nm)
                        LPSD_MTU_expcoeffs(k,0:MTU_ncoeffs(nm),q,nm) = LPSD_User_expcoeffs(k,0:MTU_ncoeffs(nm),q,nm)
                     enddo
                  endif
                  if ( Do_LinearRef ) then
                     do q = 1, nlin(nm)
                        LRFE_MTU_expcoeffs(k,0:MTU_ncoeffs(nm),q,nm) = LRFE_User_expcoeffs(k,0:MTU_ncoeffs(nm),q,nm)
                     enddo
                  endif
               enddo
            endif
         else
            if ( OpropSrc(nm) .eq. 1 ) then
               do k = 1, 6
                  MTU_expcoeffs(k,0:MTU_ncoeffs(nm),nm) = Mie_expcoeffs(k,0:MTU_ncoeffs(nm),nm)
               enddo
            else if ( OpropSrc(nm) .eq. 2 ) then
               do k = 1, 6
                  MTU_expcoeffs(k,0:MTU_ncoeffs(nm)-1,nm) = Tmat_expcoeffs(1:MTU_ncoeffs(nm),k,nm)
               enddo
               MTU_ncoeffs(nm) = MTU_ncoeffs(nm)-1
            else if ( OpropSrc(nm) .eq. 3 ) then
               do k = 1, 6
                  MTU_expcoeffs(k,0:MTU_ncoeffs(nm),nm) = User_expcoeffs(k,0:MTU_ncoeffs(nm),nm)
               enddo
            endif
         endif

!  Write aerosol optical properties to file (optional)

         if ( MakeAeroFile(nm) ) then
            if ( do_aer_Jacobians ) then
               call Write_aerosol_oprop_file_plus  &
                 ( max_user_angles, OutFileUnit(nm), w, Micron_wavelength, RefWvl,                                  & ! Inputs
                   maxAerModes, nm, nPLines, DoAeroFilePreamble, AeroFilePreamble,                                  & ! I/O
                   do_LinearRef, do_LinearPSD, NLIN(nm), NPSD(nm), ParsList(1,nm),                                  & ! Inputs
                   MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), MTU_expcoeffs(1,0,nm),           & ! Inputs
                   LRFE_MTU_bulk(1,1,nm), LRFE_MTU_asymm(1,nm), LRFE_MTU_expcoeffs(1,0,1,nm),                       & ! Inputs
                   LPSD_MTU_bulk(1,1,nm), LPSD_MTU_asymm(1,nm), LPSD_MTU_expcoeffs(1,0,1,nm),LPSD_MTU_dist(1,1,nm)  ) ! Inputs
            else
               call Write_aerosol_oprop_file  &
                 ( max_user_angles, OutFileUnit(nm), w, Micron_wavelength, RefWvl,                       & ! Inputs
                   maxAerModes, nm, nPLines, DoAeroFilePreamble, AeroFilePreamble,                       & ! I/O
                   MTU_dist(1,nm), MTU_bulk(1,nm), MTU_asymm(nm), MTU_ncoeffs(nm), MTU_expcoeffs(1,0,nm) ) ! Inputs
            endif
         endif

!  Scattering fraction

         if ( .not. jlee_method ) then
            Csca_frac(nm) = frac(nm)*MTU_bulk(2,nm)
            if ( do_aer_Jacobians ) then
               if ( Do_LinearRef) then
                  LRFE_Csca_frac(1:nlin(nm),nm) = frac(nm) * LRFE_MTU_bulk(2,1:nlin(nm),nm)
               endif
               if ( Do_LinearPSD) then
                  LPSD_Csca_frac(1:npsd(nm),nm) = frac(nm) * LPSD_MTU_bulk(2,1:npsd(nm),nm)
               endif
            endif
         endif

!  End mode loop

      enddo

!  Scattering weights

      if ( .not. jlee_method ) then
         Csca_total  = sum( Csca_frac(1:nmodes) )
         sca_wgt(1:nmodes) = Csca_frac(1:nmodes) / Csca_total 
         qm = 0 ; L_sca_wgt = 0.0_fpk
         do nm = 1, nmodes     
            if ( do_aer_Jacobians ) then
               if ( Do_LinearRef) then
                  q1 = qm+1 ; q2 = qm + nlin(nm) ; qm = qm + nlin(nm)
                  deriv(1:nlin(nm)) = frac(nm) * LRFE_MTU_bulk(2,1:nlin(nm),nm)
                  L_Csca_total(q1:q2)    = deriv(1:nlin(nm)) 
                  L_sca_wgt   (q1:q2,nm) = deriv(1:nlin(nm)) * ( 1.0_fpk - sca_wgt(nm) ) /Csca_total
               endif
               if ( Do_LinearPSD) then
                  q1 = qm+1 ; q2 = qm + npsd(nm) ; qm = qm + npsd(nm)
                  deriv(1:npsd(nm)) = frac(nm) * LPSD_MTU_bulk(2,1:npsd(nm),nm)
                  L_Csca_Total(q1:q2)    = deriv(1:npsd(nm)) 
                  L_sca_wgt   (q1:q2,nm) = frac(nm) * deriv(1:npsd(nm)) * ( 1.0_fpk - sca_wgt(nm) ) /Csca_total
               endif
               if ( nm.eq.2 .and. do_LinearBmf ) then
                  qm = qm + 1 ; L_Csca_Total(qm) = MTU_bulk(2,1) - MTU_bulk(2,2)
                  L_sca_wgt(qm,1:2) = MTU_bulk(2,1:2) * ( 1.0_fpk - sca_wgt(1:2) ) /Csca_total
               endif
            endif
         enddo
!do q = 1, qm1
!  write(*,*)' zzzz', L_sca_wgt(q,1:2), L_Csca_total(q)
!enddo
      endif

!  Set the reference quantities, if reference wavelength is in the list.
!    Values for the first (masked) wavelength
!  Composite extinction @ REF wvl produced by mode extinctions from
!    the Tmat or Mie code weighted by input mode fractions.
!    --> Case: ref wvl IS A MEMBER of the current spectral band

      if ( point_index.ne.0 .and. wc.eq.1 ) then
!mick fix 1/24/2017 - initialize "L_extinction_ref" and add comments as in the
!                     "point_index .eq. 0" if block above
         !extinction_ref = 0.0_fpk ; qm = 0

!  Reference extinction

         extinction_ref   = 0.0_fpk
         L_extinction_ref = 0.0_fpk

!  Jacobian counting

         qm = 0

!  Mode loop

         do nm = 1, nmodes
            MTU_bulk_ext_ref(nm) = MTU_bulk(1,nm)
            extinction_ref = extinction_ref + frac(nm)*MTU_bulk(1,nm)
            if ( do_aer_Jacobians ) then
               if ( Do_LinearRef) then
                  q1 = qm+1 ; q2 = qm + nlin(nm) ; qm = qm + nlin(nm)
                  L_extinction_ref(q1:q2) = L_extinction_ref(q1:q2) + frac(nm) * LRFE_MTU_bulk(1,1:nlin(nm),nm)
                  LRFE_MTU_bulk_ext_ref(1:nlin(nm),nm) = LRFE_MTU_bulk(1,1:nlin(nm),nm)
               endif
               if ( Do_LinearPSD) then
                  q1 = qm+1 ; q2 = qm + npsd(nm) ; qm = qm + npsd(nm)
                  L_extinction_ref(q1:q2) = L_extinction_ref(q1:q2) + frac(nm) * LPSD_MTU_bulk(1,1:npsd(nm),nm)
                  LPSD_MTU_bulk_ext_ref(1:npsd(nm),nm) = LPSD_MTU_bulk(1,1:npsd(nm),nm)
               endif
               if ( nm.eq.2 .and. do_LinearBmf ) then
                  qm = qm + 1 ; L_extinction_ref(qm) = MTU_bulk(1,1) - MTU_bulk(1,2)
               endif
            endif
         enddo

!  End set reference quantities (for case when reference wavelength is in the list)

      endif

!  Scattering linearization (J.Lee Method)
!  ---------------------------------------

      if ( jlee_method ) then

!  Fractions

         do nm = 1, nmodes
            div(nm) = MTU_bulk(2,nm)/MTU_bulk_ext_ref(nm)
            Csca_frac(nm) = frac(nm) * div(nm)
            if ( do_aer_Jacobians ) then
               if ( Do_LinearRef) then
                  deriv(1:nlin(nm)) = LRFE_MTU_bulk(2,1:nlin(nm),nm) - div(nm) * LRFE_MTU_bulk_ext_ref(1:nlin(nm),nm)
                  LRFE_Csca_frac(1:nlin(nm),nm) = frac(nm) * deriv(1:nlin(nm)) / MTU_bulk_ext_ref(nm)
               endif
               if ( Do_LinearPSD) then
                  deriv(1:npsd(nm)) = LPSD_MTU_bulk(2,1:npsd(nm),nm) - div(nm) * LPSD_MTU_bulk_ext_ref(1:npsd(nm),nm)
                  LPSD_Csca_frac(1:npsd(nm),nm) = frac(nm) * deriv(1:npsd(nm)) / MTU_bulk_ext_ref(nm)
               endif
            endif
         enddo

!  total and weights

         Csca_total        = sum( Csca_frac(1:nmodes) )            ! sum of normalized scattering AODs at calc. wavelen.
         sca_wgt(1:nmodes) = Csca_frac(1:nmodes) / Csca_total ! scattering AOD fraction of each mode at calc. wavelen
         qm = 0 ; L_sca_wgt = 0.0_fpk
         do nm = 1, nmodes     
            if ( do_aer_Jacobians ) then
               if ( Do_LinearRef) then
                  q1 = qm+1 ; q2 = qm + nlin(nm) ; qm = qm + nlin(nm)
                  deriv(1:nlin(nm)) = frac(nm) * LRFE_MTU_bulk(2,1:nlin(nm),nm)
                  L_Csca_total(q1:q2)    = deriv(1:nlin(nm)) 
                  L_sca_wgt   (q1:q2,nm) = deriv(1:nlin(nm)) * ( 1.0_fpk - sca_wgt(nm) ) /Csca_total
               endif
               if ( Do_LinearPSD) then
                  q1 = qm+1 ; q2 = qm + npsd(nm) ; qm = qm + npsd(nm)
                  deriv(1:npsd(nm)) = frac(nm) * LPSD_MTU_bulk(2,1:npsd(nm),nm)
                  L_Csca_Total(q1:q2)    = deriv(1:npsd(nm)) 
                  L_sca_wgt   (q1:q2,nm) = frac(nm) * deriv(1:npsd(nm)) * ( 1.0_fpk - sca_wgt(nm) ) /Csca_total
               endif
               if ( nm.eq.2 .and. do_LinearBmf ) then
                  qm = qm + 1 ; L_Csca_Total(qm) = MTU_bulk(2,1) - MTU_bulk(2,2)
                  L_sca_wgt(qm,1:2) = MTU_bulk(2,1:2) * ( 1.0_fpk - sca_wgt(1:2) ) /Csca_total
               endif
            endif
         enddo

!  End J.Lee Method

      endif

!  NOW set the optical property output
!  ===================================

!  Composite extinction @ CURRENT wvl produced by mode extinctions from
!     the Tmat or Mie code weighted by input mode fractions
!     * jlee Method  --> sum of normalized extinction AODs at calculation wavelength

!  Scaling factor. All wavelengths, Using the composite Tmat/Mie extinctions
!     @ CURRENT wvl and REF wvl, define a scale factor to make the resulting
!     aero tau @ CURRENT wvl compatible with the input aero tau @ REF wvl.

!  1a. For the interpolation case, Create output on local arrays
!  -------------------------------------------------------------
 
      if ( interpolate_aerosols ) then

!  Regular properties
!  ******************

!  Extinction and scaling

         extinction = 0.0_fpk
         if ( point_index.ne.0 .and. wc.eq.1 ) then
            extinction = Extinction_ref
            local_aodscaling(w) = 1.0_fpk
            if ( jlee_method ) local_aodscaling(w) =  extinction
         else
            if ( jlee_method ) then
               do nm = 1, nmodes
                  extinction = extinction + frac(nm)*MTU_bulk(1,nm)/MTU_bulk_ext_ref(nm) 
               enddo
               local_aodscaling(w) = extinction
            else
               extinction = dot_product(frac(1:nmodes),MTU_bulk(1,1:nmodes))
               local_aodscaling(w) = extinction / extinction_ref
            endif
         endif
 
!  SSAs and Expansion coefficients

         local_aerssalbs(w) = Csca_total / extinction
         l = 0 ; local_scatmoms(1,0,w) = 1.0d0
         do while ( local_scatmoms(1,l,w).gt.momsize_cutoff .and. l.lt.maxaermoms )
            l = l + 1 ; local_scatmoms(1,l,w) = dot_product(sca_wgt(1:nmodes),MTU_expcoeffs(1,l,1:nmodes))
         enddo
         n_scatmoms_w = l
         do l = 0, n_scatmoms_w
            do k = 2, nmuller
!mick fix 1/24/2017 - replace "1" with "k" in 1st dim of "local_scatmoms"
               !local_scatmoms(1,l,w) = dot_product(sca_wgt(1:nmodes),MTU_expcoeffs(k,l,1:nmodes))
               local_scatmoms(k,l,w) = dot_product(sca_wgt(1:nmodes),MTU_expcoeffs(k,l,1:nmodes))
            enddo
         enddo
         local_nscatmoms(w) = n_scatmoms_w

!  Linearized properties
!  *********************

         if ( do_aer_Jacobians ) then

!  Linearized of extinction, and scaling

            qm = 0
            if ( point_index.ne.0 .and. wc.eq.1 ) then
               L_extinction(1:qm1) = L_extinction_Ref(1:qm1)
               Local_L_aodscaling(w,1:qm1) = 0.0_fpk
               if ( jlee_method )  Local_L_aodscaling(w,1:qm1) =  L_extinction(1:qm1)
            else
               if ( jlee_method ) then
                  do nm = 1, nmodes
                     div(nm) = MTU_bulk(1,nm)/MTU_bulk_ext_ref(nm)
                     if ( Do_LinearRef) then
                        q1 = qm+1 ; q2 = qm + nlin(nm) ; qm = qm + nlin(nm)
                        deriv(1:nlin(nm)) = LRFE_MTU_bulk(1,1:nlin(nm),nm) - div(nm) * LRFE_MTU_bulk_ext_ref(1:nlin(nm),nm) 
                        L_extinction(q1:q2)    = frac(nm) * deriv(1:nlin(nm)) / MTU_bulk_ext_ref(nm)
                     endif
                     if ( Do_LinearPSD) then
                        q1 = qm+1 ; q2 = qm + npsd(nm) ; qm = qm + npsd(nm)
                        deriv(1:npsd(nm)) = LPSD_MTU_bulk(2,1:npsd(nm),nm) - div(nm) * LPSD_MTU_bulk_ext_ref(1:npsd(nm),nm) 
                        L_extinction(q1:q2)    = frac(nm) * deriv(1:npsd(nm)) / MTU_bulk_ext_ref(nm)
                     endif
                     if ( nm.eq.2 .and. do_LinearBmf ) then
                        qm = qm + 1 ; L_extinction(qm) = div(1) - div(2)
                     endif
                  enddo
                  Local_L_aodscaling(w,1:qm1) = L_extinction(1:qm1)
               else
                  do nm = 1, nmodes
                     if ( Do_LinearRef) then
                        q1 = qm+1 ; q2 = qm + nlin(nm) ; qm = qm + nlin(nm)
                        L_extinction(q1:q2)    = frac(nm) * LRFE_MTU_bulk(1,1:nlin(nm),nm)
                     endif
                     if ( Do_LinearPSD) then
                        q1 = qm+1 ; q2 = qm + npsd(nm) ; qm = qm + npsd(nm)
                        L_extinction(q1:q2)    = frac(nm) * LPSD_MTU_bulk(1,1:npsd(nm),nm)
                     endif
                     if ( nm.eq.2 .and. do_LinearBmf ) then
                        qm = qm + 1 ; L_extinction(qm) = MTU_bulk(1,1) - MTU_bulk(1,2)
                     endif
                  enddo
                  do q = 1, qm1
                     Local_L_aodscaling(w,q) = ( L_extinction(q) - Local_aodscaling(w) * L_extinction_ref(q) ) / extinction_ref 
                  enddo
               endif
            endif

!  Linearization of scattering properties
!  SSAs and Expansion coefficients

            do q = 1, qm1
               Local_L_aerssalbs(w,q) = ( L_Csca_total(q) - Local_aerssalbs(w) * L_extinction(q) ) / extinction
            enddo

!  Bug 12/20/16. Initializing Jacobian count must be done inside K-loop.
!            qm = 0

            do k = 1, nmuller
               qm = 0              !  This is correct. 12/20/16
               do nm = 1, nmodes
                  if ( Do_LinearRef) then
                     q1 = qm+1 ; q2 = qm + nlin(nm) ; qm = qm + nlin(nm)
                     do l = 0, n_scatmoms_w
                        Local_L_scatmoms(k,l,w,q1:q2) = sca_wgt(nm)       * LRFE_MTU_expcoeffs(k,l,1:nlin(nm),nm) &
                                                    + L_sca_wgt(q1:q2,nm) *      MTU_expcoeffs(k,l,nm)
                     enddo
                  endif
                  if ( Do_LinearPSD) then
                     q1 = qm+1 ; q2 = qm + npsd(nm) ; qm = qm + npsd(nm)
                     do l = 0, n_scatmoms_w
                        Local_L_scatmoms(k,l,w,q1:q2) = sca_wgt(nm)       * LPSD_MTU_expcoeffs(k,l,1:npsd(nm),nm) &
                                                    + L_sca_wgt(q1:q2,nm) *      MTU_expcoeffs(k,l,nm)
                     enddo
                  endif
                  if ( nm.eq.2 .and. do_LinearBmf ) then
                     qm = qm + 1
                     do l = 0, n_scatmoms_w
                        Local_L_scatmoms(k,l,w,qm) = &
                                        L_sca_wgt(qm,1) *  MTU_expcoeffs(k,l,1) - L_sca_wgt(qm,2) *  MTU_expcoeffs(k,l,2)
                     enddo
                  endif
               enddo
            enddo

!  End Jacobians clause

         endif

!  End "interpolate aerosols" if block

      endif

!  1b. For the Monochromatic case, Create output on final arrays
!  -------------------------------------------------------------

      if ( .not. interpolate_aerosols ) then

!  Regular properties
!  ******************

!  Extinction and scaling

         extinction = 0.0_fpk
         if ( point_index.ne.0 .and. wc.eq.1 ) then
            extinction = Extinction_ref
            aod_scaling(w) = 1.0_fpk
            if ( jlee_method ) aod_scaling(w) = extinction
         else
            if ( jlee_method ) then
               do nm = 1, nmodes
                  extinction = extinction + frac(nm)*MTU_bulk(1,nm)/MTU_bulk_ext_ref(nm) 
               enddo
               aod_scaling(w) = extinction
            else
               extinction = dot_product(frac(1:nmodes),MTU_bulk(1,1:nmodes))
               aod_scaling(w) = extinction / extinction_ref
            endif
         endif
 
!  SSAs and assymetry parameters (Not doing the latter)
!   Rob Fix 8/25/14. Total SSA is not scatter-weighted sum of individual mode SSAs

         aerosol_ssalbs(w) = Csca_total / extinction
!         aerosol_Asymms(w) = 0.0_fpk
!         do nm = 1, nmodes
!            aerosol_Asymms(w)   = aerosol_Asymms(w)  + sca_wgt(nm)*MTU_Asymm(nm)
!         enddo

         l = 0 ;  aerosol_scatmoms(1,0,w) = 1.0_fpk
         do while ( aerosol_scatmoms(1,l,w).gt.momsize_cutoff .and. l.lt.maxaermoms )
            l = l + 1
            do nm = 1, nmodes
               aerosol_scatmoms(1,l,w) = aerosol_scatmoms(1,l,w) + sca_wgt(nm)*MTU_expcoeffs(1,l,nm)
            enddo
         enddo
         n_scatmoms_w = l
         do l = 0, n_scatmoms_w
            do k = 2, nmuller
               do nm = 1, nmodes
                  aerosol_scatmoms(k,l,w) = aerosol_scatmoms(k,l,w) + sca_wgt(nm)*MTU_expcoeffs(k,l,nm)
               enddo
            enddo
         enddo
         n_scatmoms(w) = n_scatmoms_w

!  Linearized properties
!  *********************

         if ( do_aer_Jacobians ) then

!  Linearized extinction, and scaling

            if ( point_index .ne. 0 .and. wc.eq.1 ) then
               L_extinction(1:qm1) = L_extinction_Ref(1:qm1)
               L_aod_scaling(w,1:qm1) = 0.0_fpk
               if ( jlee_method )  L_aod_scaling(w,1:qm1) =  L_extinction(1:qm1)
            else
               if ( jlee_method ) then
                  do nm = 1, nmodes
                     div(nm) = MTU_bulk(1,nm)/MTU_bulk_ext_ref(nm)
                     if ( Do_LinearRef) then
                        q1 = qm+1 ; q2 = qm + nlin(nm) ; qm = qm + nlin(nm)
                        deriv(1:nlin(nm)) = LRFE_MTU_bulk(1,1:nlin(nm),nm) - div(nm) * LRFE_MTU_bulk_ext_ref(1:nlin(nm),nm) 
                        L_extinction(qm1:qm2)    = frac(nm) * deriv(1:nlin(nm)) / MTU_bulk_ext_ref(nm)
                     endif
                     if ( Do_LinearPSD) then
                        q1 = qm+1 ; q2 = qm + npsd(nm) ; qm = qm + npsd(nm)
                        deriv(1:npsd(nm)) = LPSD_MTU_bulk(2,1:npsd(nm),nm) - div(nm) * LPSD_MTU_bulk_ext_ref(1:npsd(nm),nm) 
                        L_extinction(qm1:qm2)    = frac(nm) * deriv(1:npsd(nm)) / MTU_bulk_ext_ref(nm)
                     endif
                     if ( nm.eq.2 .and. do_LinearBmf ) then
                        qm = qm + 1 ; L_extinction(qm) = div(1) - div(2)
                     endif
                  enddo
                  L_aod_scaling(w,1:qm1) = L_extinction(1:qm1)
               else
                  do nm = 1, nmodes
                     if ( Do_LinearRef) then
                        q1 = qm+1 ; q2 = qm + nlin(nm) ; qm = qm + nlin(nm)
                        L_extinction(qm1:qm2)    = frac(nm) * LRFE_MTU_bulk(1,1:nlin(nm),nm)
                     endif
                     if ( Do_LinearPSD) then
                        q1 = qm+1 ; q2 = qm + npsd(nm) ; qm = qm + npsd(nm)
                        L_extinction(qm1:qm2)    = frac(nm) * LPSD_MTU_bulk(1,1:npsd(nm),nm)
                     endif
                     if ( nm.eq.2 .and. do_LinearBmf ) then
                        qm = qm + 1 ; L_extinction(qm) = MTU_bulk(1,1) - MTU_bulk(1,2)
                     endif
                  enddo
                  do q = 1, qm1
                     L_aod_scaling(w,q) = ( L_extinction(q) - aod_scaling(w) * L_extinction_ref(q) ) / extinction_ref 
                  enddo
               endif
            endif

!  Linearization of scattering properties
!  SSAs and Expansion coefficients

            do q = 1, qm1
               L_aerosol_ssalbs(w,q) = ( L_Csca_total(q) - aerosol_ssalbs(w) * L_extinction(q) ) / extinction
            enddo

!  Bug 12/20/16. Jacobian count not initialized at all, for the K-loop

            do k = 1, nmuller
               qm = 0      ! 12/20/16 This inserted here.
               do nm = 1, nmodes
                  if ( Do_LinearRef) then
                     q1 = qm+1 ; q2 = qm + nlin(nm) ; qm = qm + nlin(nm)
                     do l = 0, n_scatmoms_w
                        L_aerosol_scatmoms(k,l,w,q1:q2) = sca_wgt(nm)       * LRFE_MTU_expcoeffs(k,l,1:nlin(nm),nm) &
                                                      + L_sca_wgt(q1:q2,nm) *      MTU_expcoeffs(k,l,nm)
                     enddo
                  endif
                  if ( Do_LinearPSD) then
                     q1 = qm+1 ; q2 = qm + npsd(nm) ; qm = qm + npsd(nm)
                     do l = 0, n_scatmoms_w
                        L_aerosol_scatmoms(k,l,w,q1:q2) = sca_wgt(nm)       * LPSD_MTU_expcoeffs(k,l,1:npsd(nm),nm) &
                                                      + L_sca_wgt(q1:q2,nm) *      MTU_expcoeffs(k,l,nm)
                     enddo
                  endif
                  if ( nm.eq.2 .and. do_LinearBmf ) then
                     qm = qm + 1
                     do l = 0, n_scatmoms_w
                        L_aerosol_scatmoms(k,l,w,qm) = &
                                  L_sca_wgt(qm,1) *  MTU_expcoeffs(k,l,1) - L_sca_wgt(qm,2) *  MTU_expcoeffs(k,l,2)
                     enddo
                  endif
               enddo
            enddo

!  End Jacobians clause

         endif

!  Apply scalings to loadings and linearizations
!     (i.e. define the layer aero tau @ CURRENT wvl so it's
!           compatible with the layer aero tau @ REF wvl using the
!           previously defined scale factors)

         do n = 1, nlayers
            aerosol_deltau(n,w) = aertau_unscaled(n) * aod_scaling(w)
         enddo
         if ( do_aer_jacobians ) then
            do n = 1, nlayers
               do q = 1, qmt
                  L_aerosol_deltau(n,w,q) = L_aertau_unscaled(n,q) *   aod_scaling(w) &
                                            + aertau_unscaled(n)   * L_aod_scaling(w,q)
               enddo
            enddo
         endif

!  Debug aerosol optical properties. VERSION TWO only

         if ( do_iopchecker ) then
            do n = 1, nlayers
              if (aerlayerflags(N).and.n.eq.107 ) then
                 write(999,'(i4,1p6e20.10)')n,aerosol_deltau(n,w),aerosol_ssalbs(w)
                 do l = 0, n_scatmoms_w
                    write(999,'(2i5,1p6e20.10)')n,l,(aerosol_scatmoms(k,l,w),k=1,1)
                 enddo
              endif
            enddo
         endif

!  Debug aerosol optical properties. VERSION TWO only

         if ( do_iopchecker ) then
          if ( do_aer_jacobians ) then
           do n = 1, nlayers
            if (aerlayerflags(n).and.n.eq.107 ) then
              write(777,'(i4,1p6e20.10)')n,aerosol_deltau(n,w),aerosol_ssalbs(w)
              write(777,'(i4,1p11e20.10)')n,(L_aerosol_deltau(n,w,q),q=1,11)
              write(777,'(i4,1p11e20.10)')n,(L_aerosol_ssalbs(w,q),q=1,11)
              do l = 0, n_scatmoms_w
                write(777,'(2i5,1pe20.10,1p11e15.6)')n,l,aerosol_scatmoms(1,l,w),(l_aerosol_scatmoms(1,l,w,q),q=1,11)
              enddo
            endif
           enddo
!           stop'Lin 777'  ! pause'Lin 777'
          else
           do n = 1, nlayers
            if (aerlayerflags(N).and.n.eq.107 ) then
              write(888,'(i4,1p6e20.10)')n,aerosol_deltau(n,w),aerosol_ssalbs(w)
              do l = 0, n_scatmoms_w
                write(888,'(2i5,1p6e20.10)')n,l,(aerosol_scatmoms(k,l,w),k=1,1)
              enddo
            else
              write(888,'(i4,1p6e20.10)')n,aerosol_deltau(n,w)
            endif
           enddo
!           stop'Reg 888'  ! pause'Reg 888'
          endif
         endif

! End monochromatic ("not interpolate aerosols") if block
!  (replaces continuation point 677)

      endif

!  continuation point for avoiding the exact monochromatic solution
!677   continue

!  Progress

      if ( Verbose.or.Debug ) then
         write(*,'(4x,a,1x,i3,1x,f12.5,2(1x,i3))') &
            'Aerosol calculation - Finished point, wavelength   : ',wc, wavelength
      endif

!  End wavelength/wavenumber loop (note: this loop starts @ ~line 1640)

   enddo

!  1c. carry out interpolation to get final results
!  ================================================

   if ( interpolate_aerosols ) then

!  Interpolation, wavelength regime

     if ( .not. do_wavnums ) then
       wastart = 1

!  Start wavelength loop

       do w = 1, nwav

!  Define weights to interpolate from aerosol grid to RT grid

         wa = wastart ; trawl = .true.
         do while (trawl)
           if ( wav(w).ge.aerwav(wa) .and. wav(w).le.aerwav(wa+1) ) trawl = .false.
         enddo
         wa1 = wa ; wa2 = wa + 1
         fa1 = ( aerwav(wa2) - wav(w) ) / ( aerwav(wa2) -  aerwav(wa1) ) ; fa2 = 1.0d0 - fa1
         wastart = wa1
         if ( w.lt.nwav) then
           if( wav(w+1).ge.aerwav(wa+1)) wastart = wa2
         endif

!  Interpolate regular quantities

         aod_scaling(w) = fa1 * local_aodscaling(wa1) + fa2 * local_aodscaling(wa2)
         do n = 1, nlayers
           aerosol_deltau(n,w) = aertau_unscaled(n) * aod_scaling(w) 
         enddo
         aerosol_ssalbs(w) = fa1 * local_aerssalbs(wa1) + fa2 * local_aerssalbs(wa2)
         n_scatmoms(w) = min(local_nscatmoms(wa1),local_nscatmoms(wa2))
         do l = 0, n_scatmoms(w)
!mick mod 1/24/2017 - turn off do loop since f90 vector ops already in use here
           !do k = 1, nmuller
             aerosol_scatmoms(1:nmuller,l,w) = fa1 * local_scatmoms(1:nmuller,l,wa1) &
                                             + fa2 * local_scatmoms(1:nmuller,l,wa2)
           !enddo
         enddo
!         aerosol_asymms(w) = aerosol_scatmoms(1,1,w) / 3.0d0
         aerosol_scatmoms(1,0,w) = 1.0d0

!  Interpolate linearized quantities

         if ( do_aer_Jacobians ) then
!mick fix 1/24/2017 - correct "L_aod_scaling" for elements wrt aero load pars
           !L_aod_scaling(w,  1:qmt) = fa1 * local_L_aodscaling(wa1,1:qmt) &
           !                         + fa2 * local_L_aodscaling(wa2,1:qmt)
           L_aod_scaling(w,  1:qm1) = fa1 * local_L_aodscaling(wa1,1:qm1) &
                                    + fa2 * local_L_aodscaling(wa2,1:qm1)
           L_aod_scaling(w,qm2:qmt) = 0.0d0

           do n = 1, nlayers
             L_aerosol_deltau(n,w,1:qmt) =   aertau_unscaled(n)       * L_aod_scaling(w,1:qmt) &
                                         + L_aertau_unscaled(n,1:qmt) *   aod_scaling(w)
           enddo
           L_aerosol_ssalbs(w,1:qm1) = fa1 * local_L_aerssalbs(wa1,1:qm1) + fa2 * local_L_aerssalbs(wa2,1:qm1)
           do l = 0, n_scatmoms(w)
             do k = 1, nmuller
               L_aerosol_scatmoms(k,l,w,1:qm1) = fa1 * local_L_scatmoms(k,l,wa1,1:qm1) &
                                               + fa2 * local_L_scatmoms(k,l,wa2,1:qm1)
             enddo
           enddo
           L_aerosol_scatmoms(1,0,w,1:qm1) = 0.0d0
         endif

!  End wavelength loop

       enddo

!  End wavelength regime

     endif

!  Interpolation, wavenumber regime
!   12 August 2013 --> KLUTZY CODE here,,,,,,,,,,Improve it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if ( do_wavnums ) then
       wastart = 1

!  Start wavenumber loop

       do w = 1, nwav

!  Define weights to interpolate from aerosol grid to RT grid

         wa = wastart ; trawl = .true. ; lam =  1.0d+07/wav(w)
54       continue
         do while (trawl)
           if ( lam .lt. aerwav(wa) ) then
             if ( lam .ge. aerwav(wa+1) ) then
               trawl = .false.
             else if (lam .lt. aerwav(wa+1) ) then
               wa = wa + 1 ; go to 54
             endif
           endif
         enddo
         wa1 = wa ; wa2 = wa + 1
         fa1 = ( aerwav(wa2) - lam ) / ( aerwav(wa2) -  aerwav(wa1) ) ; fa2 = 1.0d0 - fa1
         wastart = wa1

!  Interpolate regular quantities

         aod_scaling(w) = fa1 * local_aodscaling(wa1) + fa2 * local_aodscaling(wa2)
         do n = 1, nlayers
           aerosol_deltau(n,w) = aertau_unscaled(n) * aod_scaling(w)
         enddo
         aerosol_ssalbs(w) = fa1 * local_aerssalbs(wa1) + fa2 * local_aerssalbs(wa2)
         n_scatmoms(w) = min(local_nscatmoms(wa1),local_nscatmoms(wa2))
         do l = 0, n_scatmoms(w)
!mick fix 1/24/2017 - turn off do loop since f90 vector ops already in use here
           !do k = 1, nmuller
             aerosol_scatmoms(1:nmuller,l,w) = fa1 * local_scatmoms(1:nmuller,l,wa1) &
                                             + fa2 * local_scatmoms(1:nmuller,l,wa2)
           !enddo
         enddo
!         aerosol_asymms(w) = aerosol_scatmoms(1,1,w) / 3.0d0
         aerosol_scatmoms(1,0,w) = 1.0d0

!  Interpolate linearized quantities

         if ( do_aer_Jacobians ) then
!mick fix 1/24/2017 - correct "L_aod_scaling" for elements wrt aero load pars
           !L_aod_scaling(w,  1:qmt) = fa1 * local_L_aodscaling(wa1,1:qmt) &
           !                         + fa2 * local_L_aodscaling(wa2,1:qmt)
           L_aod_scaling(w,  1:qm1) = fa1 * local_L_aodscaling(wa1,1:qm1) &
                                    + fa2 * local_L_aodscaling(wa2,1:qm1)
           L_aod_scaling(w,qm2:qmt) = 0.0d0

           do n = 1, nlayers
             L_aerosol_deltau(n,w,1:qmt) =   aertau_unscaled(n)       * L_aod_scaling(w,1:qmt) &
                                         + L_aertau_unscaled(n,1:qmt) *   aod_scaling(w)
           enddo
           L_aerosol_ssalbs(w,1:qm1) = fa1 * local_L_aerssalbs(wa1,1:qm1) + fa2 * local_L_aerssalbs(wa2,1:qm1)
           do l = 0, n_scatmoms(w)
             do k = 1, nmuller
               L_aerosol_scatmoms(k,l,w,1:qm1) = fa1 * local_L_scatmoms(k,l,wa1,1:qm1) &
                                               + fa2 * local_L_scatmoms(k,l,wa2,1:qm1)
             enddo
           enddo
           L_aerosol_scatmoms(1,0,w,1:qm1) = 0.0d0
         endif

!  Debug aerosol optical properties. VERSION TWO only

         if ( do_iopchecker ) then
           if ( do_aer_jacobians ) then
              write(777,'(i4,1p6e20.10)')w,aerosol_deltau(nlayers,w),aerosol_ssalbs(w)
              write(777,'(i4,1p11e20.10)')w,(L_aerosol_deltau(nlayers,w,q),q=1,11)
              write(777,'(i4,1p11e20.10)')w,(L_aerosol_ssalbs(w,q),q=1,11)
              do l = 0, n_scatmoms_w
                write(777,'(2i5,1pe20.10,1p11e15.6)')w,l,aerosol_scatmoms(1,l,w),(l_aerosol_scatmoms(1,l,w,q),q=1,11)
              enddo
!              if (w.eq.10)stop'Lin 777'  ! pause'Lin 777'
          else
             write(888,'(i4,1p6e20.10)')w,aerosol_deltau(nlayers,w),aerosol_ssalbs(w)
             do l = 0, n_scatmoms_w
               write(888,'(2i5,1p6e20.10)')nlayers,l,(aerosol_scatmoms(k,l,w),k=1,1)
             enddo
!             if ( w.eq.10)stop'Reg 888'  ! pause'Reg 888'
          endif
         endif

!  End wavenumber loop

       enddo

!  End wavenumber regime

    endif

!  Finish aerosol interpolation 

   endif

!  Safety - File Closures, before leaving routine
!  (note: these files opened @ ~line 1320)

   if ( GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) then
      do nm = 1, nmodes
         close(InFileUnit(nm))
      enddo
   else
      do nm = 1, nmodes
         if ( MakeAeroFile(nm) ) close(OutFileUnit(nm))
      enddo
   endif

!  End subroutine

   return

end subroutine GEMSTOOL_AER_PROPERTIES_PLUS

!  End module

end Module GEMSTOOL_AerProperties_Plus_m

