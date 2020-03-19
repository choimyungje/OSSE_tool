module RTSMie_master_bimodal_m

!  This is the Bimodal master #1 for RTS MIE code
!    ** RT Solutions, Version 1.4, 30 June     2011 (Bimodal control Tmatrix)
!    ** RT Solutions, Version 1.5, 17 August   2011 (Bimodal control Mie)

  use RTSMie_parameters_m
  use RTSMie_sourcecode_m

!  Everything PUBLIC here
!  ----------------------

public

contains

subroutine RTSMie_master_bimodal                             &
       ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
         PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,   & ! I
         nblocks, nweights, xparticle_limit, R1R2_cutoff,    & ! I
         n_Fmatrix_angles, Fmatrix_angles,                   & ! I
         lambda, n_real, n_imag, fraction,                   & ! I
         BMie_bulk, BMie_asymm, BMie_ncoeffs,                & ! O
         BMie_expcoeffs, BMie_Fmatrix, BMie_dist,            & ! O
         fail, istatus, Bmessages, trace_3 )                   ! O

   USE RTSMie_parameters_m

!  implicit none statement

   IMPLICIT NONE

!  List of Inputs
!  ==============

!  Flag inputs
!  -----------

!      Do_Expcoeffs      - Boolean flag for computing Expansion Coefficients
!      Do_Fmatrix        - Boolean flag for computing F-matrix at equal-angles

   logical  , intent(in)  :: Do_Expcoeffs
   logical  , intent(in)  :: Do_Fmatrix

!      Do_Monodisperse   - Boolean flag for Doing a Monodisperse calculation
!                          If set, the PSD stuff will be turned off internally

   LOGICAL  , INTENT (IN) :: do_Monodisperse

!  PSD inputs
!  ----------

!  FixR1R2 : If  set, Use "R1R2_cutoff" for smallest particle size
!                     then Internal routine to calculate R1 and R2 (outputs)
!            If Not set, Use Input R1 and R2 for PSD limits.

   logical, intent(inout)  :: FixR1R2(2)

!  R1, R2         - Minimum and Maximum radii (Microns)

   real    (KIND=dp), intent(inout)  :: R1(2), R2(2)

!  Limiting particle size value. Set to 10000.0 default.
!   If you exceed this, program will tell you to increase dimensioning.

   REAL    (KIND=dp), INTENT (IN) :: xparticle_limit

!      Monoradius     - Monodisperse radius size (Microns)

   real    (KIND=dp), intent(in)  :: Monoradius

!      psd_Index      - Index for particle size distribution of spheres
!      psd_pars       - Parameters characterizing PSD (up to 3 allowed)

!  Mie inputs (distribution index, PSD parameters)
!    PSD_index = 1 : TWO-PARAMETER GAMMA with alpha and b given
!    PSD_index = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given
!    PSD_index = 3 : BIMODAL GAMMA with equal mode weights
!    PSD_index = 4 : LOG-NORMAL with rg and sigma given
!    PSD_index = 5 : LOG-NORMAL with reff and veff given
!    PSD_index = 6 : POWER LAW
!    PSD_index = 7 : MODIFIED GAMMA with alpha, rc and gamma given
!    PSD_index = 8 : MODIFIED GAMMA with alpha, b and gamma given

   integer          , intent(in)  :: PSD_Index(2)
   real    (KIND=dp), intent(in)  :: PSD_pars (3,2)

!  PSD quadrature control
!  ----------------------

!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.
!    R1R2_cutoff particle size for setting R1 and R2 internally

   INTEGER          , INTENT (IN) :: nblocks(2)
   INTEGER          , INTENT (IN) :: nweights(2)
   REAL    (KIND=dp), INTENT (IN) :: R1R2_cutoff(2)

!  Optical: Wavelength, refractive index
!  -------------------------------------

!      LAMBDA         - wavelength of light (microns)
!      N_REAL, N_IMAG - real and imaginary parts, refractive index (N-i.GE.0)

   real    (KIND=dp), intent(in)  :: lambda, n_real(2), n_imag(2)

!  Fraction (by volume weight, First PSD mode)

   real    (KIND=dp), intent(in)  :: fraction

!  F-matrix Angular control input
!  ------------------------------

!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER          , INTENT (IN) :: n_Fmatrix_angles
   REAL    (KIND=dp), INTENT (IN) :: Fmatrix_angles(max_Mie_angles)

!  Output arguments
!  ================

!  Bulk distribution parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(KIND=dp), intent(out) :: BMie_bulk (3)

!  Expansion coefficients and Asymmetry parameter, optional output

   integer      , intent(out) :: BMie_ncoeffs
   real(KIND=dp), intent(out) :: BMie_expcoeffs (6,0:max_Mie_angles)
   real(KIND=dp), intent(out) :: BMie_asymm

!  F-matrix, optional output

   real(KIND=dp), intent(out) :: BMie_Fmatrix(4,max_Mie_angles)

!  Distribution parameters. BOTH MODES !!!!!!!!!!
!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(KIND=dp), intent(out) :: BMie_dist (5,2)

!  Exception handling

   LOGICAL          , INTENT (OUT)   :: fail
   INTEGER          , INTENT (OUT)   :: istatus
   CHARACTER*(*)    , INTENT (OUT)   :: Bmessages(3), trace_3

!  Local Arrays
!  ------------

!  Bulk distribution parameters
!  Expansion coefficients and Asymmetry parameter
!  F-matrix,  optional output

   real(KIND=dp) :: Mie1_bulk (3)
   integer       :: Mie1_ncoeffs
   real(KIND=dp) :: Mie1_expcoeffs (6,0:max_Mie_angles)
   real(KIND=dp) :: Mie1_asymm
   real(KIND=dp) :: Mie1_Fmatrix (4,max_Mie_angles)

   real(KIND=dp) :: Mie2_bulk (3)
   integer       :: Mie2_ncoeffs
   real(KIND=dp) :: Mie2_expcoeffs (6,0:max_Mie_angles)
   real(KIND=dp) :: Mie2_asymm
   real(KIND=dp) :: Mie2_Fmatrix (4,max_Mie_angles)

!  Other local variables
!  ---------------------

   integer       :: k, L
   real(KIND=dp) :: FF1, FF2, WW1, WW2
   real(KIND=dp) :: Csca1_FF1, Csca2_FF2, Csca_total

!  Zero the output
!  ---------------

   BMie_bulk      = d_zero
   BMie_Fmatrix   = d_zero
   BMie_expcoeffs = d_zero
   BMie_asymm     = d_zero
   BMie_ncoeffs   = 0
   BMie_dist      = d_zero

   FF1 = fraction
   FF2 = d_one - FF1

   fail           = .false.
   istatus        = 0
   Bmessages(1:3) = ' '
   trace_3        = ' '

!  Check: No Monodisperse here !
!  -----------------------------

   if ( do_monodisperse ) then
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_bimodal module: Input error: MONODISPERSE FLAG must be Off!'
      return
   endif

!  First Call
!  ---------

   k = 1 ; write(*,*)' ** Doing Mie Regular for PSD # 1 ----------------------'
   CALL RTSMie_main                                          & 
       ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
         PSD_Index(k), PSD_pars(:,k), MonoRadius,            & ! I
         R1(k), R2(k), FixR1R2(k), nblocks(k), nweights(k),  & ! I
         xparticle_limit, R1R2_cutoff(k),                    & ! I
         n_Fmatrix_angles, Fmatrix_angles,                   & ! I
         lambda, n_real(k), n_imag(k),                       & ! I
         Mie1_bulk, Mie1_asymm, Mie1_ncoeffs,                & ! O
         Mie1_expcoeffs, Mie1_Fmatrix, BMie_dist(:,k),       & ! O
         fail, istatus, Bmessages(1), Bmessages(2), Bmessages(3) ) ! O

!  Exception handling

   if ( fail ) then
      trace_3 = 'RTSMie_master_bimodal module: First PSD call, Error'
      if ( Istatus .eq. 1 ) return
   endif

!  Second call
!  -----------

   k = 2 ; write(*,*)' ** Doing Mie Regular for PSD # 2 ----------------------'
   CALL RTSMie_main                                          & 
       ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
         PSD_Index(k), PSD_pars(:,k), MonoRadius,            & ! I
         R1(k), R2(k), FixR1R2(k), nblocks(k), nweights(k),  & ! I
         xparticle_limit, R1R2_cutoff(k),                    & ! I
         n_Fmatrix_angles, Fmatrix_angles,                   & ! I
         lambda, n_real(k), n_imag(k),                       & ! I
         Mie2_bulk, Mie2_asymm, Mie2_ncoeffs,                & ! O
         Mie2_expcoeffs, Mie2_Fmatrix, BMie_dist(:,k),       & ! O
         fail, istatus, Bmessages(1), Bmessages(2), Bmessages(3) ) ! O

!  Exception handling

   if ( fail ) then
      trace_3 = 'RTSMie_master_bimodal module: Second PSD call, Error'
      if ( Istatus .eq. 1 ) return
   endif

!  Bimodal determination
!  ---------------------

!  Revision 20 September 2011
!    Correct definition for Expcoeffs/Fmatrix: WW1/WW2 in place of FF1/FF2

   Csca1_FF1  = FF1 * Mie1_bulk(2)   
   Csca2_FF2  = FF2 * Mie2_bulk(2)   
   Csca_total =  Csca1_FF1 + Csca2_FF2
   WW1   = Csca1_FF1 / Csca_total
   WW2   = Csca2_FF2 / Csca_total

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong
   BMie_bulk(1:2) = FF1 * Mie1_bulk(1:2) + FF2 * Mie2_bulk(1:2)
   BMie_bulk(3)   = BMie_bulk(2) / BMie_bulk(1)

!      write(18,68) 'Bulk-B-pert    ',(BMie_bulk(k),k=1,3)
!68 format(a15,1p3e20.10)

!  original code
!   BMie_bulk(1:3) = FF1 * Mie1_bulk(1:3) + FF2 * Mie2_bulk(1:3)

   if ( Do_Expcoeffs ) then
      BMie_asymm = WW1 * Mie1_asymm + WW2 * Mie2_asymm
      BMie_ncoeffs = max(Mie1_ncoeffs,Mie2_ncoeffs)      
      do L = 0, min(Mie1_ncoeffs,Mie2_ncoeffs)
         BMie_expcoeffs(1:6,L) = WW1 * Mie1_expcoeffs(1:6,L) + & 
                                 WW2 * Mie2_expcoeffs(1:6,L)
      enddo
      if ( Mie1_ncoeffs .lt. Mie2_ncoeffs ) then
          do L = Mie1_ncoeffs + 1,Mie2_ncoeffs
              BMie_expcoeffs(1:6,L) = WW2 * Mie2_expcoeffs(1:6,L)
          enddo
      else if ( Mie1_ncoeffs .gt. Mie2_ncoeffs ) then
          do L = Mie2_ncoeffs + 1,Mie1_ncoeffs
              BMie_expcoeffs(1:6,L) = WW1 * Mie1_expcoeffs(1:6,L)
          enddo
      endif
   endif
   if ( Do_Fmatrix ) then
      do L = 1, n_Fmatrix_angles
         BMie_Fmatrix(1:4,L) = WW1 * Mie1_Fmatrix(1:4,L) + & 
                               WW2 * Mie2_Fmatrix(1:4,L)
      enddo
   endif   

!  Finish

   return
end subroutine RTSMie_master_bimodal

!  End module

end module RTSMie_master_bimodal_m
