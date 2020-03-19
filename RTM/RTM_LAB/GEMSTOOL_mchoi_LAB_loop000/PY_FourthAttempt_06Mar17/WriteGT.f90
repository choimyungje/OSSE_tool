module Write_TriAxaerosol_file_m

!  This is the module that writes to GEMSTOOL User-aerosol file
!  Based on subroutine "Write_aerosol_oprop_file".

public

contains

subroutine Write_TriAxaerosol_file ( &
   max_Coeffs, OutFileUnit, w, Micron_wavelength, RefWvl, &
   TriAx_dist, TriAx_bulk, TriAx_asymm, TriAx_ncoeffs, TriAx_expcoeffs )

!  Introduced 10/19/16. Adapted from CODECTOOL routine of the same name.
!    Major new thing: aerosol distribution characteristics for Ref Wavl.

!  This routine created 12/28/16 for writing TriAx aerosols to file.
!  This routine amended 03/03/17 for Mk3 attempt. No file-name

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   integer, intent(in)       :: max_coeffs
   integer, intent(in)       :: OutFileUnit, w
   real(fpk), intent(in)     :: Micron_wavelength
   logical, intent(in)       :: RefWvl

!  Distribution parameters

   real(fpk), intent(in) :: TriAx_dist(5)

!  Bulk optical parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk), intent(in)  :: TriAx_bulk(3)

!  Asymmetry parameter and number of expansion coefficients

   real(fpk), intent(in)  :: TriAx_asymm

!  Number and value of expansion coefficients

   integer, intent(in)    :: TriAx_ncoeffs
   real(fpk), intent(in)  :: TriAx_expcoeffs(6,0:max_coeffs)

!  Local variables
!  ---------------

   integer :: i,l


!  Write output to file

   write(OutFileUnit,'(3x,a)') '--------------------------------------------' // &
                               '--------------------------------------------' // &
                               '--------------------------------------------'
   write(OutFileUnit,'(6x,a,7x,a)') 'Wvl idx','Wvl (um)'
   if ( .not.RefWvl ) then
      write(OutFileUnit,'(8x,i4.4,e20.11)') w,Micron_wavelength

      !Sample:
      !         Ext Coeff           Scat Coeff              SSA               Asym Par       Num of Expan Coeffs
      !      0.30186852870E-01   0.29730300500E-01   0.98487578775E+00   0.59090185717E+00           0
      write(OutFileUnit,'(9x,a,11x,a,14x,a,15x,a,7x,a)') 'Ext Coeff','Scat Coeff','SSA','Asym Par','Num of Expan Coeffs'
      write(OutFileUnit,'(3x,4e20.11,9x,i4)') &  ! AMS 09 JAN 2015 Updated expansion coeffs from i3 to i4 as sometimes we have >1000.
         TriAx_bulk(1),& ! Extinction coefficient
         TriAx_bulk(2),& ! Scattering coefficient
         TriAx_bulk(3),& ! Single scattering albedo
         TriAx_asymm  ,& ! Asymmetry parameter
         TriAx_ncoeffs   ! Number of Expansion coefficients

      !Sample:
      !       Moment                                             Expan Coeffs (a1,a2,a3,a4,b1,b2)
      !        0000   0.10000000000E+01   0.00000000000E+00   0.00000000000E+00   0.85192295256E+00   0.00000000000E+00  -0.00000000000E+00
      write(OutFileUnit,'(7x,a,45x,a)') 'Moment','Expan Coeffs (a1,a2,a3,a4,b1,b2)'
      do l=0,TriAx_ncoeffs
         write(OutFileUnit,'(8x,i4.4,6e20.11)') &
            l,(TriAx_expcoeffs(i,l),i=1,6) ! Index, Expansion coefficient output (a1,a2,a3,a4,b1,b2)
      enddo
   else
      write(OutFileUnit,'(8x,a,e20.11)') 'Ref#',Micron_wavelength

      !Sample:
      !         Ext Coeff           Scat Coeff              SSA
      !      0.30186852870E-01   0.29730300500E-01   0.98487578775E+00
      write(OutFileUnit,'(9x,a,11x,a,14x,a)') 'Ext Coeff','Scat Coeff','SSA'
      write(OutFileUnit,'(3x,4e20.11,9x,i3)') &
         TriAx_bulk(1),& ! Extinction coefficient
         TriAx_bulk(2),& ! Scattering coefficient
         TriAx_bulk(3)   ! Single scattering albedo

!  Here is the Mie/Tmat specification
      !Sample:
      !                 ----------------Particle Size Distribution Characteristics----------------
      !        Normalization       Geom Xsection        Geom Volume      Effective radius    Effective variance
      !      0.30186852870E-01   0.29730300500E-01   0.98487578775E+00   0.98487578775E+00   0.98487578775E+00
!      write(OutFileUnit,'(17x,a)') '----------------Particle Size Distribution Characteristics----------------'
!      write(OutFileUnit,'(8x,a,7x,a,8x,a,5x,a,4x,a)') 'Normalization','Geom Xsection','Geom. Volume',&
!                                                      'Effective radius','Effective variance'
!      write(OutFileUnit,'(3x,5e20.11)') MTU_dist(1:5)

!  Diagnostic information is different
      !Sample:
      !                 ----------------Particle Size Distribution Characteristics----------------
      !        RadiusSEsphere      RadiusVEsphere      Projected Area     Particle Volume     (Entry not set)
      !      0.30186852870E-01   0.29730300500E-01   0.98487578775E+00   0.98487578775E+00   0.000000000000E+00
      write(OutFileUnit,'(17x,a)') '----------------Particle Size Distribution Characteristics----------------'
      write(OutFileUnit,'(8x,a,6x,a,6x,a,5x,a,5x,a)') 'RadiusSEsphere','RadiusVEsphere','ProjectedArea',&
                                                      'Particle Volume','(Entry not set) '
      write(OutFileUnit,'(3x,5e20.11)') TriAx_dist(1:5)

   endif

!  Finish

   return
end subroutine Write_TriAxaerosol_file

!  End module

end Module Write_TriAxaerosol_file_m
