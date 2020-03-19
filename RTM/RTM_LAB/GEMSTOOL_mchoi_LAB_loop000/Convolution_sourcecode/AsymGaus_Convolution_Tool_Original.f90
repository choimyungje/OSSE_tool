Module FASTPCA_Convolution_Tool_m

!  Convolution Tool.

implicit none

!  Precision

integer, parameter :: fpk = SELECTED_REAL_KIND(15)

public  asym_gauss_f2c 
private trapz, fpk

contains

SUBROUTINE asym_gauss_f2c (fwave, fspec, nf, nspec, nc, hw1e, asym, cwave, cspec)

! =========================================================================
!
! Convolves input spectrum with an asymmetric Gaussian slit function of
! specified HW1E (half-width at 1/e intensity)
!
!  g1=exp(-((lam_iter-lam_ch1(i)).^2)./(hwe(i).*(1+sign(lam_iter-lam_ch1(i)).*asym(i))).^2);
!
! The :asymmetric Gaussian g(x) is defined as
!                   _                                      _
!                  |               (x-x0)^2                 |
!      g(x) =  EXP | - ------------------------------------ |
!                  |_     (hw1e*(1+sign(x-x0)*asym))^2     _|
!
!  FWHM = 2.0 * sqrt(ln(2.0)) * hw1e = 1.66551 * hw1e
! =========================================================================

! convolve high-resolution spectra to low resolution spectra
! fwave: high-resolution wavelength grid
! fspec: high-resolution reference spectra (possibly more than 1)
! nf:    number of wavelengths at high-resolution
! nspec: number of spectra
! fwhm:  slit function in terms of FWHM (nm)
! cwave: low-resolution wavelength grid (could be the same as high-resolution grid)
! cspec: spectra at fwhm
! nc:    number of wavelengths for the low-resolution grid

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                         INTENT (IN)         :: nc, nf, nspec
  REAL (KIND=FPK), DIMENSION (nf), INTENT (IN)         :: fwave
  REAL (KIND=FPK), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=FPK), DIMENSION (nc), INTENT (IN)         :: cwave,hw1e,asym

  REAL (KIND=FPK), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                         :: i, j, midx, sidx, eidx, nhalf, qqq
  REAL (KIND=FPK)                 :: hw1esq, dfw, ssum, intf, cspeci
  REAL (KIND=FPK), DIMENSION (nf) :: slit, sign_asym, asym_factor2, fspeci

  dfw  = fwave(2) - fwave(1)
  nhalf  = INT((sum(hw1e) / MAX(1,size(hw1e))) / dfw * 4.0)

  ! Integration done at the fine wavelength grid (fwave) but slit function centred at the coarse wavelength grid (cwave)

  DO i = 1, nc
     ! Find the closest pixel
     midx = MINVAL(MAXLOC(fwave, MASK=(fwave <= cwave(i)))) + 1

     sidx = MAX(midx - nhalf, 1)
     eidx = MIN(nf, midx + nhalf)

     hw1esq = hw1e(i) ** 2

!mick fix 10/16/2013 - added loop & if block to avoid singularity
     !sign_asym(sidx:eidx) =  (fwave(sidx:eidx) - cwave(i)) / SQRT((fwave(sidx:eidx) - cwave(i))**2)
     DO j=sidx,eidx
       IF ((fwave(j) - cwave(i)) < 1.0E-16_FPK) THEN
         sign_asym(j) = 0.0_FPK
       ELSE
         sign_asym(j) = (fwave(j) - cwave(i)) / SQRT((fwave(j) - cwave(i))**2)
       END IF
       !write(*,*) 'j = ',j,' sign_asym(j) = ',sign_asym(j)
     ENDDO

!  Zero point (code by R. Spurr to avoid trouble)

     qqq = MINVAL(MAXLOC(fwave, MASK=(fwave .eq. cwave(i))))
     if ( qqq.gt.0) sign_asym(MINVAL(MAXLOC(fwave, MASK=(fwave .eq. cwave(i))))) = 1.0

     asym_factor2(sidx:eidx) = (1.0 + sign_asym(sidx:eidx) * asym(i)) ** 2

     slit(sidx:eidx) = EXP(-(cwave(i) - fwave(sidx:eidx))**2 / (hw1esq*asym_factor2(sidx:eidx)) )

     ! Trapezoidal integration
     call trapz(fwave(sidx:eidx),slit(sidx:eidx),eidx-sidx+1,ssum)

     DO j = 1, nspec
         fspeci(sidx:eidx) = fspec(sidx:eidx, j)

         ! Trapezoidal integration
         call trapz(fwave(sidx:eidx), fspeci(sidx:eidx)*slit(sidx:eidx), eidx-sidx+1, intf)

         ! Normalization by slit function area
         if (ssum.gt.0.0) then
             cspeci = intf / ssum
             cspec(i, j) = cspeci
         endif
     ENDDO
  ENDDO

END SUBROUTINE asym_gauss_f2c

! Same as above, except for correcting solar I0 effect using high resolution
! solar reference spectrum following the Beer's Law
! i0: high resolution solar reference
! scalex: scaling, normally number of molecules
SUBROUTINE asym_gauss_f2ci0(fwave, fspec, i0, nf, nspec, nc, scalex, hw1e, asym, cwave, cspec)

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                                INTENT (IN)  :: nc, nf, nspec
  REAL (KIND=FPK),                        INTENT (IN)  :: scalex
  REAL (KIND=FPK), DIMENSION (nf), INTENT (IN)         :: fwave, i0
  REAL (KIND=FPK), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=FPK), DIMENSION (nc), INTENT (IN)         :: cwave, hw1e, asym
  REAL (KIND=FPK), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                             :: i, ntemp
  REAL (KIND=FPK), DIMENSION(nc)        :: ci0
  REAL (KIND=FPK), DIMENSION(nc, nspec) :: cabspec
  REAL (KIND=FPK), DIMENSION(nf, nspec) :: abspec

  ntemp = 1
  CALL asym_gauss_f2c (fwave, i0, nf, ntemp, nc, hw1e, asym, cwave, ci0)

  ! Follow Beer's Law
  DO i = 1, nspec
     abspec(:, i) = i0 * EXP(-fspec(:, i) * scalex)
  ENDDO

  CALL asym_gauss_f2c (fwave, abspec, nf, nspec, nc, hw1e, asym, cwave, cabspec)

  DO i = 1, nspec
     cspec(:, i) = - LOG(cabspec(:, i)/ ci0) / scalex
  ENDDO

END SUBROUTINE asym_gauss_f2ci0


SUBROUTINE trapz (x,y,np,aire)

  ! Integration par la methode des trapezes
  !
  ! Dominique Lefebvre janvier 2007

  ! a = borne inferieure d'integration
  ! b = borne superieure d'integration
  ! n = nombre de pas

  ! aire = surface retournee
  IMPLICIT NONE

  INTEGER,                         INTENT (IN)  :: np
  REAL (KIND=FPK), DIMENSION (np), INTENT (IN)  :: x,y
  REAL (KIND=FPK),                 INTENT (OUT) :: aire

  real(KIND=FPK) a,b,h

  INTEGER i

  ! Boucle d'integration
  ! Initialisation des variables

  aire = 0.0

  ! Boucle de calcul
  ! h = (b-a)/n. Je sais que l'aire de chaque petit trapèze est Ai = (h/2)*(f(a+ih) + f(a+(i-1)h)).

  DO i = 1, np-1
    b = x(i+1)
    a = x(i)
    h = (b-a) 
    if (h.gt.0.0) aire = aire + (h / 2.0D0) * (y(i+1) + y(i))
  ENDDO
END SUBROUTINE trapz


!  End module

end module FASTPCA_Convolution_Tool_m


