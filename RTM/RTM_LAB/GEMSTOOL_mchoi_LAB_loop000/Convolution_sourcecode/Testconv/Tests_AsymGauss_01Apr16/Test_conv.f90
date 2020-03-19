program testconv
   implicit none

   integer, parameter :: dpm = selected_real_kind(15)
   integer, parameter :: E_ndat      = 5201     ! MW1 maximum size
   integer, parameter :: E_nscenes   = 3        ! Clear/Cloud/Total
   integer, parameter :: E_ndat_c    = 301      ! MW1 maximum size
   integer, parameter :: E_ngeoms    = 1

   integer        :: ngeoms, ndat
   real(kind=dpm) :: lambdas(E_ndat)
   real(kind=dpm) :: intensity_Exact   (E_ndat,E_ngeoms,E_nscenes)

!  Coarse-grid variables. Introduced, 3/29/16

   integer        :: ndat_c
   real(kind=dpm) :: lambdas_c   (E_ndat_c)
   real(kind=dpm) :: solarflux_c (E_ndat_c)
   real(kind=dpm) :: asym_c      (E_ndat_c)
   real(kind=dpm) :: hw1e_c      (E_ndat_c)
   real(kind=dpm) :: IExact_c    (E_ndat_c,E_ngeoms,E_nscenes)

!  help

   integer        :: idum, i, w, iout
   real(kind=dpm) :: dum

!  Get Fine data

   ndat = 2871 ; ngeoms = 1
   open(20,file='ExactScalar_D04_Ray_RS_Obsg_S30V00A010_MW2_CF022_Tr106.Out',status='old')
   do i = 1, ndat
       read(20,*)idum,lambdas(i), dum, intensity_Exact(i,1,1:3), dum,dum,dum,dum,dum,dum,dum,dum,dum
   enddo
   close(20)

!  Get Coarse data

   open(23,file= 'OMI_MW2_Coarse.dat_Fiddled', status = 'old')
   read(23,*)ndat_c
   do w = 1, ndat_c
      read(23,*) lambdas_c(w), solarflux_c(w)
   enddo
   close(23)

!  Compute HWMS = spectral resolution

   hw1e_c(1) = lambdas_c(2) - lambdas_c(1)
   do w = 2, ndat_c - 1
      hw1e_c(w) = 0.5d0* ( (lambdas_c(w+1)-lambdas_c(w)) + (lambdas_c(w)-lambdas_c(w-1)) )
   enddo
   hw1e_c(ndat_c) = lambdas_c(ndat_c) - lambdas_c(ndat_c-1)
   asym_c = 0.0d0

!  Convolve

   IExact_C = 0.0d0
   do iout = 3, 3
     Call asym_gauss_f2c (lambdas(1:ndat), intensity_Exact(1:ndat,1:ngeoms,iout), ndat, ngeoms, ndat_c, &
                          hw1e_c(1:ndat_c), asym_c(1:ndat_c), lambdas_c(1:ndat_c), IExact_C(1:ndat_c,1:ngeoms,iout) )
   enddo

!  Write results for show

   open(24,file= 'ExactScalar_D04_Ray_RS_Obsg_S30V00A010_MW2_CF022_Tr106.COA', status = 'unknown')
   do w = 1, ndat_c
      write(24,'(f10.5,1p3e17.7)')lambdas_c(w), IExact_C(w,1,1:3)
   enddo
   close(24)

!  Done test

   stop
end program testconv


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

!  Precision

integer, parameter :: fpk = SELECTED_REAL_KIND(15)

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


     write(100,*)i,nhalf,midx,sidx,eidx, cwave(i)
     write(*,*)i,nhalf,midx,sidx,eidx, cwave(i)

!mick fix 10/16/2013 - added loop & if block to avoid singularity
     !sign_asym(sidx:eidx) =  (fwave(sidx:eidx) - cwave(i)) / SQRT((fwave(sidx:eidx) - cwave(i))**2)
     DO j=sidx,eidx
       IF ((fwave(j) - cwave(i)) < 1.0E-16_FPK) THEN
         sign_asym(j) = 0.0_FPK
       ELSE
         sign_asym(j) = (fwave(j) - cwave(i)) / SQRT((fwave(j) - cwave(i))**2)
       END IF
       !write(*,*) 'j = ',j,' sign_asym(j) = ',sign_asym(j)
!        write(*,*)j,fwave(j),fwave(j) - cwave(i),sign_asym(j)
     ENDDO

!  Zero point (code by R. Spurr to avoid trouble)

     qqq = MINVAL(MAXLOC(fwave, MASK=(fwave .eq. cwave(i))))
     if ( qqq.gt.0) sign_asym(MINVAL(MAXLOC(fwave, MASK=(fwave .eq. cwave(i))))) = 1.0

     asym_factor2(sidx:eidx) = (1.0 + sign_asym(sidx:eidx) * asym(i)) ** 2

     slit(sidx:eidx) = EXP(-(cwave(i) - fwave(sidx:eidx))**2 / (hw1esq*asym_factor2(sidx:eidx)) )

     DO j=sidx,eidx
        write(101,*)j,fwave(j),slit(j)
     enddo
     stop'first point'

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

SUBROUTINE trapz (x,y,np,aire)

  ! Integration par la methode des trapezes
  !
  ! Dominique Lefebvre janvier 2007

  ! a = borne inferieure d'integration
  ! b = borne superieure d'integration
  ! n = nombre de pas

  ! aire = surface retournee
  IMPLICIT NONE
!  Precision

integer, parameter :: fpk = SELECTED_REAL_KIND(15)

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


