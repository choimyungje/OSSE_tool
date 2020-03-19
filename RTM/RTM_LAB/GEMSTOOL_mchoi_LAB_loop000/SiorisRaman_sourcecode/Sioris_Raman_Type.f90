module Sioris_Raman_Specdata_m

IMPLICIT NONE

INTEGER, PARAMETER        :: dpr = KIND(1.0D0)
INTEGER, PARAMETER        :: N2Jmax=28, O2maxJ=53, O2max=94

TYPE Sioris_Raman_Specdata

  REAL (KIND=dpr), DIMENSION(0:N2Jmax)      :: N2E
  REAL (KIND=dpr), DIMENSION(0:2*N2Jmax-3)  :: N2b
  REAL (KIND=dpr), DIMENSION(0:O2maxJ)      :: O2EnZ
  REAL (KIND=dpr), DIMENSION(0:2*O2max-7)   :: O2E,  O2b
  INTEGER, DIMENSION (0:O2maxJ)            :: O2JZ
  INTEGER, DIMENSION (0:2*O2max-7)         :: O2J2, O2shift
  INTEGER, DIMENSION (0:2*N2Jmax-3)        :: N2shift  

END TYPE Sioris_Raman_Specdata

end module Sioris_Raman_Specdata_m
