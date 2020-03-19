module tmat_parameters

      implicit none

!  Type kinds

      integer, parameter :: tmat_dkind  = SELECTED_REAL_KIND(15)
      integer, parameter :: tmat_qkind  = SELECTED_REAL_KIND(29)
      integer, parameter :: tmat_rkind  = SELECTED_REAL_KIND(7)
      integer, parameter :: tmat_fpkind = tmat_dkind
      integer, parameter :: tmat_cpkind = tmat_dkind

!  Dimensioning file for the T-matrix code

      integer npn1, npng1, npn4
      integer npng2, npl, npl1, npn2, npn3, npn5, npn6

!  fundamental numbers

      parameter ( npn1  = 150 )
      parameter ( npng1 = 400 )
      parameter ( npn4  = 100 )

!  original settings

!      parameter ( npn1  = 100 )
!      parameter ( npng1 = 300 )
!      parameter ( npn4  = 80  )

!  R2 limit settings

!      parameter ( npn1  = 100 )
!      parameter ( npng1 = 350 )
!      parameter ( npn4  = 100  )

!  derived dimensioning

      parameter ( npn2  = 2 * npn1  )
      parameter ( npn5  = 2 * npn4  )
      parameter ( npng2 = 2 * npng1 )

      parameter ( npn3  = npn1 + 1  )
      parameter ( npn6  = npn4 + 1  )

      parameter ( npl   = npn2 + 1  )
      parameter ( npl1  = npn5 + 1  )

!  New parameters

      integer   maxnpa
      parameter ( maxnpa = 91 )

!  original code

!      PARAMETER (NPN1=100, NPNG1=300, NPNG2=2*NPNG1, NPN2=2*NPN1,  
!     &           NPL=NPN2+1, NPN3=NPN1+1,  
!     &           NPN4=80, NPN5=2*NPN4, NPN6=NPN4+1, NPL1=NPN5+1)

!  End of file.

! Numbers such as 1, 2, pie, etc.......

  REAL (KIND=8), PARAMETER :: d_zero  = 0.0d0
  REAL (KIND=8), PARAMETER :: d_one   = 1.0d0
  REAL (KIND=8), PARAMETER :: d_two   = 2.0d0
  REAL (KIND=8), PARAMETER :: d_three = 3.0d0
  REAL (KIND=8), PARAMETER :: d_four  = 4.0d0
  REAL (KIND=8), PARAMETER :: d_half  = 0.5d0

!  everything public here

public

end module tmat_parameters
