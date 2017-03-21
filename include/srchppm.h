! file name- srchppm.h
! variables to be used while searching the amount of boron or power value to achieve criticality
      logical :: srchppm
      logical :: srchmwt
      real :: ppm,targetk

      common /srchppml/ srchppm, srchmwt
      common /srchppmf/ ppm,targetk