! file name : thcntl.inc - parameters to control t/h
      integer :: ntth, npint
      common /thcntli/ ntth, npint

      logical fdbk,isflatpower,flagth
      common /thcntl/ fdbk,isflatpower,flagth
      
      real :: epstf    !fuel temp convergence criterion
      common /thcntld/ epstf
      
      real :: thetac,rthetac,thetaf,thetacb,thetafb
      common /thcntld2/ thetac,rthetac,thetaf,thetacb,thetafb

      real :: kgap,kgapb,kgap2,kgap4,kgap4b
      common /dtrcntl/kgap, &                    !effective gap conductivity (see thinit)
                      kgapb,kgap2,kgap4,kgap4b   !see thinit for definition
      