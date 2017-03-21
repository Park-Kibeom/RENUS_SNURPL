! added in ARTOS ver. 0.2 ( for Thermal Expansion ). 2012_07_04 by SCB

! volumetric thermal expansion coefficient
! fvolexp in 1/K, T in K
    function fvcool(t)
      use param
      include 'global.h'
      real :: fvcool

      fvcool = 1./(8383.2-t)
      return
    end function

!=======================================================================================!
! linear thermal expansion coefficient of stainless steel
    function fagrid(t)
      use param
      include 'global.h'
      real :: fagrid
      real :: a0, a1, a2

      a0 = 1.7887e-5
      a1 = 2.3977e-9
      a2 = 3.2692e-13

      fagrid = a0 + a1*t + a2*t*t
		
      return
    end function

!=======================================================================================!
! linear thermal expansion coefficient of HT-9
    function faclad(t)
      use param
      include 'global.h'
      real :: faclad
      real :: a0, a1, a2, a3

      a0 = -0.16256
      a1 = 1.62307e-4
      a2 = 1.42357e-6
      a3 = -5.50344e-10

      faclad = 0.01 * (a0 + a1*t + a2*t*t + a3*t*t*t)

      return
    end function
!=======================================================================================!
! linear thermal expansion coefficient of U-Zr-Pu Alloys
    function fafuel(t)
      use param
      include 'global.h'
      real :: fafuel
      real :: a0, a1, a2, a3

      a0 = -0.424
      a1 = 1.658e-3
      a2 = -1.052e-6
      a3 = 1.115e-9

      fafuel = 0.01 * (a0 + a1*t + a2*t*t + a3*t*t*t)

      return
    end function
