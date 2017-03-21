! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      subroutine cheby(iffront,iout)

      use fdm_solver
      implicit none
!
! extrapolate fission source using the chebyshev acceleration parameters
!
      integer :: iout
      logical :: iffront
!
      integer :: mcpm1, m
      integer :: ioutcb0
      integer :: ninfix
      real :: gamma1, theta, acosh, cosh
      character*132 mesg
 
      if(.not.iffront) go to 1000
!
! frontend cheby
!
      if(chebyon) then
! set new value of dominance ratio to be used in the extrapolation
         if(newcheby) then
            sigma=sigbar
            if(iout.le.6) sigma=min(sigma,0.900)
            if(iout.le.9) sigma=min(sigma,0.950)
            if(iout.le.12) sigma=min(sigma,0.985)
            alphacb=2/(2-sigma)
            betacb=0.0
            newcheby=FALSE
            mcp=1
            erfcb=1
         else
            mcpm1=mcp
            mcp=mcp+1
            gamma1=2/sigma-1
            theta=acosh(gamma1)
            alphacb=4*cosh(mcpm1*theta)/(sigma*cosh(mcp*theta))
            betacb=(1-sigma*0.5)*alphacb-1
         endif
      else
         mcp=0
      endif
!
      return
!
! backend cheby
!
 1000 continue
!
! turn on chebyshev
      if(.not.chebyon .and. ncheby.gt.0) then
         if(domr.gt.domrlo .and. domr.lt.domrhi .and. iout-ioutcb0.ge.icheby0) then
            chebyon=TRUE
            newcheby=TRUE
         endif
      endif
      sigbar=domr       
!
! monitor the effectiveness of the current chebyshev cycle 
! and estimate new dominance ratio
      if(chebyon .and. mcp.gt.1) then
         erfcb=erfcb*domr
         gamma1=2/sigma-1
         theta=acosh(gamma1)
         erfcb0=1/cosh(mcpm1*theta)  !theoretical error reduction factor
         erfratio=erfcb/erfcb0       !relative 
         mcpm1=mcp-1
         if(erfratio.lt.1.0) then
            sigbar=sigma*(1+cos(acos(erfratio)/mcpm1))*0.5
         else
            sigbar=sigma*(1+cosh(acosh(erfratio)/mcpm1))*0.5
         endif
!
! reset polynomial if no further improvement is expected
         if(mcp.ge.3) then
!
            if(mcp.ge.ncheby .or. erfratio.ge.1.05) newcheby=TRUE
! turn off chebyshev if the new estimate of dominance ratio becomes too high           
            if(sigbar.ge.domrhi) then
               newcheby=TRUE              
               chebyon=FALSE
               do m=1,ng
!                 ninfix(m)=nint(1.2*ninfix(m))
!                 write(mesg,'(a,i5,a,i5)') 'Increasing Ninners to',ninfix(m),' for Group',m
!                 call message(FALSE,FALSE,mesg)
               enddo
               ioutcb0=iout
               firstdom=TRUE
            endif
         endif
      endif
!
      return
      end subroutine
!
!=================================================
!
      function acosh(x)
! compute inverse cosine hyperbolic of x
      real :: acosh, x
      acosh=log(x+sqrt(x*x-1.d0))
      end function
