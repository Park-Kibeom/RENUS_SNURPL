! functions for water properties at 15.5 MPa
! cubic polynomial for thermal conductivity, max err=0.0204%
!  15.5 Mpa,  280 < T < 340, k in w/m-C, T in C
    function fcond(t)
      use param
      include 'global.h'
      real :: fcond

! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
      if (iflfr) then
        fcond=3.61+1.517e-02*(t+CKELVIN)-1.741e-06*((t+CKELVIN)**2)
! added end
      else
        fcond=8.9182016e-01+t*(-2.1996892e-03+t*(9.9347652e-06+t*(-2.0862471e-08)))
      endif

      return
    end function
!=======================================================================================!
! cubic polynomial for density as a function of temperature, max err=0.0692%
!  15.5 Mpa,  280 < T < 340, rho in Kg/M^3, T in C
    function fdens(t)
      use param
      include 'global.h'
      real :: fdens

! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
      if (iflfr) then
        fdens=11096.0-1.3236*(t+CKELVIN)
! added end
      else
        fdens=5.9901166e+03+t*(-5.1618182e+01+t*(1.7541848e-01+t*(-2.0613054e-04)))
      endif
      
      return
    end function
!=======================================================================================!
! cubic polynomial for density as a function of enthalpy, max err=0.0112%
!  15.5 Mpa,  280 < T(h) < 340, rho in Kg/M^3, input h in J/Kg
    function fdensh(h)
      use param
      include 'global.h'
      real :: fdensh

! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
      if (iflfr) then
        fdensh=3.4970e-15*h*h*h-9.5718e-09*h*h-7.3728e-03*h+10946.6480
! added end
      else
        hk=0.001*h
        fdensh=1.4532039e+03+hk*(-1.1243975e+00+hk*(7.4502004e-04+hk*(-2.3216531e-07)))
      endif
      
      return
    end function
!=======================================================================================!
! cubic polynomial for enthalpy as a function of temperature, max err=0.0306%
!  15.5 Mpa,  280 < T < 340, output h in J/Kg, T in C
    function fenthal(t)
      use param
      include 'global.h'
      real :: fenthal

! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
      if (iflfr) then
        fenthal=48076.9231+159.1466*(t+CKELVIN-397.7)-2.7225e-2*((t+CKELVIN-397.7)**2)+7.1264e-6*((t+CKELVIN-397.7)**3)
! added end
      else
        y=-5.9301427e+03+t*(6.5488800e+01+t*(-2.1237562e-01+t*(2.4941725e-04)))
        fenthal=y*1000
      endif

      return
    end function
!=======================================================================================!
! quartic polynomial for heat cappacity, max error=0.3053%
!  15.5 Mpa,  280 < T < 340, output h in J/Kg-C, T in C
    function fhcap(t)
      use param
      include 'global.h'
      real :: fhcap

! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
      if (iflfr) then
        fhcap=159.0-2.72e-02*(t+CKELVIN)+7.12e-06*((t+CKELVIN)**2)
! added end
      else
        y=3.0455749e+03+t*(-4.0684599e+01+t*(2.0411250e-01+t*(-4.5526705e-04+t*(3.8115453e-07))))
        fhcap=y*1000
      endif
      
      return
    end function
!=======================================================================================!
! wall-to-coolant heat transfer coeffcient in w/m^2-C
    function fhtcoef(t,deq,rhou)
      use param
      include 'global.h'
      real :: fhtcoef
      real k,mu,cp,pr,re
       
      k=fcond(t)
      mu=fvisco(t)
      cp=fhcap(t)
      pr=cp*mu/k
      re=deq*rhou/mu
      
! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
      if (iflfr) then
        fhtcoef=k/deq*(4.5+0.0156*(re**0.85)*(pr**0.86))
! added end
      else
        fhtcoef=0.023*k/deq*(pr**0.4)*(re**0.8)
      endif     
      
      return
    end function
!=======================================================================================!
! cubic polynomial for temperature as a function of enthalpy, max err=0.0055%
!  15.5 Mpa,  280 < T < 340, T in C, input h in J/Kg
    function ftemp(h)
      use param
      include 'global.h'
      real :: ftemp

! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
      if (iflfr) then
        ftemp=-2.6420e-15*h*h*h+7.5228e-09*h*h+5.5703e-03*h+112.8377-CKELVIN
! added end
      else
        hk=0.001*h
        ftemp=1.4851739e+02+hk*(-1.2764991e-01+hk*(3.0781294e-04+hk*(-9.5429959e-08)))
      endif
      
      return
    end function
!=======================================================================================!
! cubic polynomial for viscosity, max err=0.0641%
!  15.5 Mpa,  280 < T < 340, mu in Pa-sec, T in C
    function fvisco(t)
      use param
      include 'global.h'
      real :: fvisco

! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
      if (iflfr) then
        fvisco=4.94e-04*exp(754.1/(t+CKELVIN))
! added end
      else
        fvisco=9.0836878e-04+t*(-7.4542195e-06+t*(2.3658072e-08+t*(-2.6398601e-11)))
      endif
      
      return
    end function
!=======================================================================================!
! find enthalpy for given temperature, input h is the initial guess
! use fenthal instead of this unless it is really necessary
    function findh(t,h)
      use param
      include 'global.h'
      include 'files.h'

      real :: findh

      x1=1.01*h
      x2=h
      y1=ftemp(x1)
      y2=ftemp(x2)
      do i=1,20
         if(abs(x2-x1).lt.0.01) go to 100
         slope=(y2-y1)/(x2-x1)
         xn=x2+(t-y2)/slope
         x1=x2
         y1=y2
         x2=xn
         y2=ftemp(x2)
      enddo
  100 findh=xn
      if(i.eq.20) then
         write(ioscn,'("Fail to find enthalpy for temperature",e12.5)') t
      endif
      
      return
    end function
!=======================================================================================!