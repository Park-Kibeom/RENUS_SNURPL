! functions for thermal properties of fuel regions
! thermal conductivity of uo2 in w/m-C, t in K
    function fkf(t)
      use param
      
      include 'global.h'
      include 'thfuel.inc'

      real ::  fkf

! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
      if (iflfr) then
        fkf=17.5*((1-2.23*wzr)/(1+1.61*wzr)-2.62*wpu)+1.54e-2*((1+0.061*wzr)/(1+1.61*wzr)+0.9*wpu)*t+9.38e-6*(1-2.7*wpu)*t*t
! added end
      else
        fkf=akfuel(0)+t*(akfuel(1)+t*(akfuel(2)+t*akfuel(3)))+akfuel(4)/(t-akfuel(5))  
      endif
          
      return      
    end function
!=======================================================================================!
! thermal conductivity of Zr in w/m-C, t in K
    function fkc(t)
      use param
      
      include 'global.h'
      include 'thfuel.inc'
!
      real ::  fkc

! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
      if (iflfr) then
        fkc=17.622+2.428e-2*t-1.696e-5*t*t
! added end
      else
        fkc=akclad(0)+t*(akclad(1)+t*(akclad(2)+t*akclad(3)))
      endif
      
      return
    end function
!=======================================================================================! 
! volumetric heat capacity of uo2 in J/m^3-C, t in K
    function frhocpf(t)
      use param
      
      include 'global.h'
      include 'thfuel.inc'

      real ::  frhocpf
      
      frhocpf=arcpfuel(0)+t*(arcpfuel(1)+t*(arcpfuel(2)+t*arcpfuel(3)))
      return
    end function
!=======================================================================================! 
! volumetric heat capacity of Zr in J/m^3-C, t in K
    function frhocpc(t)
      use param
      
      include 'global.h'
      include 'thfuel.inc'

      real ::  frhocpc

      frhocpc=arcpclad(0)+t*(arcpclad(1)+t*(arcpclad(2)+t*arcpclad(3)))
      return
    end function
!=======================================================================================! 
! fuel volume avg. temperature
    function ftfavg(x,nr,dr2)
      use param
      
      include 'global.h'
!
      real ::  ftfavg

      dimension x(*) !dmm
      xavg=0
      i=1
      a=(x(i+1)-x(i))/dr2
      xavg=(7*x(i)+x(i+1))*0.25
      do l=2,nr
         area=2/3.*l*x(l+1)+44/3.*i*x(l)+2/3.*(i-1)*x(l-1)
         xavg=xavg+area
         i=l
      enddo
      nrm1=nr-1
      area=(16/3.*nrm1+4+5/24.)*x(nr+1)+(10/3.*nrm1+2.25)*x(nr)-(2/3.*nrm1+11/24.)*x(nr-1)
      ftfavg=(xavg+area)*dr2/(8.*(nr*nr*dr2))
      return
    end function
!=======================================================================================! 
    function fenthalf(tcel)
      use param
      
      include 'global.h'

      real ::  fenthalf

      t=tcel+273.15
      fenthalf=(162.3+t*(0.1519+t*(-7.970e-5+t*1.601e-8)))*t
!
      return
    end function
! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
!=======================================================================================!
! average linear thermal expansion of HT9
    function falphabc(t)
      use param
      
      include 'global.h'
      real ::  falphabc, tref

      tref=282.2649

      falphabc=(-2.191e-3+5.678e-6*t+8.111e-9*t*t-2.576e-12*t*t*t)/(t-tref)

      return
    end function
!=======================================================================================!
    function falphac(t0,t1)
      use param
      
      include 'global.h'
      real ::  falphac

      tref=282.2649

      falphac=falphabc(t1)+(t0-tref)/(t1-t0)*(falphabc(t1)-falphabc(t0))

      return
    end function
! added end
