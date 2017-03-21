! added in ARTOS ver. 0.2 ( for Thermal Expansion ). 2012_07_04 by SCB	
    subroutine updexpan(atype)

      use param
      use trinx_cntl

      include 'global.h'      
      include 'thcntl.inc'
      include 'thop.inc'
      include 'thexpan.inc'
      include 'thlink.inc'

      integer, intent(in) :: atype ! 1 = FUEL ASSEMBLY
      ! 2 = CONTROL ASSEMBLY
      ! 3 = REFLECTOR		
      real :: lfuel, lclad, lpin
      real, parameter :: TBASE = 293.0    ! added in ARTOS ver. 0.2 ( for Thermal Expansion ). 2012_07_24 by SCB

      if(ifreact) then
        tfuelavg = base_ft + del_ax
        tcoolavg = base_ct + del_rad
        tfueln   = base_ft + del_ax
        tcooln   = base_ct + del_rad

        ! radial reflector coolant
        tc       = base_ct + del_den
        bcoolr   = fvcool(tc)

        ! core coolant
        tc       = base_ct + del_den
        bcoola   = fvcool(tc)
      endif

      bcool0   = fvcool(623.15)		

! update geometries, atom densties and volume fractions.
!
! thermal expansion	
      if(ilink.lt.0) then
        tinavg = tcoolavg
      endif	
      
      select case(atype)	
        case(0) ! update the core geometry
          ! core
          dtfuel = tfuelavg-TBASE
          dtcool = tinavg-TBASE

          !ar = 1.+agrid*dtcool ! radial expansion
          !azf = 1.+aclad*dtfuel ! axial expansion

          ar = 1.+fagrid(tinavg)*dtcool
          azf = 1.+faclad(tfuelavg)

          azc = 1.+faclad(base_ft)

          ! radial reflector
          if(.not.ifreact) then
            bcoolr  = fvcool(tinavg)
          endif
          adback = bcool0/bcoolr
        case(COREM) ! update the fuel assembly geometry
          ! fuel 	
          dtfuel = tfueln-TBASE
          dtcool = tcooln-TBASE

          !lfuel = 1.+afuel*dtfuel ! fuel expansion
          !lclad = 1.+aclad*dtfuel ! cladding expansion

          lfuel = 1.+fafuel(tfueln) ! fuel expansion
          lclad = 1.+faclad(tfueln) ! cladding expansion

          vnode = ar*ar*azf
          vfuel = lfuel*lfuel*azf
          vfuel = vfuel/vnode
          vclad = lclad*lclad*azf
          vclad = vclad/vnode

          adfuel = 1./(vfuel*vnode)
          adclad = 1./(vclad*vnode)
          if(.not.ifreact) then
            bcoola  = fvcool(tcooln)
          endif
          adcool = bcool0/bcoola
        case(COREL) ! update uppde plenum geometry
          dtcool = tcooln-TBASE

          !lpin = 1.+aclad*dtcool ! fuel expansion

          lpin = 1.+faclad(tcooln)

          vnode = ar*ar
          vpin = lpin*lpin
          vpin = vpin/vnode

          adpin = 1./(vpin*vnode)
          if(.not.ifreact) then
            bcoola  = fvcool(tcooln)
          endif
          adcool = bcool0/bcoola				
        case(COREU) ! update uppde plenum geometry
          dtcool = tcooln-TBASE

          !lpin = 1.+aclad*dtcool ! fuel expansion

          lpin = 1.+faclad(tcooln)

          vnode = ar*ar
          vpin = lpin*lpin
          vpin = vpin/vnode

          adpin = 1./(vpin*vnode)
          if(.not.ifreact) then
            bcoola  = fvcool(tcooln)
          endif
          adcool = bcool0/bcoola								
        case(CONTM)
          return											
        case(AREFL)
          ! axial refelctor
          if(.not.ifreact) then
            bcoola  = fvcool(tcooln)
          endif
          adcool = bcool0/bcoola	
          return
        case default
          return
      end select

      return
    end subroutine
