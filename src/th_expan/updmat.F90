! added in ARTOS ver. 0.2 ( for Thermal Expansion ). 2012_07_04 by SCB	
    subroutine updmat(comp, typea, rodfrac)

      use param
      use BasicParam
      use MaterialMod
      use trinx_cntl

      implicit none

      include 'thexpan.inc'

      type(Composition) :: comp

      real, intent(in) :: rodfrac
      integer, intent(in) :: typea  ! 1 = FUEL ASSEMBLY
                                      ! 2 = CONTROL ASSEMBLY
                                      ! 3 = REFLECTOR ASSEMBLY

      integer :: m, i, niso, nmat, ifa
      real :: fuelv, cladv, rod_vol
      integer,parameter :: FUEL=1, CLAD=2, COOL=3

      nmat = comp%no_material

      select case(typea)
        case(COREM) 
          if(.not.thexpan) return
          goto 100 
        case(CONTM)
          goto 200
        case(RREFL)
          if(.not.thexpan) return
          goto 300
        case(AREFL)
          if(.not.thexpan) return
          goto 400
        case(COREL) 
          if(.not.thexpan) return
          goto 500 
        case(COREU) 
          if(.not.thexpan) return
          goto 600 
      end select

100   continue ! fuel assembly

      fuelv = 0.0
      cladv = 0.0 
      do m = 1, nmat
        niso = comp%material(m)%no_isotope
        do ifa=1, nfuelassm	
          if(trim(comp%mat_name(m)).eq.trim(fuelmat(ifa,FUEL))) then
            do i = 1, niso
              comp%material(m)%iso_density(i) = comp%material(m)%iso_density0(i)*adfuel
            enddo
            comp%vol_fraction(m) = comp%vol_fraction0(m)*vfuel
            fuelv = comp%vol_fraction(m)
            call make_temp_coef(comp%material(m))
          elseif(trim(comp%mat_name(m)).eq.trim(fuelmat(ifa,CLAD))) then
            do i = 1, niso
              comp%material(m)%iso_density(i) = comp%material(m)%iso_density0(i)*adclad
            enddo
            comp%vol_fraction(m) = comp%vol_fraction0(m)*vclad
            cladv = comp%vol_fraction(m)
            call make_temp_coef(comp%material(m))
          endif
        enddo
      enddo

      do m = 1, nmat
        niso = comp%material(m)%no_isotope
        do ifa=1, nfuelassm				
          if(trim(comp%mat_name(m)).eq.trim(fuelmat(ifa,COOL))) then 
            comp%vol_fraction(m) = 1-fuelv-cladv
#ifdef EXPANSION
            adcool=1./(comp%vol_fraction(m)*vnode)
#endif
            do i = 1, niso
              comp%material(m)%iso_density(i) = comp%material(m)%iso_density0(i)*adcool
            enddo
            call make_temp_coef(comp%material(m))
          endif
        enddo
      enddo
      
      return

200   continue ! control assembly
      do m=1, nmat
        niso = comp%material(m)%no_isotope	
        ! check the coolant in CA
        if(comp%mat_name(m).eq.contmat(2)) then
          comp%vol_fraction(m)=comp%vol_fraction0(m)*(1.-rodfrac)
          rod_vol=comp%vol_fraction0(m)*rodfrac
          call make_temp_coef(comp%material(m))
          exit
        endif
      enddo

      do m=1, nmat
        ! check the control rod in CA
        if(comp%mat_name(m).eq.contmat(1)) then
          comp%vol_fraction(m)=rod_vol
          call make_temp_coef(comp%material(m))
          exit
        endif
      enddo	
      
      return

300   continue ! radial reflector 
      do m = 1, nmat
#ifdef EXPANSION
        adback = 1./(ar*ar*azf)
#endif
        niso = comp%material(m)%no_isotope
        if(trim(comp%material(m)%name).eq.trim(coolmat(1))) then
          do i = 1, niso
            comp%material(m)%iso_density(i) = comp%material(m)%iso_density0(i)*adback
          enddo
          comp%vol_fraction(m) = comp%vol_fraction0(m)
          call make_temp_coef(comp%material(m))
        endif
      enddo		
      
      return

400   continue ! axial reflector 
      do m = 1, nmat
#ifdef EXPANSION
        adcool = 1./(ar*ar)
#endif
        niso = comp%material(m)%no_isotope
        if(trim(comp%material(m)%name).eq.trim(coolmat(1))) then
          do i = 1, niso
            comp%material(m)%iso_density(i) = comp%material(m)%iso_density0(i)*adcool
          enddo
          comp%vol_fraction(m) = comp%vol_fraction0(m)
          call make_temp_coef(comp%material(m))
        endif
      enddo		
      return

500   continue ! lower plenum
600   continue ! upper plenum
      m=1
      niso = comp%material(m)%no_isotope	
      comp%vol_fraction(m) = comp%vol_fraction0(m)*vpin	
      do i = 1, niso
        comp%material(m)%iso_density(i) = comp%material(m)%iso_density0(i)*adpin
      enddo
      call make_temp_coef(comp%material(m))

      m=2
#ifdef EXPANSION
      adcool = 1./(ar*ar*azf)
#endif
      niso = comp%material(m)%no_isotope	
      comp%vol_fraction(m) = 1.-comp%vol_fraction(1)
      do i = 1, niso
        comp%material(m)%iso_density(i) = comp%material(m)%iso_density0(i)*adcool
      enddo
      call make_temp_coef(comp%material(m))
        
      return
    end subroutine
