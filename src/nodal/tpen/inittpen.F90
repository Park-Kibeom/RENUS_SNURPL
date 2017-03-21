! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine inittpen

! variable initialization

	   use const
	   use geomhex
	   use tpen
	   use xsec
	   use sfam,  only : phif, phi

	   implicit none

	   real :: phi0(ng2),phi0f(ng2),cnt0(ng2),phifc(ng2)
	   real :: rdat, fflx
	   integer :: m, mf, ig, iz, ih, it, ip, nnz

!
! coarse group index for each fine group
      do m=1,ng2
         do mf=mgb(m),mge(m)
            igc(mf)=m
         enddo
      enddo 
!	phi(1,:,:)=40.
!	phi(2,:,:)=10.  
!
      if( ng .le. 12 ) then
          nitrpfc=3
      else
          nitrpfc=7
      endif
!
      do ig=1,ng2
         phi0(ig)=phi(ig,1,1)
         phi0f(ig)=phi0(ig)/(mge(ig)-mgb(ig)+1)
         cnt0(ig)=0.25*phi0f(ig)
      enddo
!
      do iz=1,nz        
         do ih=1,nassy
            do ig=1,ng2
               hflx(ig,ih,iz)=phi0(ig)
            enddo
            do it=1,ntph
               do ig=1,ng 
                  cnto(ig,it,ih,iz)=cnt0(igc(ig))
               enddo
            enddo
            do ig=1,ng
               cntzo(ig,1,ih,iz)=cnt0(igc(ig)) 
               cntzo(ig,2,ih,iz)=cnt0(igc(ig)) 
            enddo
         enddo 
         do ip=1,ncorn
            do ig=1,ng 
               pflx(ig,ip,iz)=phi0f(igc(ig))    
            enddo
         enddo
      enddo 

! set values commonly needed for calculating nodal coupling coeff.
      r9hs2=2.d0/9.d0/(hside*hside)
      if(nz.eq.1) go to 100
      do iz=1,nz 
        do it=1,2
          nnz=neigz(it,iz) 
          if(nnz.eq.0) then
            rhzbar2(it,iz)=1./hz(iz)/hz(iz)
          else
            rhzbar2(it,iz)=2./hz(iz)/(hz(iz)+hz(nnz))
          endif
        enddo
      enddo 

100   continue

      rdat=1./74165.d0
      chlval(1)=-22230.*rdat
      chlval(2)=2590.*rdat
      chlval(3)=-2187.*rdat
      chlval(4)=58968.*rdat
      chlval(5)=-600.*rdat
      chlval(6)=4704.*rdat
      chlval(7)=-10512.*rdat

      if(ng.ne.2) then
        do iz=1,nz 
          do ih=1,nassy 
            do ig=1,ng
              fflx=phi0f(igc(ig))
              hflxf(ig,ih,iz)=fflx  
!              phiadjf(ig,ih,iz)=fflx                 !temp 
              fhflx(ig,ih,iz)=fflx/phi0(igc(ig))
              fohflx(ig,ih,iz)=fhflx(ig,ih,iz)
              do it=1,6
                aflx(ig,it,ih,iz)=fflx
                fcnto(ig,it,ih,iz)=fhflx(ig,ih,iz)
                focnto(ig,it,ih,iz)=fhflx(ig,ih,iz)
              enddo
              do it=1,2
                fcntzo(ig,it,ih,iz)=fhflx(ig,ih,iz)
                focntzo(ig,it,ih,iz)=fhflx(ig,ih,iz)
              enddo
            enddo
          enddo
        enddo
      endif

      return
      end subroutine

