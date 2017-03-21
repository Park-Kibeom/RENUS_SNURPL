! added in ARTOS ver. 0.2 . 2012_07_24 by SCB     
   subroutine inittpen_sp3

! variable initialization

	   use const
	   use geomhex
      use tpen,  only : fhflx, fohflx, fcnto, focnto, fcntzo, focntzo
	   use tpen_sp3
	   use xsec
	   use sfam,  only : phi_p2

	   implicit none

	   real :: phi0(2,ng2),phi0f(2,ng2),cnt0(2,ng2),phifc(ng2)
	   real :: rdat, fflx
	   integer :: m, mf, ig, iz, ih, it, ip, nnz, im

! coarse group index for each fine group
      do m=1,ng2
         do mf=mgb(m),mge(m)
            igc(mf)=m
         enddo
      enddo 
      r9hs2=2./9.d0/(hside*hside)

      rdat=1./74165.d0
      chlval(1)=-22230.*rdat
      chlval(2)=2590.*rdat
      chlval(3)=-2187.*rdat
      chlval(4)=58968.*rdat
      chlval(5)=-600.*rdat
      chlval(6)=4704.*rdat
      chlval(7)=-10512.*rdat

      if( ng .le. 12 ) then
          nitrpfc=3
      else
          nitrpfc=7
      endif

	   phi_p2(1,:,:,:)=1.0
	   phi_p2(2,:,:,:)=0.1
 
      do ig=1,ng2
         do im=1,2
            phi0(im,ig)=phi_p2(im,ig,1,1)
            phi0f(im,ig)=phi0(im,ig)/(mge(ig)-mgb(ig)+1)
            cnt0(im,ig)=0.25*phi0f(im,ig)
         enddo
      enddo

      r9hs2=2./9.d0/(hside*hside)
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


      do iz=1,nz        
         do ih=1,nassy
            do ig=1,ng2
               do im=1,2
                  hflx2(im,ig,ih,iz)=phi0(im,ig)
               enddo
            enddo
            do it=1,ntph
               do ig=1,ng
                  do im=1,2
                     cnto2(im,ig,it,ih,iz)=cnt0(im,igc(ig))
                  enddo
               enddo
            enddo
            do ig=1,ng
               do im=1,2
                  cntzo2(im,ig,1,ih,iz)=cnt0(im,igc(ig)) 
                  cntzo2(im,ig,2,ih,iz)=cnt0(im,igc(ig))
               enddo 
            enddo
         enddo 
         do ip=1,ncorn
            do ig=1,ng 
               do im=1,2
                  pflx2(im,ig,ip,iz)=phi0f(im,igc(ig)) 
               enddo   
            enddo
         enddo
      enddo 

!      if(ng.ne.2) then
         do iz=1,nz 
            do ih=1,nassy 
               do ig=1,ng
                  do im=1,2
                     fflx=phi0f(im,igc(ig))
                     hflxf2(im,ig,ih,iz)=fflx  
                     fhflx2(im,ig,ih,iz)=fflx/phi0(im,igc(ig))
                     fohflx2(im,ig,ih,iz)=fhflx2(im,ig,ih,iz)
                     do it=1,6
                        aflx2(im,ig,it,ih,iz)=fflx
                        fcnto2(im,ig,it,ih,iz)=fhflx2(im,ig,ih,iz)
                        focnto2(im,ig,it,ih,iz)=fhflx2(im,ig,ih,iz)
                     enddo
                     do it=1,2
                        fcntzo2(im,ig,it,ih,iz)=fhflx2(im,ig,ih,iz)
                        focntzo2(im,ig,it,ih,iz)=fhflx2(im,ig,ih,iz)
                     enddo
                  enddo ! im
               enddo ! ig
            enddo ! ih
         enddo ! iz
!      endif ! ng

      do iz=1,nz 
         do ih=1,nassy 
            do ig=1,ng
               fflx=phi0f(1,igc(ig))
               hflxf2(1,ig,ih,iz)=phi0f(1,igc(ig))
               hflxf2(2,ig,ih,iz)=phi0f(2,igc(ig))
               fhflx(ig,ih,iz)=fflx/phi0(1,igc(ig))
               fohflx(ig,ih,iz)=fhflx(ig,ih,iz)
               do it=1,6
                  aflx2(:,ig,it,ih,iz)=fflx
                  fcnto(ig,it,ih,iz)=fhflx(ig,ih,iz)
                  focnto(ig,it,ih,iz)=fhflx(ig,ih,iz)
               enddo
               do it=1,2
                  fcntzo(ig,it,ih,iz)=fhflx(ig,ih,iz)
                  focntzo(ig,it,ih,iz)=fhflx(ig,ih,iz)
               enddo
            enddo ! ig
         enddo ! ih
      enddo ! iz


      return
   end subroutine

