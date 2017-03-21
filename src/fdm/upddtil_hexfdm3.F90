! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine upddtil_hexfdm3

! Calculate D_tilde for hex that is defined by ordinary FDM
 
	   use const
	   use geomhex
	   use hexfdm3
	   use xsec,    only : xsdf
	   implicit none

	   integer :: isfc, ind, is1, is2, iz, ig, iz1, iz2, ih
	   real :: dc1, dc2

      real :: beta1, beta2, rt3rh, betabc(ng), beta
      real :: rhz1, rhz2, rhz


      rt3rh=rt3/hdiv
      do isfc=1,nxsfc
         ind=neigsnd(3,isfc)
         if(ind.eq.12) then
            ! inner node
            is1=neigsnd(1,isfc)
            is2=neigsnd(2,isfc)
            do iz=1,nz 
               do ig=1,ng
                  ! 0-th
                  beta1=rt3rh*xsdf(ig,is1,iz)
                  beta2=rt3rh*xsdf(ig,is2,iz)
                  dtilrfdm(ig,isfc,iz)=2.*beta1*beta2/(beta1+beta2)
                  ! 2-nd
                  beta1=rt3rh*xsdf2(ig,is1,iz)
                  beta2=rt3rh*xsdf2(ig,is2,iz)
                  dtilrfdm2(ig,isfc,iz)=2.*beta1*beta2/(beta1+beta2)
               enddo
            enddo
         else
            ! boundary node
            do iz=1,nz 
               do ig=1,ng
                  dtilrfdm(ig,isfc,iz)=0.
                  dtilrfdm2(ig,isfc,iz)=0.
               enddo
            enddo
         endif
      enddo
      if(nz.eq.1) return   !if 2D

      do iz=1,nz+1 
         ind=neigsndz(3,iz)
         if(ind.eq.12) then
            iz1=neigsndz(1,iz)
            iz2=neigsndz(2,iz)
            rhz1=1./hz(iz1)
            rhz2=1./hz(iz2)
            do ih=1,nassy 
               do ig=1,ng
                  beta1=xsdf(ig,ih,iz1)*rhz1
                  beta2=xsdf(ig,ih,iz2)*rhz2
                  dtilzfdm(ig,ih,iz)=2.*beta1*beta2/(beta1+beta2)

                  beta1=xsdf2(ig,ih,iz1)*rhz1
                  beta2=xsdf2(ig,ih,iz2)*rhz2
                  dtilzfdm2(ig,ih,iz)=2.*beta1*beta2/(beta1+beta2)
               enddo
            enddo
         else
            do ih=1,nassy 
               do ig=1,ng
                  dtilzfdm(ig,ih,iz)=0.
                  dtilzfdm2(ig,ih,iz)=0.
               enddo
            enddo
         endif
      enddo

      return
      end subroutine
