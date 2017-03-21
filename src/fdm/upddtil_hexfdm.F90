! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine upddtil_hexfdm

! Calculate D_tilde for hex that is defined by ordinary FDM
 
	   use const
	   use geomhex
	   use hexfdm
	   use xsec,    only : xsdf
	   implicit none

      real ::  cofbd(ng),cofbdz(ng)
	   integer :: isfc, ind, is1, is2, iz, ig, iz1, iz2, ih
	   real :: dc1, dc2, hzr

      real :: beta1, beta2, rt3rh, betabc, beta
      real :: rhz1, rhz2, rhz


      rt3rh=rt3/hdiv
      do isfc=1,nxsfc
         ind=neigsnd(3,isfc)
         if(ind.eq.12) then
            is1=neigsnd(1,isfc)
            is2=neigsnd(2,isfc)
            do iz=1,nz 
               do ig=1,ng
                  beta1=rt3rh*xsdf(ig,is1,iz)
                  beta2=rt3rh*xsdf(ig,is2,iz)
                  dtilrfdm(ig,isfc,iz)=2.*beta1*beta2/(beta1+beta2)
               enddo
            enddo
         else
            is1=neigsnd(ind,isfc)
            do iz=1,nz 
               do ig=1,ng
                  beta1=rt3rh*xsdf(ig,is1,iz)
                  beta2=0.5*alxr
                  dtilrfdm(ig,isfc,iz)=2.*beta1*beta2/(beta1+beta2)
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
               enddo
            enddo
         else
            iz1=neigsndz(ind,iz)
            rhz=1./hz(iz1)
            if(iz.eq.1) then
               betabc=0.5*alzl   
            else
               betabc=0.5*alzr
            endif
            do ih=1,nassy 
               do ig=1,ng
                  beta=xsdf(ig,ih,iz1)*rhz
                  dtilzfdm(ig,ih,iz)=2.*beta*betabc/(beta+betabc)
               enddo
            enddo
         endif
      enddo

      return
      end
