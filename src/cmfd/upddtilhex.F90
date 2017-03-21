! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    subroutine upddtilhex

! Calculate D_tilde for hex that is defined by ordinary FDM
 
      use const
      use geomhex
      use tpen
      use xsec,    only : xsd
      use cmfdhex2g, only : dfd, dhat, dfdz, dhatz
      implicit none

      real ::  cofbd(ng2),cofbdz(ng2)
	    integer :: isfc, ind, is1, is2, iz, ig, iz1, iz2, ih
	    real :: dc1, dc2, hzr

      cofbd(1)=rt3*alxr*hside
      cofbd(2)=rt3*alxr*hside

      do isfc=1,nxsfc
        ind=neigsnd(3,isfc)
        if(ind.eq.12) then
          is1=neigsnd(1,isfc)
          is2=neigsnd(2,isfc)
          do iz=1,nz 
            do ig=1,ng2
              dc1=xsd(ig,is1,iz)
              dc2=xsd(ig,is2,iz)
              dfd(ig,isfc,iz)=2.*dc1*dc2/(dc1+dc2)
            enddo
          enddo
        else
          is1=neigsnd(ind,isfc)
          do iz=1,nz 
            do ig=1,ng2
              dc1=xsd(ig,is1,iz)
              dfd(ig,isfc,iz)=dc1*2.*cofbd(ig)/(cofbd(ig)+2.*dc1)
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
          hzr=hz(iz2)/hz(iz1)
          do ih=1,nassy 
            do ig=1,ng2
              dc1=xsd(ig,ih,iz1)
              dc2=xsd(ig,ih,iz2)
              dfdz(ig,ih,iz)=dc2*(1.+hzr)/(hzr+dc2/dc1)
            enddo
          enddo
        else
          iz1=neigsndz(ind,iz)
          if(iz.eq.1) then
             cofbdz(1)=hz(iz1)*alzl
             cofbdz(2)=hz(iz1)*alzl
          else
             cofbdz(1)=hz(iz1)*alzr
             cofbdz(2)=hz(iz1)*alzr
          endif
          do ih=1,nassy 
            do ig=1,ng2
              dc1=xsd(ig,ih,iz1)
              dfdz(ig,ih,iz)=dc1*2.*cofbdz(ig)/(cofbdz(ig)+2.*dc1)
            enddo
          enddo
        endif
      enddo

      return
      end
