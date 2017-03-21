! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    subroutine tpenbc(iftran)

#define p1_nodal_correction   ! 2014_03_12 . scb

! Update Nodal Solution from CMFD results

	    use const
	    use geomhex
	    use tpen
	    use xsec
	    use sfam,     only : phi, psi
	    use tran,     only : betap, rvelotm, srctrf
      use cmfdhex2g, only : dfd, dhat, betaphis, dfdz, dhatz, betaphisz
	    implicit none

	    logical :: iftran
	    real :: alphaz, sumf(ng2),sump(ng2),fnorm(ng2)
	    integer :: m, k, l, iz, ih, ig, ip, it, nn, isfc,   &
			          ind, iz1, iz2, id1, id2
	    real :: vol, dfdm, dhatm, bt, oneflxm, w,    &
			        oneflxn, phis, cjn, hzr, flx1, flx2, &
			        rvol


      do m=1,ng2
         sumf(m)=0
         sump(m)=0
      enddo
      do k=1,nz
         do l=1,nxy
            vol=volnode(l,k)
            do m=1,ng2
               sumf(m)=sumf(m)+hflx(m,l,k)*vol
               sump(m)=sump(m)+phi(m,l,k)*vol
            enddo
         enddo
      enddo
      do m=1,ng2
         fnorm(m)=sump(m)/sumf(m)
      enddo
!      
      if(ng.ne.2) then
         call tpenbcmg(fnorm)
         return
      endif
!
      do iz=1,nz
        do ih=1,nassy
          do ig=1,2
            hflx(ig,ih,iz)=phi(ig,ih,iz)
          enddo
        enddo
        
! 2014_03_11 . scb comment        
#ifndef p1_nodal_correction   ! 2014_03_12 . scb
        do ip=1,nxpnt
          do ig=1,ng2
            pflx(ig,ip,iz)=pflx(ig,ip,iz)*fnorm(ig)
          enddo
        enddo
#endif        
! comment end        
      enddo

! surface flux cal. from phis=w*phin+(1-w)phim+beta*(phim+phin)/2
      do iz=1,nz
        do ih=1,nassy
          do it=1,6
            nn=neignd(it,ih)
            isfc=neigsfc(it,ih)
            do ig=1,2
              dfdm=dfd(ig,isfc,iz)
              dhatm=wtdhat(it,ih)*dhat(ig,isfc,iz) 
              bt=betaphis(ig,isfc,iz)
              oneflxm=hflx(ig,ih,iz)
              if(nn.ne.0) then
                w=dfdm/xsd(ig,ih,iz)/2.
                oneflxn=hflx(ig,nn,iz) 
                phis=w*oneflxn+(1.-w)*oneflxm+bt*(oneflxm+oneflxn)/2.
                cjn=-((dfdm+dhatm)*oneflxn-(dfdm-dhatm)*oneflxm)   &
                      /rt3/hside 
                cnto(ig,it,ih,iz)=(phis+2.*cjn)/4.
              else
                cjn=(dfdm-dhatm)*oneflxm/rt3/hside 
                if(bt.eq.0) then
                   phis=oneflxm/(1+0.25*hf2f/xsd(ig,ih,iz)*alxr)
                else       
                   phis=bt*oneflxm
                endif
                cnto(ig,it,ih,iz)=(phis+2.*cjn)/4.
              endif
            enddo
          enddo
        enddo
      enddo

! update z-direction out-current
      !if(nz.gt.1) then   ! 2012_12_18 . scb
        do iz=1,nz+1 
          ind=neigsndz(3,iz)
          if(ind.eq.12) then
            iz1=neigsndz(1,iz)
            iz2=neigsndz(2,iz)
            id1=neigsndz(4,iz)
            id2=neigsndz(5,iz)
            hzr=hz(iz2)/hz(iz1)
            do ih=1,nassy
              do ig=1,2
                dfdm=dfdz(ig,ih,iz)
                dhatm=dhatz(ig,ih,iz) 
                bt=betaphisz(ig,ih,iz)
                flx1=hflx(ig,ih,iz1)
                flx2=hflx(ig,ih,iz2)
                w=dfdm/xsd(ig,ih,iz1)/(1.+hzr)
                phis=w*flx2+(1.-w)*flx1+bt*(flx1+flx2)/2.
                cjn=-((dfdm+dhatm)*flx2-(dfdm-dhatm)*flx1)*2. &
                      /(hz(iz1)+hz(iz2)) 
                cntzo(ig,id1,ih,iz1)=(phis+2.*cjn)/4.
                cntzo(ig,id2,ih,iz2)=(phis-2.*cjn)/4.
              enddo
            enddo
          else
            iz1=neigsndz(ind,iz)
            id1=neigsndz(ind+3,iz)
            if(iz.eq.1) then
               alphaz=alzl
            else
               alphaz=alzr
            endif
            do ih=1,nassy 
              do ig=1,ng2
                dfdm=dfdz(ig,ih,iz)
                dhatm=dhatz(ig,ih,iz)
                bt=betaphisz(ig,ih,iz)
                if(ind.eq.2) dhatm=-dhatm
                flx1=hflx(ig,ih,iz1)
                cjn=(dfdm-dhatm)*flx1/hz(iz1)
                if(bt.eq.0) then
                   phis=flx1/(1+0.25*hz(iz1)/xsd(ig,ih,iz1)*alphaz)
                else
                   phis=bt*flx1
                endif
                cntzo(ig,id1,ih,iz1)=(phis+2.*cjn)/4.
              enddo
            enddo
          endif
        enddo
      !endif

      if(.not.iftran) return

! calculate effective trasient fixed source
#ifndef FIXED_SOURCE
      !do k=1,nz
      !   do l=1,nxy
      !      rvol=1/volnode(l,k)
      !      sefve(1,l,k)=srctrf(1,l,k)*rvol-rvelotm(1,l,k)*hflx(1,l,k) &
      !                  -(1-betap(l,k))*psi(l,k)*rvol
      !      sefve(2,l,k)=srctrf(2,l,k)*rvol-rvelotm(2,l,k)*hflx(2,l,k)
      !   enddo
      !enddo
#endif

      return
    end
