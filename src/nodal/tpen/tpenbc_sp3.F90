! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    subroutine tpenbc_sp3(iftran)

#define sp3_nodal_correction   ! 2014_03_12 . scb    
! Update Nodal Solution from CMFD results

	    use const
	    use geomhex
	    use tpen_sp3
	    use xsec,     only : xsd2, xsd
      use cmfdhex2g !, only : dfd, dhat, betaphis, dfdz, dhatz, betaphisz, dfd2, betaphis2

      use sfam,     only : phi, phif   ! 2014_02_27 . phif added
      use sfam,     only : psi                                ! 2012_12_07 . scb
      use bdf                                                 ! 2012_12_06 . scb
      use tran,     only : betap, rvelotm, srctrf, srctrf2    ! 2012_12_07 . scb 
      
	    implicit none

	    logical :: iftran
	    integer :: m, k, l, iz, ih, ig, ip, it, nn, isfc,   &
			        ind, iz1, iz2, id1, id2, im

! 2014_02_27 . scb . modified array size to ng      
	    real :: alphaz, sumf(ng), sump(ng), fnorm(ng), &
              sumf2(ng), sump2(ng), fnorm2(ng)
	    real :: vol, dfdm, dhatm, bt, oneflxm, w,    &
			        oneflxn, hzr, flx1, flx2,            &
			        rvol
      real :: dc, phis, cjn, phis2
      
	    logical :: transient ! 2014_03_02 . scb      
      common / sp3dbg / transient ! 2014_03_02 . scb
      
! 2014_02_27 . scb : sumf routines are commented because facscb covers them.

      transient=iftran  ! 2014_03_02 . scb
           
! 2014_03_12 . scb      
#ifndef sp3_nodal_correction   ! 2014_03_12 . scb      
      do m=1,ng2
         sumf(m)=0.d0
         sumf2(m)=0.d0   ! 2012_12_13 . scb
         sump(m)=0.d0
      enddo
      
! 2014_02_27 . scb . hflx -> hflxf , ng2 -> ng
      do k=1,nz
         do l=1,nxy
            vol=volnode(l,k)
            do m=1,ng
               sumf(m)=sumf(m)+hflxf2(1,m,l,k)*vol
               sumf2(m)=sumf2(m)+hflxf2(2,m,l,k)*vol  ! 2012_12_13 . scb
               sump(m)=sump(m)+phif(m,l,k)*vol
            enddo
         enddo
      enddo
      do m=1,ng
         fnorm(m)=sump(m)/sumf(m)
         fnorm2(m)=sump(m)/sumf2(m)   ! 2012_12_13 . scb
      enddo
#endif      
! added end

      if(ng.ne.2) then
         call tpenbcmg_sp3(fnorm,fnorm2)
         return
      endif

!! 2012_12_13 . scb      
      do iz=1,nz
         do ih=1,nassy
            do ig=1,2
               hflx2(1,ig,ih,iz)=phi(ig,ih,iz)+2.*hflx2(2,ig,ih,iz)
            enddo
         enddo
#ifndef sp3_nodal_correction   ! 2014_03_12 . scb               
         do ip=1,nxpnt
            do ig=1,2
               pflx2(1,ig,ip,iz)=pflx2(1,ig,ip,iz)*fnorm(ig)
               !pflx2(2,ig,ip,iz)=pflx2(2,ig,ip,iz)*fnorm2(ig)
               pflx2(2,ig,ip,iz)=pflx2(2,ig,ip,iz)*fnorm(ig)  ! 2013_04_25 . scb
            enddo
         enddo
#endif         
      enddo
!! added end      

! surface flux cal. from phis=w*phin+(1-w)phim+beta*(phim+phin)/2
      do iz=1,nz
         do ih=1,nassy
            do it=1,6
               nn=neignd(it,ih)
               isfc=neigsfc(it,ih)
               do ig=1,2
                  dc=xsd(ig,ih,iz)
                  dfdm=dfd(ig,isfc,iz)
                  dhatm=wtdhat(it,ih)*dhat(ig,isfc,iz)
                  bt=betaphis(ig,isfc,iz)
                  oneflxm=phi(ig,ih,iz)
                  if(nn.ne.0) then
                     w=dfdm/dc*0.5d0

                     oneflxn=phi(ig,nn,iz)
                     phis=w*oneflxn+(1.d0-w)*oneflxm+bt*(oneflxm+oneflxn)*0.5d0
            
                     w=dfd2(ig,isfc,iz)/xsd2(ig,ih,iz)*0.5d0
                     bt=betaphis2(ig,isfc,iz)
                     phis2=w*hflx2(2,ig,nn,iz)+(1.d0-w)*hflx2(2,ig,ih,iz) +bt*(hflx2(2,ig,ih,iz)+hflx2(2,ig,nn,iz))*0.5d0

                     cjn=-((dfdm+dhatm)*oneflxn-(dfdm-dhatm)*oneflxm)/rt3/hside            
                  else
                     cjn=(dfdm-dhatm)*oneflxm/rt3/hside 
                     if(bt.eq.0) then
                        phis=oneflxm/(1.d0+0.25d0*hf2f/dc*alxr)
                        phis2=hflx2(2,ig,ih,iz)/(1.d0+0.25d0*hf2f/xsd2(ig,ih,iz)*alxr)   ! 2012_12_13 . scb
                     else
                        phis=bt*oneflxm

                        bt=betaphis2(ig,isfc,iz)
                        phis2=bt*hflx2(2,ig,ih,iz)
                     endif
                  endif ! nn 
                                
                  !cnto2(1,ig,it,ih,iz)=0.5d0*cjn+0.25d0*phis +5.d0/16.d0*sfli2(2,ig,it,ih,iz)   ! 2012_12_13 . scb
                   
                  cnto2(1,ig,it,ih,iz)=0.5d0*cjn+0.25d0*phis +5.d0/16.d0*phis2
               enddo ! ig
            enddo ! it
         enddo ! ih
      enddo ! iz

      if(nz.eq.1) return   ! 2012_12_18 . scb 
! update z-direction out-current
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
                  flx1=phi(ig,ih,iz1)
                  flx2=phi(ig,ih,iz2)
                  w=dfdm/xsd(ig,ih,iz1)/(1.+hzr)
                  phis=w*flx2+(1.-w)*flx1+bt*(flx1+flx2)*0.5d0

                  cjn=-((dfdm+dhatm)*flx2-(dfdm-dhatm)*flx1)*2.d0 &
                        /(hz(iz1)+hz(iz2))
                         
 
                  w=dfdz2(ig,ih,iz)/xsd2(ig,ih,iz1)/(1.+hzr) 
                  bt=betaphisz2(ig,ih,iz)
                  phis2=w*hflx2(2,ig,ih,iz2)+(1.-w)*hflx2(2,ig,ih,iz1)+bt*(hflx2(2,ig,ih,iz1)+hflx2(2,ig,ih,iz2))*0.5d0
                          
                  !cntzo2(1,ig,id1,ih,iz1)=(phis+2.d0*cjn)/4.d0 +5.d0/16.d0*phisz(2,ig,2,ih,iz1)   ! 2012_12_13 . scb
                  !cntzo2(1,ig,id2,ih,iz2)=(phis-2.d0*cjn)/4.d0 +5.d0/16.d0*phisz(2,ig,1,ih,iz2)   ! 2012_12_13 . scb

                  cntzo2(1,ig,id1,ih,iz1)=(phis+2.d0*cjn)/4.d0 +5.d0/16.d0*phis2
                  cntzo2(1,ig,id2,ih,iz2)=(phis-2.d0*cjn)/4.d0 +5.d0/16.d0*phis2
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
                  flx1=phi(ig,ih,iz1)
                  cjn=(dfdm-dhatm)*flx1/hz(iz1)
                  if(bt.eq.0) then
                     phis=flx1/(1.d0+0.25d0*hz(iz1)/xsd(ig,ih,iz1)*alphaz)
                  else
                     phis=bt*flx1

                     bt=betaphisz2(ig,ih,iz)
                     phis2=bt*hflx2(2,ig,ih,iz1)
                  endif
                  !cntzo2(1,ig,id1,ih,iz1)=(phis+2.d0*cjn)/4.d0 +5.d0/16.d0*phisz(2,ig,2,ih,iz1)   ! 2012_12_13 . scb
                  cntzo2(1,ig,id1,ih,iz1)=(phis+2.d0*cjn)/4.d0 +5.d0/16.d0*phis2
               enddo
            enddo
         endif
      enddo

      if(.not.iftran) return

! 2012_12_05 . scb      
      !if(.not.flagbdf) stop "Time-dependent SP3 TPEN is available only for BDF method !"
      !do k=1,nz
      !  do l=1,nxy
      !    rvol=1/volnode(l,k)
      !    sefve2(1,1,l,k)=srctrf(1,l,k)*rvol-rvelotm(1,l,k)*phi(1,l,k) &
      !                -(1-betap(l,k))*psi(l,k)*rvol
      !    sefve2(1,2,l,k)=srctrf(2,l,k)*rvol-rvelotm(2,l,k)*phi(2,l,k)
      !    
      !    sefve2(2,1,l,k)=-2.d0/3.d0*sefve2(1,1,l,k)
      !    sefve2(2,2,l,k)=-2.d0/3.d0*sefve2(1,2,l,k)
      !  enddo
      !enddo
      !
      !if(flagbdf) then
      !  do k=1,nz
      !    do l=1,nxy
      !      sefve2(2,1,l,k)=sefve2(2,1,l,k)+srctrf2(1,l,k)-5.d0/3.d0*rvelotm(1,l,k)*hflx2(2,1,l,k)
      !      sefve2(2,2,l,k)=sefve2(2,2,l,k)+srctrf2(2,l,k)-5.d0/3.d0*rvelotm(2,l,k)*hflx2(2,2,l,k)
      !    enddo
      !  enddo
      !endif
! added end

      return
    end subroutine

