! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine tpenbcmg(fnorm)

#define p1_nodal_correction   ! 2014_03_12 . scb

! Update Nodal Solution from CMFD results

	use const
	use geomhex
	use tpen
	use xsec
	use sfam,     only : phi, phif
   use cmfdhex2g, only : dfd, dhat, betaphis, dfdz, dhatz, betaphisz

	implicit none

    real :: fnorm(ng2)
	real :: alphaz
	integer :: iz, ih, ig, ip, igf, ig2, it, nn, isfc, ind, iz1, &
	           iz2, id1, id2
	real :: wf, onewf, rrt3h, ttt, ttt2, dfdm, dhatm, bt, oneflxm, w, &
	        phis, cjn, cntoavg, hzr, flx1, flx2, cntoavg1, cntoavg2,  &
	        tt1, tt2, oneflxn

	integer :: iptype=0
	logical :: transient


! iptype  
! 0 : eigenvlaue problem
! ! : fixed source problem

      do iz=1,nz
        do ih=1,nassy
          do ig=1,ng2
            hflx(ig,ih,iz)=phi(ig,ih,iz)
          enddo
          
! 2014_03_12 . scb          
#ifdef p1_nodal_correction   ! 2014_03_12 . scb
          do ig=1,ng
            hflxf(ig,ih,iz)=phif(ig,ih,iz)
          enddo
#endif          
! added end
        enddo
        
#ifndef p1_nodal_correction   ! 2014_03_12 . scb
        do ip=1,nxpnt
          do ig=1,ng2
            do igf=mgb(ig),mge(ig)
              pflx(igf,ip,iz)=pflx(igf,ip,iz)*fnorm(ig)  
            enddo
          enddo
        enddo
#endif        
      enddo

#ifdef p1_nodal_correction   ! 2014_03_12 . scb
      wf=1.d0
#else
! 2014_02_26 . scb -> remove below lines...
      if(iptype.eq.0) then
         wf=1.15D0
      else
         wf=1
      endif      
#endif      

      onewf=1.-wf
      rrt3h=1./rt3/hside
            
#ifndef p1_nodal_correction   ! 2014_03_12 . scb
      do iz=1,nz
        do ih=1,nassy
          do ig=1,2
            do ig2=mgb(ig),mge(ig)
              ttt=wf*fhflx(ig2,ih,iz)+(1.-wf)*fohflx(ig2,ih,iz)
              ttt2=ttt*hflx(ig,ih,iz)/hflxf(ig2,ih,iz)
              hflxf(ig2,ih,iz)=ttt*hflx(ig,ih,iz)
              phif(ig2,ih,iz)=hflxf(ig2,ih,iz)
              do it=1,6
                aflx(ig2,it,ih,iz)=aflx(ig2,it,ih,iz)*ttt2
                xmom(ig2,it,ih,iz)=xmom(ig2,it,ih,iz)*ttt2
                ymom(ig2,it,ih,iz)=ymom(ig2,it,ih,iz)*ttt2
              enddo
              zmom(1,ig2,ih,iz)=zmom(1,ig2,ih,iz)*ttt2
              zmom(2,ig2,ih,iz)=zmom(2,ig2,ih,iz)*ttt2
            enddo
          enddo
        enddo
      enddo
#endif      

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
                phis=w*oneflxn+(1.-w)*oneflxm+0.5*bt*(oneflxm+oneflxn)
                cjn=-((dfdm+dhatm)*oneflxn-(dfdm-dhatm)*oneflxm) &
                      /rt3/hside 
                if(bt.eq.0) then
!                   cntoavg=0.25*(phis*xsadf(mgb(ig),ih,iz)+2.*cjn)
                   cntoavg=0.25*(phis+2.*cjn)
                else
                   cntoavg=0.25*(phis+2.*cjn)
                endif
              else
                cjn=(dfdm-dhatm)*oneflxm/rt3/hside 
                if(bt.eq.0) then
                   phis=oneflxm/(1+0.25*hf2f/xsd(ig,ih,iz)*alxr)
!     +                 *xsadf(mgb(ig),ih,iz)
                else
                   phis=bt*oneflxm
                endif
                cntoavg=0.25*(phis+2.*cjn)
              endif
              do ig2=mgb(ig),mge(ig)
                ttt=wf*fcnto(ig2,it,ih,iz)+(1.-wf)*focnto(ig2,it,ih,iz)
                cnto(ig2,it,ih,iz)=ttt*cntoavg
              enddo
            enddo
          enddo
        enddo
      enddo
! update z-direction out-current
      if(nz.gt.1) then
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
                phis=w*flx2+(1.-w)*flx1+0.5*bt*(flx1+flx2)
                cjn=-((dfdm+dhatm)*flx2-(dfdm-dhatm)*flx1)*2. &
                      /(hz(iz1)+hz(iz2)) 
                cntoavg1=0.25*(phis+2.*cjn)
                cntoavg2=0.25*(phis-2.*cjn)
                do ig2=mgb(ig),mge(ig)
                  tt1=wf*fcntzo(ig2,id1,ih,iz1) &
                     +(1.-wf)*focntzo(ig2,id1,ih,iz1)
                  tt2=wf*fcntzo(ig2,id2,ih,iz2) &
                     +(1.-wf)*focntzo(ig2,id2,ih,iz2)
                  cntzo(ig2,id1,ih,iz1)=tt1*cntoavg1
                  cntzo(ig2,id2,ih,iz2)=tt2*cntoavg2
                enddo
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
                cntoavg1=0.25*(phis+2.*cjn)
                do ig2=mgb(ig),mge(ig)
                  tt1=wf*fcntzo(ig2,id1,ih,iz1) &
                     +(1.-wf)*focntzo(ig2,id1,ih,iz1)
                  cntzo(ig2,id1,ih,iz1)=tt1*cntoavg1
                enddo
              enddo
            enddo
          endif
        enddo
      endif
!
      if(.not.(iptype.eq.1 .and. transient)) return

#ifdef FIXED_SOURCE
! calculate effective trasient fixed source
      do k=1,nz
         do l=1,nxy
            rvol=1/volnode(l,k)
            betamw=1-betap(l,k)
            do m=1,ng
!dj+ : nemmgtr
!o             sefve(m,l,k)=aphif(m,l,k)*rvol
               sefve(m,l,k)=srcf(m,l,k)*rvol
!dj- : nemmgtr
     +                     -rvdeltf(m,l,k)*hflxf(m,l,k)  
     +                   -betamw*xschidf(m,l,k)*psi(l,k)*rvol
            enddo            
         enddo
      enddo
#endif
!
      return
      end

