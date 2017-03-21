! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    subroutine tpenbcmg_sp3(fnorm,fnorm2)
    
#define sp3_nodal_correction   ! 2014_03_12 . scb        

! Update Nodal Solution from CMFD results

	    use const
	    use geomhex
	    use tpen_sp3
	    use xsec,     only : xsd2, xsd, mgb, mge
	    use sfam,     only : phi, phif, fphi, fphi2
      use tpen, only : fohflx, fhflx, focnto, fcnto, fcntzo, focntzo

      use cmfdhex2g, only : dfd, dfd2, dhat, betaphis, betaphis2, dfdz, dhatz, betaphisz, phisz, &
                            dfdz2, betaphisz2


	    implicit none

      real :: fnorm(ng2), fnorm2(ng2)
	    real :: alphaz
	    integer :: iz, ih, ig, ip, igf, ig2, it, nn, isfc, ind, iz1, &
	              iz2, id1, id2, im
	    real :: wf, onewf, rrt3h, ttt, ttt2, dfdm, dhatm, bt, oneflxm, w, &
	            hzr, flx1, flx2, cntoavg,                  &
	            tt1, tt2, oneflxn, ttt3	           
      real :: dc, phis, cjn, phis2
      real :: cntoavg1, cntoavg2

	    integer :: iptype=0
	    logical :: transient ! 2014_03_02 . scb
      
      integer,save :: istep=0   ! 2014_02_27 . scb for dbg
      
      common / sp3dbg / transient ! 2014_03_02 . scb
      
      istep=istep+1   ! 2014_02_27 . scb
      

! 2014_02_27 . scb : comment here to iptype if statement      
! iptype  
! 0 : eigenvlaue problem
! ! : fixed source problem
#ifndef sp3_nodal_correction   ! 2014_03_12 . scb    
      do iz=1,nz
         do ih=1,nassy
            do ig=1,ng2        
               hflx2(2,ig,ih,iz)=0.
               do ig2=mgb(ig),mge(ig)
                  hflx2(2,ig,ih,iz)=hflx2(2,ig,ih,iz)+hflxf2(2,ig2,ih,iz)
               enddo ! ig2
               hflx2(1,ig,ih,iz)=phi(ig,ih,iz) 
            enddo
         enddo
         do ip=1,nxpnt
            do ig=1,ng2
               do ig2=mgb(ig),mge(ig)
                  pflx2(1,ig2,ip,iz)=pflx2(1,ig2,ip,iz)*fnorm(ig)   
               enddo
            enddo
         enddo
      enddo ! iz
      
      if(iptype.eq.0) then
         wf=1.5
      else
         wf=1.
      endif
#else
      wf=1.d0  ! 2014_02_26 . scb    
#endif      

      do iz=1,nz
         do ih=1,nassy
! 2014_03_12 . scb           
            do ig=1,2
               hflx2(1,ig,ih,iz)=phi(ig,ih,iz)
            enddo
! added end
            do ig=1,ng
               hflxf2(1,ig,ih,iz)=phif(ig,ih,iz)+2.*hflxf2(2,ig,ih,iz)
            enddo
         enddo
      enddo      
! added end   
   

      !if(istep.gt.3) return   ! 2014_02_27 . scb
      
      !if(transient) return   ! 2014_03_02 . scb for dbg
            
      onewf=1.-wf
      rrt3h=1./rt3/hside
      
! 2014_02_26 . scb -> remove below lines... 
#ifndef sp3_nodal_correction   ! 2014_03_12 . scb    
      do iz=1,nz
         do ih=1,nassy
            do ig=1,2
               do ig2=mgb(ig),mge(ig)
                  ttt=wf*fhflx(ig2,ih,iz)+(1.-wf)*fohflx(ig2,ih,iz)
                  ttt=fphi(ig2,ih,iz)
                  ttt2=ttt*hflx2(1,ig,ih,iz)/hflxf2(1,ig2,ih,iz)
                  hflxf2(1,ig2,ih,iz)=ttt*hflx2(1,ig,ih,iz)
                  phif(ig2,ih,iz)=hflxf2(1,ig2,ih,iz)
      
                  ttt=fphi2(ig2,ih,iz)
                  ttt3=ttt*hflx2(2,ig,ih,iz)/hflxf2(2,ig2,ih,iz)
                  hflxf2(2,ig2,ih,iz)=ttt*hflx2(2,ig,ih,iz)
                  do it=1,6
                     aflx2(1,ig2,it,ih,iz)=aflx2(1,ig2,it,ih,iz)*ttt2
                     xmom2(1,ig2,it,ih,iz)=xmom2(1,ig2,it,ih,iz)*ttt2
                     ymom2(1,ig2,it,ih,iz)=ymom2(1,ig2,it,ih,iz)*ttt2
      
                     aflx2(2,ig2,it,ih,iz)=aflx2(2,ig2,it,ih,iz)*ttt3
                     xmom2(2,ig2,it,ih,iz)=xmom2(2,ig2,it,ih,iz)*ttt3
                     ymom2(2,ig2,it,ih,iz)=ymom2(2,ig2,it,ih,iz)*ttt3
                  enddo ! it
                  zmom1(1,ig2,ih,iz)=zmom1(1,ig2,ih,iz)*ttt2
                  zmom2(1,ig2,ih,iz)=zmom2(1,ig2,ih,iz)*ttt2
      
                  zmom1(2,ig2,ih,iz)=zmom1(2,ig2,ih,iz)*ttt3
                  zmom2(2,ig2,ih,iz)=zmom2(2,ig2,ih,iz)*ttt3
               enddo ! ig2
            enddo ! ig
         enddo ! ih
      enddo ! iz
#endif      
! added end      
      
      !if(iftran) 
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
                  oneflxm=hflx2(1,ig,ih,iz)
                  if(nn.ne.0) then
                     w=dfdm/dc/2.
                     oneflxn=hflx2(1,ig,nn,iz) 
                     phis=w*oneflxn+(1.-w)*oneflxm+0.5*bt*(oneflxm+oneflxn)

                     w=dfd2(ig,isfc,iz)/xsd2(ig,ih,iz)/2.
                     bt=betaphis2(ig,isfc,iz)
                     phis2=w*hflx2(2,ig,nn,iz)+(1.-w)*hflx2(2,ig,ih,iz) +bt*(hflx2(2,ig,ih,iz)+hflx2(2,ig,nn,iz))/2.

                     cjn=-((dfdm+dhatm)*oneflxn-(dfdm-dhatm)*oneflxm)/rt3/hside                                              
                  else
                     cjn=(dfdm-dhatm)*oneflxm/rt3/hside 
                     if(bt.eq.0) then
                        phis=oneflxm/(1+0.25*hf2f/dc*alxr)
                     else
                        phis=bt*oneflxm

                        bt=betaphis2(ig,isfc,iz)
                        phis2=bt*hflx2(2,ig,ih,iz)
                     endif
                  endif ! nn

#ifndef DEBUG
                  phis2=0.
                  do ig2=mgb(ig),mge(ig)
                     phis2=phis2+sfli2(2,ig2,it,ih,iz)  
                  enddo ! ig2
#endif

                  cntoavg=0.5*cjn+0.25*phis +5./16.*phis2
                  do ig2=mgb(ig),mge(ig)
                     ttt=wf*fcnto(ig2,it,ih,iz)+(1.-wf)*focnto(ig2,it,ih,iz)
                     cnto2(1,ig2,it,ih,iz)=ttt*cntoavg 
                  enddo ! ig2
               enddo ! ig
            enddo ! it
         enddo ! ih
      enddo ! iz
!
      if(nz.eq.1) return
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
                  flx1=hflx2(1,ig,ih,iz1)
                  flx2=hflx2(1,ig,ih,iz2)
                  w=dfdm/xsd(ig,ih,iz1)/(1.+hzr)
                  phis=w*flx2+(1.-w)*flx1+0.5*bt*(flx1+flx2)
                  cjn=-((dfdm+dhatm)*flx2-(dfdm-dhatm)*flx1)*2. &
                        /(hz(iz1)+hz(iz2))

                  w=dfdz2(ig,ih,iz)/xsd2(ig,ih,iz1)/(1.+hzr) 
                  bt=betaphisz2(ig,ih,iz)
                  phis2=w*hflx2(2,ig,ih,iz2)+(1.-w)*hflx2(2,ig,ih,iz1)+bt*(hflx2(2,ig,ih,iz1)+hflx2(2,ig,ih,iz2))/2.  
                  
#ifndef DEBUG
                  phis2=0.
                  do ig2=mgb(ig),mge(ig)
                     phis2=phis2+phisz(2,ig2,id1,ih,iz1)
                  enddo ! ig2
#endif
                                                                             
                  cntoavg1=0.25*(phis+2.*cjn) +5./16.*phis2 
                  cntoavg2=0.25*(phis-2.*cjn) +5./16.*phis2 
                  do ig2=mgb(ig),mge(ig)
                     tt1=wf*fcntzo(ig2,id1,ih,iz1) &
                              +(1.-wf)*focntzo(ig2,id1,ih,iz1)
                     tt2=wf*fcntzo(ig2,id2,ih,iz2) &
                              +(1.-wf)*focntzo(ig2,id2,ih,iz2)
                     cntzo2(1,ig2,id1,ih,iz1)=tt1*cntoavg1
                     cntzo2(1,ig2,id2,ih,iz2)=tt2*cntoavg2
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
                  flx1=hflx2(1,ig,ih,iz1)
                  cjn=(dfdm-dhatm)*flx1/hz(iz1)
                  if(bt.eq.0) then
                     phis=flx1/(1+0.25*hz(iz1)/xsd(ig,ih,iz1)*alphaz)
                  else
                     phis=bt*flx1

                     bt=betaphisz2(ig,ih,iz)
                     flx1=hflx2(2,ig,ih,iz1)
                     phis2=bt*flx1
                  endif

#ifndef DEBUG
                  phis2=0.
                  do ig2=mgb(ig),mge(ig)
                     phis2=phis2+phisz(2,ig2,id1,ih,iz1)  
                  enddo ! ig2
#endif

                  cntoavg1=0.25*(phis+2.*cjn) +5./16.*phis2
                  do ig2=mgb(ig),mge(ig)
                     tt1=wf*fcntzo(ig2,id1,ih,iz1) &
                           +(1.-wf)*focntzo(ig2,id1,ih,iz1)
                     cntzo2(1,ig2,id1,ih,iz1)=tt1*cntoavg1 
                  enddo
               enddo
            enddo
         endif
      enddo

      return
    end subroutine

