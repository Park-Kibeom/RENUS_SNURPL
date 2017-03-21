! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    subroutine swptpen(errmax,ifrom,ito,istp,ifromz,itoz,istpz)
 
#define paral_p1_tpen

#ifdef paral_p1_tpen
	    use const
	    use geomhex
	    use tpen
	    use xsec
	    use sfam, only : reigv
	    implicit none

	    real :: errmax
	    integer :: ifrom,ito,istp,ifromz,itoz,istpz

! variables for hopen solver in r-direction
	    real :: cnti(ng,ntph),pbflx(ng,ntph)

! variables for nem solver in z-direction
! atleak : transverse leakage for axial direction
! atleaku,atleakl : transverse leakage for upper and lower node
! dcu,dcl : xsd(diffusion coefficients) for upper and lower node
! cntz0i,cntz1i : incoming currents from lower and upper surface 
! 2013_06_11 . scb for parallel
      real :: srczn(ng,ntph),srcbal(ng,ntph),adfl(ng)
      real ::  srcmomx(ng,ntph),srcmomy(ng,ntph) 
! added end      
      real ::  radfl(ng), bsflx(2,6), refl0, adfl0
      integer ::  mv(6)
	    integer :: it, ig, ih, iz, m, nn
	    real :: reflrhom, cntihet, cntohet, srczm, errh
	    real :: cntihom, cntohom
      data mv/4,5,6,1,2,3/

! initialize source
      errmax=0.d0
      do it=1,6 
        do ig=1,ng
          srcbal(ig,it)=0.d0
          srcmomx(ig,it)=0.d0
          srcmomy(ig,it)=0.d0
        enddo
      enddo

! prepare boundary condition(incoming J(cnti) and point flux(pbflx))
!$OMP PARALLEL DO  private(iz,ih,m,radfl,it,nn,ig,cnti,srczn,cntihet,cntohet,pbflx,srczm,srcbal,srcmomx,srcmomy)    
	    do iz=ifromz,itoz,istpz
        do ih=ifrom,ito,istp
          do m=1,ng
!             radfl(m)=1/xsadf(m,ih,iz)
            radfl(m)=1.
          enddo
          do it=1,6
            nn=neignd(it,ih)
            if(nn.eq.0) then
              do ig=1,ng
!                reflrhom=(1-2*alxr(ig))/(xsadf(ig,ih,iz)+2*alxr(ig))
	              !reflrhom=(1-2*alxr)/(1+2*alxr)
                cnti(ig,it)=cnto(ig,it,ih,iz)*(1-2*alxr)/(1+2*alxr)
                if(alxr.eq.0) then
                  srczn(ig,it)=srcz(ig,ih,iz)
                else 
                  srczn(ig,it)=0.
				          !srczn(ig,it)=-srcz(ig,ih,iz) ! fix_tpen     ! commented by 2012_10_23 . scb
                endif
              enddo
            else
              do ig=1,ng
! 2014_05_21 . scb adf hex                
                !cntihet=cnto(ig,neigjin(it,ih),nn,iz)
                !cntohet=cnto(ig,it,ih,iz) 
                !cnti(ig,it)=0.5*((cntohet+cntihet)*radfl(ig)-(cntohet-cntihet))
                !srczn(ig,it)=srcz(ig,nn,iz)
                cnti(ig,it)=cnto(ig,neigjin(it,ih),nn,iz)
                srczn(ig,it)=srcz(ig,nn,iz)
! added end                
              enddo
            endif
          enddo

! 2014_05_21 . scb adf hex          
          do it=1,6
            do ig=1,ng
              cntihet=cnti(ig,it)
              cntohet=cnto(ig,it,ih,iz) 
              cnti(ig,it)=0.5*((cntohet+cntihet)/xsadf(it,ig,ih,iz)-(cntohet-cntihet))
            enddo
          enddo
! added end          
          
          do it=1,6
            nn=neigpt(it,ih)
            do ig=1,ng
              !pbflx(ig,it)=pflx(ig,nn,iz)*radfl(ig)
              pbflx(ig,it)=pflx(ig,nn,iz)/xscdf(it,ig,ih,iz)   ! 2014_05_21 . scb adf hex
            enddo
          enddo

          !if(nz.gt.1) then   ! 2012_12_18 . scb
            do it=1,6
              do ig=1,ng
                srczm=srcz(ig,ih,iz) 
                srcbal(ig,it)=srczm                                 &
                      +(83.*srczn(ig,it)+17.*srczn(ig,mp1(it))        &
                      -37.*srczn(ig,mp2(it))-43.*srczn(ig,mp3(it))   &
                      -37.*srczn(ig,mp4(it))+17.*srczn(ig,mp5(it)))/540. 
                srcmomx(ig,it)=(-60.*srczm                          &
                      +59.*srczn(ig,it)+14.*srczn(ig,mp1(it))        &
                      -10.*srczn(ig,mp2(it))-7.*srczn(ig,mp3(it))    & 
                      -10.*srczn(ig,mp4(it))+14.*srczn(ig,mp5(it)))/3240. 
                srcmomy(ig,it)=-(srczn(ig,mp1(it))-srczn(ig,mp5(it)))/40.  &  !m-g
                        -(srczn(ig,mp2(it))-srczn(ig,mp4(it)))/360.          
              enddo
            enddo
          !endif   ! 2012_12_18 . scb

! Solve the problem by TPEN
          
!#define SC_dbg
#ifndef SC_dbg
          call onetpen                                                     &
             (reigv,xst(1,ih,iz),xsnf(1,ih,iz),xss(1,2,ih,iz),xsd(1,ih,iz) &
             ,pbflx,cnti,hside,srcbal,srcmomx,srcmomy                      &
             ,cnto(1,1,ih,iz),hflx(1,ih,iz),errh)
#endif
#ifdef SC_dbg
! 2015_02_11 . scb for SC input & output          
          call onetpen_dbg                                                 &
             (ih,iz,reigv,xst(1,ih,iz),xsnf(1,ih,iz),xss(1,2,ih,iz),xsd(1,ih,iz) &
             ,pbflx,cnti,hside,srcbal,srcmomx,srcmomy                      &
             ,cnto(1,1,ih,iz),hflx(1,ih,iz),errh,xskp(1,ih,iz))
#endif          
          
! 2014_05_21 . scb adf hex          
          do it=1,6
            do ig=1,ng
              cntihom=cnti(ig,it)
              cntohom=cnto(ig,it,ih,iz) 
              cnto(ig,it,ih,iz)=0.5*((cntohom+cntihom)*xsadf(it,ig,ih,iz)+(cntohom-cntihom))
            enddo
          enddo
! added end                    

	      enddo
      enddo
!$OMP END PARALLEL DO      
      
      return
    end

    
#else

	    use const
	    use geomhex
	    use tpen
	    use xsec
	    use sfam, only : reigv
	    implicit none

	    real :: errmax
	    integer :: ifrom,ito,istp,ifromz,itoz,istpz

! variables for hopen solver in r-direction
	    real :: cnti(ng,ntph),pbflx(ng,ntph)

! variables for nem solver in z-direction
! atleak : transverse leakage for axial direction
! atleaku,atleakl : transverse leakage for upper and lower node
! dcu,dcl : xsd(diffusion coefficients) for upper and lower node
! cntz0i,cntz1i : incoming currents from lower and upper surface 
! 2013_06_11 . scb for parallel
      real :: srczn(ng,ntph),srcbal(ng,ntph),adfl(ng)
      real ::  srcmomx(ng,ntph),srcmomy(ng,ntph) 
! added end      
      real ::  radfl(ng), bsflx(2,6), refl0, adfl0
      integer ::  mv(6)
	    integer :: it, ig, ih, iz, m, nn
	    real :: reflrhom, cntihet, cntohet, srczm, errh
	    real :: cntihom, cntohom
      data mv/4,5,6,1,2,3/

! initialize source
      errmax=0.d0
      do it=1,6 
        do ig=1,ng
          srcbal(ig,it)=0.d0
          srcmomx(ig,it)=0.d0
          srcmomy(ig,it)=0.d0
        enddo
      enddo

! prepare boundary condition(incoming J(cnti) and point flux(pbflx))
	    do iz=ifromz,itoz,istpz
        do ih=ifrom,ito,istp
          do m=1,ng
!             radfl(m)=1/xsadf(m,ih,iz)
            radfl(m)=1.
          enddo
          do it=1,6
            nn=neignd(it,ih)
            if(nn.eq.0) then
              do ig=1,ng
!                reflrhom=(1-2*alxr(ig))/(xsadf(ig,ih,iz)+2*alxr(ig))
	              !reflrhom=(1-2*alxr)/(1+2*alxr)
                cnti(ig,it)=cnto(ig,it,ih,iz)*(1-2*alxr)/(1+2*alxr)
                if(alxr.eq.0) then
                  srczn(ig,it)=srcz(ig,ih,iz)
                else 
                  srczn(ig,it)=0.
				          !srczn(ig,it)=-srcz(ig,ih,iz) ! fix_tpen     ! commented by 2012_10_23 . scb
                endif
              enddo
            else
              do ig=1,ng
                cntihet=cnto(ig,neigjin(it,ih),nn,iz)
                cntohet=cnto(ig,it,ih,iz) 
                cnti(ig,it)=0.5*((cntohet+cntihet)*radfl(ig)-(cntohet-cntihet))
                srczn(ig,it)=srcz(ig,nn,iz)
              enddo
            endif
          enddo

          do it=1,6
            nn=neigpt(it,ih)
            do ig=1,ng
              pbflx(ig,it)=pflx(ig,nn,iz)*radfl(ig)
            enddo
          enddo

          !if(nz.gt.1) then   ! 2012_12_18 . scb
            do it=1,6
              do ig=1,ng
                srczm=srcz(ig,ih,iz) 
                srcbal(ig,it)=srczm                                 &
                      +(83.*srczn(ig,it)+17.*srczn(ig,mp1(it))        &
                      -37.*srczn(ig,mp2(it))-43.*srczn(ig,mp3(it))   &
                      -37.*srczn(ig,mp4(it))+17.*srczn(ig,mp5(it)))/540. 
                srcmomx(ig,it)=(-60.*srczm                          &
                      +59.*srczn(ig,it)+14.*srczn(ig,mp1(it))        &
                      -10.*srczn(ig,mp2(it))-7.*srczn(ig,mp3(it))    & 
                      -10.*srczn(ig,mp4(it))+14.*srczn(ig,mp5(it)))/3240. 
                srcmomy(ig,it)=-(srczn(ig,mp1(it))-srczn(ig,mp5(it)))/40.  &  !m-g
                        -(srczn(ig,mp2(it))-srczn(ig,mp4(it)))/360.          
              enddo
            enddo
          !endif   ! 2012_12_18 . scb

! Solve the problem by TPEN

          call onetpen                                                     &
             (reigv,xst(1,ih,iz),xsnf(1,ih,iz),xss(1,2,ih,iz),xsd(1,ih,iz) &
             ,pbflx,cnti,hside,srcbal,srcmomx,srcmomy                      &
             ,cnto(1,1,ih,iz),hflx(1,ih,iz),errh)

	      enddo
      enddo
      
      return
    end
#endif      