! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine swptpenmg_senm(errmax,ifrom,ito,istp,ifromz,itoz,istpz)

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
	real :: srczn(ng,ntph),srcbal(ng,ntph),adfl(ng)
	real ::  srcmomx(ng,ntph),srcmomy(ng,ntph) 
	real ::  radfl(ng)
	integer ::  mv(6)
	integer :: it, ig, ih, iz, m, nn
	real :: reflrhom, cntihet, cntohet, srczm, errh
	real :: cntihom, cntohom
	data mv/4,5,6,1,2,3/

	logical, save :: first=.TRUE.

	if(first) then
	xmom=1.
	ymom=1.
	sfli=1.
	cpfl=1.
	first=.FALSE.
	endif

! initialize source
      errmax=0.
      do it=1,6
        do ig=1,ng
          srcbal(ig,it)=0
          srcmomx(ig,it)=0
          srcmomy(ig,it)=0
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
	            reflrhom=(1-2.*alxr)/(1+2.*alxr)
                cnti(ig,it)=cnto(ig,it,ih,iz)*reflrhom
                if(alxr.eq.0) then
                  srczn(ig,it)=srcz(ig,ih,iz)
                else 
                  srczn(ig,it)=0.
				      srczn(ig,it)=-srcz(ig,ih,iz)
                endif
              enddo
            else
              do ig=1,ng
                cntihet=cnto(ig,neigjin(it,ih),nn,iz)
                cntohet=cnto(ig,it,ih,iz) 
                cnti(ig,it)=0.5*((cntohet+cntihet)*radfl(ig) &
                               -(cntohet-cntihet))
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

          if(nz.gt.1) then
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
          endif

! Solve the problem by TPEN
          call onetpenmg_senm                                              &
             (ng,reigv,xstf(:,ih,iz),xsnff(:,ih,iz),xssf(:,:,ih,iz)        &
			    ,xsdf(:,ih,iz),xschif(:,ih,iz)                                &
             ,pbflx,cnti,hside,srcbal,srcmomx,srcmomy                      &
             ,cnto(:,:,ih,iz),hflxf(:,ih,iz),errh                          &
!			 ,aflx(:,:,ih,iz),xmom(:,:,ih,iz),ymom(:,:,ih,iz),cpfl(:,ih,iz),sfli(:,:,ih,iz)  &
             ,sfli(:,:,ih,iz) &
			 )
!
          do it=1,6
             do ig=1,ng
                cntihom=cnti(ig,it)
                cntohom=cnto(ig,it,ih,iz)
                cnto(ig,it,ih,iz)=          &
!     +             0.5*((cntohom+cntihom)*xsadf(ig,ih,iz)
                   0.5*((cntohom+cntihom)   &
                       +(cntohom-cntihom))
             enddo
          enddo
!
          if(errh.gt.errmax) errmax=errh
	  enddo
	enddo

      return
      end
