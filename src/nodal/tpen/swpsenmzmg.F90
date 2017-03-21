! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine swpsenmzmg(ifupd,iftran)

! control SENM Solver for z-direction 

	
	use const
	use allocs
	use geomhex
	use tpen
	use xsec
    use sfam,     only : reigv, psi, jnet, phisfc, phif, &
                           curil,curir,curol,curor             
	use nodal,    only : updjpart, avgjnet_h 
	use senm2n,   only : drivesenm2n 
	use senm2n,   only : resetsenm2n, initsenm2n
	use anm2n,    only : driveanm2n 
	use cmfdmg,   only : upddhatmg

	implicit none

	integer :: ih, iz, ig, it, nn
	real :: sum, adfl, reflratl, sumpsi
	integer, save :: isenm=0
	logical :: ifupd, iftran

	logical, save :: first=TRUE
	real,pointer :: form(:,:,:)
	real :: fnorm(2), sumf(2), sump(2), vol, rhz
	integer:: m, l, k, ip

	real :: time0, time1
	real, save :: time_senm=0.

	if(first) then
		call dmalloc(form,ng,nxy,nz)
		first=FALSE
	endif

      if(nz.eq.1) return  !if 2D

	isenm=isenm+1
! average transverse leakage cal. for z-direction 
      do iz=1,nz
	  do ih=1,nassy
	    sumpsi=0
          do ig=1,ng
            sum=0
            do it=1,6
              nn=neignd(it,ih)
              if(nn.eq.0) then
!                adfl=xsadf(ig,ih,iz)
                adfl=1.
                reflratl=(adfl-2*alxr)/(adfl+2*alxr)
                sum=sum+(1.-reflratl)*cnto(ig,it,ih,iz)
              else
                sum=sum+cnto(ig,it,ih,iz)-cnto(ig,neigjin(it,ih),nn,iz)
              endif
            enddo
!            atleak(ig,ih,iz)=sum*2*rt3/9./hside-sefve(ig,ih,iz)  !1aug00
	      avgjnet_h(ig,ih,iz)=sum*twort3o9h
!		  if(iftran) avgjnet_h(ig,ih,iz)=avgjnet_h(ig,ih,iz)-sefve(ig,ih,iz)
!		  sumpsi=sumpsi+xsnff(ig,ih,iz)*hflx(ig,ih,iz)			
		  phif(ig,ih,iz)=hflxf(ig,ih,iz)
          enddo 
!		if(iftran) psi(ih,iz)=sumpsi*volnode(ih,iz)  
        enddo   
      enddo   
	
      if(nz.eq.1) return  !if 2D

	call cpu_time(time0)
	call drivesenm2n(iftran,reigv,phif,psi,jnet,phisfc)
!	call driveanm2n(iftran,reigv,phif,psi,jnet,phisfc)
	call cpu_time(time1)
	time_senm=time_senm+time1-time0

!	print *, 'SENM', time_senm
	do iz=1,nz
		rhz=1./hz(iz)	 
		do ih=1,nassy
			do ig=1,ng
				srcz(ig,ih,iz)=rhz*(jnet(LEFT,ig,ih,iz,ZDIR)-jnet(RIGHT,ig,ih,iz,ZDIR))
				if(iftran) srcz(ig,ih,iz)=srcz(ig,ih,iz)+sefve(ig,ih,iz)
			enddo
		enddo
	enddo
	if(.not.ifupd) return

#ifdef DEBUG
!	print '(a10,10(1pe12.2))', 'SENM SrcZ',srcz(1,1,1:nz)
!	print '(a10,10(1pe15.5))', 'SENM jout',cnto(1,1,1,1:nz)
!	print '(a10,10(1pe15.5))', 'SENM Phif',phif(1,1,1:nz)
!	print '(a10,10(1pe15.5))', '         ',srcz(2,1,1:nz)
	print '(a10,10(1pe15.5))', 'SENM Jnet',jnet(1,1,1,1:nz,3), jnet(2,1,1,nz,3)
!	print '(a10,10(1pe15.5))', 'SENM JoutZ',curol(3,1,1:nz,1)
!	print '(a10,10(1pe15.5))', '         ',jnet(1,2,1,1:nz,3), jnet(2,2,1,nz,3)
!	print '(a10,10F9.5)', 'Leakage', avgjnet_h(1,1,1:nz)
#endif

	call updjpart(jnet,phisfc)
	do iz=1,nz
		do ih=1,nassy
			do ig=1,ng
				cntzo(ig,1,ih,iz)=curol(3,ih,iz,ig)
				cntzo(ig,2,ih,iz)=curor(3,ih,iz,ig)
			enddo
		enddo
	enddo

      return
      end

