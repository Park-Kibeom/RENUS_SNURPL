! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine swpsenmz(ifrom,ito,istp)
!
! control SENM Solver for z-direction 
!
	use param
	use allocs
    use sfam,     only : reigv, psi, jnet, phisfc, phif, &
                         curil,curir,curol,curor, resetsfam
	use senm2n,   only : drivesenm2n
	use anm2n,    only : driveanm2n
	use nodal,    only : updjpart, drivenodal
    use cmfdmg,   only : dhatzf, dtilzf 
	use cmfd2g,   only : reigvs
	use nodal,    only : avgjnet_h                   

    include 'global.h'
	include 'geom.h'
	include 'geomh.h'
	include 'xsec.h'
	include 'defhex.h'
	include 'defpnt.h'
	include 'lscoefh.h'
	include 'deffg.h'
	include 'defsfc.h'

	include 'dummy.h'
!
      if(nz.eq.1) return  !if 2D
!
! average transverse leakage cal. for z-direction 
      do ih=1,nassy
        do iz=1,nz
          do ig=1,ng2
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
	      avgjnet_h(ig,ih,iz)=sum*2*rt3/9./hside
          enddo   
        enddo   
      enddo   
!
! variables for nem solver in z-direction
! atleak : transverse leakage for axial direction
! atleaku,atleakl : transverse leakage for upper and lower node
! dcu,dcl : xsd(diffusion coefficients) for upper and lower node
! cntz0i,cntz1i : incoming currents from lower and upper surface 
!
      if(nz.eq.1) return  !if 2D
!
	dhatzf=dhatz
	dtilzf=dfdz 

	do ih=1,nassy
		do iz=1,nz
			sum=0.
			do ig=1,ng2
				phif(ig,ih,iz)=hflx(ig,ih,iz)
				sum=sum+xsnff(ig,ih,iz)*phif(ig,ih,iz)
			enddo
!			psi(ih,iz)=sum*volnode(ih,iz)
		enddo
	enddo

	call drivesenm2n(FALSE,reigv,phif,psi,jnet,phisfc)
!	call driveanm2n(FALSE,reigv,phif,psi,jnet,phisfc)	 
	do ih=1,nassy
		do iz=ifrom,ito,istp
			do ig=1,ng2
				srcz(ig,ih,iz)=1./hz(iz)*(jnet(LEFT,ig,ih,iz,ZDIR)-jnet(RIGHT,ig,ih,iz,ZDIR))
			enddo
		enddo
	enddo

#ifndef DEBUG
	print '(a10,10F9.5)', 'SENM SrcZ',srcz(1,1,1:10)
#endif

	call updjpart(jnet,phisfc)
	do iz=ifrom,ito,istp
		do ih=1,nassy
			do ig=1,ng2
				cntzo(ig,1,ih,iz)=curol(3,ih,iz,ig)
				cntzo(ig,2,ih,iz)=curor(3,ih,iz,ig)
			enddo
		enddo
	enddo


      return
      end

