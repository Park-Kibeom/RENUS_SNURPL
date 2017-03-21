! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    subroutine swpnemz(iftran)

! control NEM Solver for z-direction 

	    use const
	    use allocs
	    use geomhex
	    use tpen
	    use xsec
      use sfam,     only : reigv  

	    implicit none

! variables for nem solver in z-direction
! atleak : transverse leakage for axial direction
! atleaku,atleakl : transverse leakage for upper and lower node
! dcu,dcl : xsd(diffusion coefficients) for upper and lower node
! cntz0i,cntz1i : incoming currents from lower and upper surface 
#define scb_correction_tl
#define paral_p1_nem

#ifdef paral_p1_nem
	    logical :: iftran

      real :: atleakl(ng2),atleaku(ng2)
      real :: dcu(ng2),dcl(ng2)
      real :: cntz0i(ng2),cntz1i(ng2)
      integer :: mv(6)
      data mv/4,5,6,1,2,3/

	    integer :: ih, iz, ig, it, nn, neigdn, neigup
	    real :: sum, adfl, reflratl, atleakm, dcm, hzm, hzu, hzl

      !if(nz.eq.1) return  !if 2D     ! commented by 2012_12_18 . scb
!      if(nassy.eq.1) goto 1000       ! commented by 2012_10_23 . scb

! average transverse leakage cal. for z-direction 
      do ih=1,nassy
!        do iz=ifrom,ito,istp
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
            atleak(ig,ih,iz)=sum*2*rt3/9./hside !1aug00
	          if(iftran) then
	            atleak(ig,ih,iz)=atleak(ig,ih,iz)-sefve(ig,ih,iz)
	          endif
          enddo   
        enddo   
      enddo   
      
1000  continue
      
! solve z-direction using nem solver     
!$OMP PARALLEL DO  private(ih,iz,neigdn,neigup,ig,atleakm,atleakl,atleaku,dcm,dcl,dcu,hzm,hzu,hzl,cntz0i,cntz1i)
      do ih=1,nassy
! ih : index for each hexagonal node 
!        do iz=ifrom,ito,istp
        do iz=1,nz
          neigdn=neigz(1,iz) 
          neigup=neigz(2,iz)
          do ig=1,ng2
            atleakm=atleak(ig,ih,iz) 
            atleakl(ig)=atleakm 
            atleaku(ig)=atleakm 
            dcm=xsd(ig,ih,iz) 
            dcl(ig)=dcm 
            dcu(ig)=dcm
            hzm=hz(iz) 
            hzu=hzm
            hzl=hzm 
            if(neigdn.ne.0) then 
              cntz0i(ig)=cntzo(ig,2,ih,neigdn)
              atleakl(ig)=atleak(ig,ih,neigdn)
              dcl(ig)=xsd(ig,ih,neigdn)
              hzl=hz(neigdn) 
            else
              cntz0i(ig)=reflratzb*cntzo(ig,1,ih,iz)
              if(alzl.ne.0) atleakl(ig)=0 
! 2013_07_10 . scb            
#ifdef scb_correction_tl
              if(alzl.ne.0) hzl=0.d0
#endif              
            endif  
            if(neigup.ne.0) then 
              cntz1i(ig)=cntzo(ig,1,ih,neigup)
              atleaku(ig)=atleak(ig,ih,neigup) 
              dcu(ig)=xsd(ig,ih,neigup) 
              hzu=hz(neigup) 
            else
              cntz1i(ig)=reflratzt*cntzo(ig,2,ih,iz)
              if(alzr.ne.0) atleaku(ig)=0
! 2013_07_10 . scb            
#ifdef scb_correction_tl
              if(alzl.ne.0) hzu=0.d0
#endif              
            endif
	        enddo

          call onenemz                                                     &
            (reigv,cntz0i,cntz1i,atleak(1,ih,iz),atleaku,atleakl,          &
             xsd(1,ih,iz),dcu,dcl,                                         &
             hzm,hzu,hzl,hside,xst(1,ih,iz),xsnf(1,ih,iz),xss(1,2,ih,iz),  &
             hflx(1,ih,iz),                                                &
             cntzo(1,1,ih,iz),cntzo(1,2,ih,iz),srcz(1,ih,iz))

          if(iftran) then
            srcz(1,ih,iz)=srcz(1,ih,iz)+sefve(1,ih,iz)
            srcz(2,ih,iz)=srcz(2,ih,iz)+sefve(2,ih,iz)
	        endif
	      enddo 
	    enddo
!$OMP END PARALLEL DO      
 
      return
    end
#else
	    logical :: iftran

      real :: atleakl(ng2),atleaku(ng2)
      real :: dcu(ng2),dcl(ng2)
      real :: cntz0i(ng2),cntz1i(ng2)
      integer :: mv(6)
      data mv/4,5,6,1,2,3/

	    integer :: ih, iz, ig, it, nn, neigdn, neigup
	    real :: sum, adfl, reflratl, atleakm, dcm, hzm, hzu, hzl

      !if(nz.eq.1) return  !if 2D     ! commented by 2012_12_18 . scb
!      if(nassy.eq.1) goto 1000       ! commented by 2012_10_23 . scb

! average transverse leakage cal. for z-direction 
      do ih=1,nassy
!        do iz=ifrom,ito,istp
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
            atleak(ig,ih,iz)=sum*2*rt3/9./hside !1aug00
	          if(iftran) then
	            atleak(ig,ih,iz)=atleak(ig,ih,iz)-sefve(ig,ih,iz)
	          endif
          enddo   
        enddo   
      enddo   
      
1000  continue
      
! solve z-direction using nem solver     
      do ih=1,nassy
! ih : index for each hexagonal node 
!        do iz=ifrom,ito,istp
        do iz=1,nz
          neigdn=neigz(1,iz) 
          neigup=neigz(2,iz)
          do ig=1,ng2
            atleakm=atleak(ig,ih,iz) 
            atleakl(ig)=atleakm 
            atleaku(ig)=atleakm 
            dcm=xsd(ig,ih,iz) 
            dcl(ig)=dcm 
            dcu(ig)=dcm
            hzm=hz(iz) 
            hzu=hzm
            hzl=hzm 
            if(neigdn.ne.0) then 
              cntz0i(ig)=cntzo(ig,2,ih,neigdn)
              atleakl(ig)=atleak(ig,ih,neigdn)
              dcl(ig)=xsd(ig,ih,neigdn)
              hzl=hz(neigdn) 
            else
              cntz0i(ig)=reflratzb*cntzo(ig,1,ih,iz)
              if(alzl.ne.0) atleakl(ig)=0 
            endif  
            if(neigup.ne.0) then 
              cntz1i(ig)=cntzo(ig,1,ih,neigup)
              atleaku(ig)=atleak(ig,ih,neigup) 
              dcu(ig)=xsd(ig,ih,neigup) 
              hzu=hz(neigup) 
            else
              cntz1i(ig)=reflratzt*cntzo(ig,2,ih,iz)
              if(alzr.ne.0) atleaku(ig)=0
            endif
	        enddo

          call onenemz                                                     &
            (reigv,cntz0i,cntz1i,atleak(1,ih,iz),atleaku,atleakl,          &
             xsd(1,ih,iz),dcu,dcl,                                         &
             hzm,hzu,hzl,hside,xst(1,ih,iz),xsnf(1,ih,iz),xss(1,2,ih,iz),  &
             hflx(1,ih,iz),                                                &
             cntzo(1,1,ih,iz),cntzo(1,2,ih,iz),srcz(1,ih,iz))

          if(iftran) then
            srcz(1,ih,iz)=srcz(1,ih,iz)+sefve(1,ih,iz)
            srcz(2,ih,iz)=srcz(2,ih,iz)+sefve(2,ih,iz)
	        endif
	      enddo 
	    enddo
 
      return
    end
#endif