! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  !
      subroutine swptpenmg(iftran,irfrom,irto,irstp,iafrom,iato,iastp)

! control NEM Solver for z-direction 

	use const
	use geomhex
	use tpen
	use xsec
	use sfam, only : reigv
	implicit none

	integer :: irfrom,irto,irstp,iafrom,iato,iastp
	logical :: iftran

! variables for nem solver in z-direction
! atleak : transverse leakage for axial direction
! atleaku,atleakl : transverse leakage for upper and lower node
! xsdu,xsdl : xsd(diffusion coefficients) for upper and lower node
! cntz0i,cntz1i : incoming currents from lower and upper surface 

      real,allocatable,save,dimension(:) :: atleakl,atleaku, &
                                            xsdu,xsdl,cntz0i,cntz1i,chieff,radfl
      real,allocatable,save,dimension(:,:) :: xssmg,cnti,pbflx,srczn

      integer :: mv(6)
	data mv/4,5,6,1,2,3/
	logical,save :: first=TRUE
      real, save :: areavol

	integer :: ig, ih, iz, it, nn, m, neigdn, neigup, md, jupstt
	real :: sum, adfl, reflratl, rhz, errl2tpen, flxl2tpen, reflrhom, &
	        cntihet, cntohet, hzm, hzu, hzl, atleakm, xsdm,           &
	        cntihom, cntohom


#ifdef DEBUG
	integer :: iscatib(7), iscatie(7)
	data iscatib/1,1,1,3,3,3,3/
	data iscatie/0,1,2,5,7,7,6/
#endif

      if(first) then
        areavol=2.*rt3/9.D0/hside
        allocate(atleakl(ng),atleaku(ng),xsdu(ng),xsdl(ng), &
                 cntz0i(ng),cntz1i(ng),chieff(ng),radfl(ng))
        atleakl=0
        atleaku=0
        xsdu=0
        xsdl=0
        cntz0i=0
        cntz1i=0
        chieff=0
        radfl=0
        allocate(xssmg(ng,ng),cnti(ng,ntph),pbflx(ng,ntph), &
                 srczn(ng,ntph))
        xssmg=0
        cnti=0
        pbflx=0
        srczn=0
        first=.FALSE.
      endif


! average transverse leakage cal. for z-direction 
      do ih=1,nassy
        do iz=1,nz
          do ig=1,ng
            atleak(ig,ih,iz)=-sefve(ig,ih,iz)
            srcz(ig,ih,iz)=sefve(ig,ih,iz)
          enddo   
        enddo
      enddo   
      if(nassy.eq.1) goto 100
      do ih=1,nassy
        do iz=1,nz
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
            atleak(ig,ih,iz)=atleak(ig,ih,iz)+sum*areavol
          enddo   
        enddo   
      enddo   
! Axial Source for TPEN cal.
  100 continue
      if(nz.eq.1) goto 200
      do iz=1,nz
        rhz=1./hz(iz)
        if(iz.eq.1) then
          do ih=1,nassy
            do ig=1,ng
              srcz(ig,ih,iz)=srcz(ig,ih,iz)                         &
                            +((reflratzb-1.)*cntzo(ig,1,ih,iz) &
                            +cntzo(ig,1,ih,iz+1)-cntzo(ig,2,ih,iz))*rhz
            enddo
          enddo

        else
          if(iz.eq.nz) then
            do ih=1,nassy
              do ig=1,ng
                srcz(ig,ih,iz)=srcz(ig,ih,iz)                       &
                            +(cntzo(ig,2,ih,iz-1)-cntzo(ig,1,ih,iz) &
                            +(reflratzt-1.)*cntzo(ig,2,ih,iz))*rhz
              enddo
            enddo
          else
            do ih=1,nassy
              do ig=1,ng
                srcz(ig,ih,iz)=srcz(ig,ih,iz)                         &
                            +(cntzo(ig,2,ih,iz-1)+cntzo(ig,1,ih,iz+1) &
                            -cntzo(ig,1,ih,iz)-cntzo(ig,2,ih,iz))*rhz
              enddo
            enddo
          endif
        endif
      enddo
! solve from irfrom to irto in radial direction
!   and from iafrom to iato in axial direction
  200 continue
      errl2tpen=0
      flxl2tpen=0
      do iz=iafrom,iato,iastp
        do ih=irfrom,irto,irstp
! preparation for radial solver(cnti,srczn,pbflx)
          do m=1,ng
!             radfl(m)=1/xsadf(m,ih,iz)
             radfl(m)=1.
          enddo
          if(nassy.gt.1) then
            do it=1,6
              nn=neignd(it,ih)
              if(nn.eq.0) then
                do ig=1,ng
                  reflrhom=(1-2*alxr)/(1.+2*alxr)
                  cnti(ig,it)=cnto(ig,it,ih,iz)*reflrhom
                  if(ifbcref) then
                    srczn(ig,it)=srcz(ig,ih,iz)
                  else 
                    srczn(ig,it)=0.
					srczn(ig,it)=-srcz(ig,ih,iz) ! fix_tpen
                  endif
                enddo
              else
                do ig=1,ng
! 2014_05_21 . scb adf hex
                  !cntihet=cnto(ig,neigjin(it,ih),nn,iz)
                  !cntohet=cnto(ig,it,ih,iz) 
                  !cnti(ig,it)=0.5*((cntohet+cntihet)*radfl(ig) &
                  !                -(cntohet-cntihet))
                  cnti(ig,it)=cnto(ig,neigjin(it,ih),nn,iz)
                  srczn(ig,it)=srcz(ig,nn,iz)
! added end                  
                enddo
              endif
            enddo
! 2014_05_21 .  scb adf hex
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
                pbflx(ig,it)=pflx(ig,nn,iz)/xscdf(it,ig,ih,iz)   ! 2014_05_21 .  scb adf hex
              enddo
            enddo
          endif
! preparation for axial solver(cnti,srczn,pbflx)
          hzm=hz(iz) 
          hzu=hzm
          hzl=hzm 
          if(nz.gt.1) then
            neigdn=neigz(1,iz) 
            neigup=neigz(2,iz)
            do ig=1,ng
              atleakm=atleak(ig,ih,iz) 
              atleakl(ig)=atleakm 
              atleaku(ig)=atleakm 
              xsdm=xsdf(ig,ih,iz) 
              xsdl(ig)=xsdm 
              xsdu(ig)=xsdm
            enddo
            if(neigdn.ne.0) then 
              do ig=1,ng
                cntz0i(ig)=cntzo(ig,2,ih,neigdn)
                atleakl(ig)=atleak(ig,ih,neigdn)
                xsdl(ig)=xsdf(ig,ih,neigdn)
                hzl=hz(neigdn) 
              enddo
            else
              do ig=1,ng
                cntz0i(ig)=reflratzb*cntzo(ig,1,ih,iz)
                if(.not.ifbcrefzb) atleakl(ig)=0 
              enddo
            endif
            if(neigup.ne.0) then 
              do ig=1,ng
                cntz1i(ig)=cntzo(ig,1,ih,neigup)
                atleaku(ig)=atleak(ig,ih,neigup) 
                xsdu(ig)=xsdf(ig,ih,neigup) 
                hzu=hz(neigup) 
              enddo
            else
              do ig=1,ng
                cntz1i(ig)=reflratzt*cntzo(ig,2,ih,iz)
                if(.not.ifbcrefzt) atleaku(ig)=0
              enddo
            endif
          endif
          do md=2,ng
             do m=1,md-1
                xssmg(m,md)=xssf(m,md,ih,iz)
             enddo
             do m=md+1,ng  !temp
                xssmg(m,md)=xssf(m,md,ih,iz)
             enddo
          enddo
          do m=1,ng
             chieff(m)=xschif(m,ih,iz)
!     +             +betat(ih,iz)*(xschidf(m,ih,iz)-xschif(m,ih,iz))
          enddo
!          jupstt=lupscat
           jupstt=1
!          jupstt=mge(1)
          call onetpenmg(ng,nassy,nz,ih,iz,jupstt, &
!           iscatib,iscatie, &
           xssfs(:,ih,iz), xssfe(:,ih,iz), &
           xstf(1,ih,iz),xsnff(1,ih,iz),xssmg,xsdf(1,ih,iz), &
           chieff,hside,hzm,hzu,hzl, &
           pbflx,cnti,cntz0i,cntz1i,sefve(1,ih,iz),srczn, &
           atleaku,atleakl,reigv, &
           aflx(1,1,ih,iz),xmom(1,1,ih,iz),ymom(1,1,ih,iz), &
           zmom(1,:,ih,iz),zmom(2,:,ih,iz), &
           cnto(1,1,ih,iz),cntzo(1,1,ih,iz),cntzo(1,2,ih,iz), &
           hflxf(1,ih,iz))
           do it=1,6
              do ig=1,ng
                 cntihom=cnti(ig,it)
                 cntohom=cnto(ig,it,ih,iz)
                 cnto(ig,it,ih,iz)=0.5*((cntohom+cntihom)*xsadf(it,ig,ih,iz)+(cntohom-cntihom))  ! 2014_05_21 . scb adf hex
             enddo
           enddo
	   enddo
	enddo

#ifdef DEBUG 
	print '(a10,10F9.5)', 'NEM SrcZ',srcz(1,1,1:nz)
!	print '(a10,10F9.5)', 'Leakage', atleak(1,1,1:10)
#endif

      return
      end
