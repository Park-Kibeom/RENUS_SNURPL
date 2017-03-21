    subroutine initth
! initialize t/h condition and calculation parameters
      use param
      use allocs
!      
      include 'global.h'
      include 'geom.h'
      include 'xsec.h'
      include 'pow.h'
      include 'thop.inc'
      include 'thgeom.inc'
      include 'thfuel.inc'
      include 'frzfdbk.inc'
      include 'thcntl.inc'
      include 'thfdbk.inc'
      include 'thcool.inc'
!      include 'mslb.inc'  ! MSLB
      
      integer :: meshx(0:nxa), meshy(0:nya)

      thetafb=1-thetaf
      thetacb=1-thetac
      rthetac=thetacb/thetac
            
!! count the number of channel
!      nchan=0
!      do ja=1,nya
!        nthy1=nthy(ja)
!        nthx1=0
!        do ia=nxsfa(ja),nxefa(ja)
!          nthx1=nthx1+nthx(ia)
!        enddo
!        nchan=nchan+nthx1*nthy1
!      enddo          
!
!! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
!! determine nchan considering control assemblies
!      if(iflfr .or. .not.rect) then    ! FAST REACTOR
!        ichan=0                         ! FAST REACTOR
!        do l=1,nxy                      ! FAST REACTOR
!          la=ltola(l)                   ! FAST REACTOR
!          iat=iassytyp(la)              ! FAST REACTOR
!          if(iffuela(iat)) then        ! FAST REACTOR
!            ichan=ichan+1               ! FAST REACTOR
!          endif                        ! FAST REACTOR
!        enddo                          ! FAST REACTOR
!        nchan=ichan                     ! FAST REACTOR
!      endif                            ! FAST REACTOR
!! added end      

! 2014_12_29 . scb
! count the number of channel
      if(rect) then
        nchan=0
        do ja=1,nya
          nthy1=nthy(ja)
          nthx1=0
          do ia=nxsfa(ja),nxefa(ja)
            nthx1=nthx1+nthx(ia)
          enddo
          nchan=nchan+nthx1*nthy1
        enddo          
      else
        ichan=0                         ! FAST REACTOR
        do l=1,nxy                      ! FAST REACTOR
          la=ltola(l)                   ! FAST REACTOR
          iat=iassytyp(la)              ! FAST REACTOR
          if(iffuela(iat)) then        ! FAST REACTOR
            ichan=ichan+1               ! FAST REACTOR
          endif                        ! FAST REACTOR
        enddo                          ! FAST REACTOR
        nchan=ichan                     ! FAST REACTOR
      endif                            ! FAST REACTOR
! added end      
	
! t/h allocation
      call allocth
                           
! channel assignment
      if(rect) then  ! 2014_12_29 . scb added if statement
        meshx(0)=0
        meshy(0)=0
        do ia=1,nxa
          meshx(ia)=meshx(ia-1)+nmeshx(ia)
        enddo
        do ja=1,nya
          meshy(ja)=meshy(ja-1)+nmeshy(ja)
        enddo

        ichan=0;lc=0;npthy=0
        do ja=1,nya
          npthy=nmeshy(ja)/nthy(ja)
          do ji=1,nthy(ja)
            jb=meshy(ja-1)+(ji-1)*npthy
            do ia=nxsfa(ja),nxefa(ja)
              npthx=nmeshx(ia)/nthx(ia)
              do ii=1,nthx(ia)
                ichan=ichan+1
                ib=meshx(ia-1)+(ii-1)*npthx
                j=jb
                do jn=1,npthy
                  j=j+1
                  i=ib
                  do in=1,npthx
                    i=i+1
                    l=nodel(i,j)
  ! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
  !                  ltochan(l)=ichan
                    if(.not.iflfr)  ltochan(l)=ichan
  ! added end                  
                    lc=lc+1
                    !lchantol(ic)=l   ! 2015_08_03 . scb commented. 
                    lchantol(ichan)=l   ! 2015_08_03 . 
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo      
      endif
! added end      
      
! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
!      if(ichan .ne. nchan) then
      if(.not.iflfr .and. rect .and. ichan.ne.nchan) then
! added end      
        call terminate("ERROR DURING ASSIGNING CHANNEL");
      endif
      
      
! channel to node correspondence for the radial reflector 
! added in ARTOS ver. 0.2 ( for LFR, CR ). 2012_07_03 by SCB
!      nchanp1=nchan+1
!      if(fdbk) then
!        do l=1,nxy
!          la=ltola(l)
!          iat=iassytyp(la)
!          if(.not.iffuela(iat)) ltochan(l)=nchanp1
!        enddo
!      endif
      nchanp1=nchan+1
      nchanp2=nchan+1
      ichan=0
      if(fdbk) then
        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          if(.not.iffuela(iat)) then
            if(irodtyp(la).ne.0) then
              ltochan(l)=nchanp2 ! control rod
            else
              ltochan(l)=nchanp1 ! reflector
            endif 
!          else                          ! 2012_08_22 . scb
          elseif(.not.rect) then        ! 2012_08_22 . scb
            ichan=ichan+1     ! FAST REACTOR
            ltochan(l)=ichan  ! FAST REACTOR
          endif
        enddo
      endif
! added end      

! in case of no feedback or fixed t/h, set t/h vectors with constants
      if(.not.fdbk .or. frozendm.ne.0 .or. frozentf.ne.0) then
         dmfix=dmref
         tmfix=tmref
         tffix=tfref
         if(frozendm.ne.0) dmfix=frozendm
         if(frozentf.ne.0) tffix=frozentf
         do kth=1,nzth
            do lth=1,nchanp1
               dcool(kth,lth)=dmfix
               tcool(kth,lth)=tmfix
               tdopl(kth,lth)=tffix
            enddo
         enddo
      endif
      !if(.not.fdbk) return   ! 2014_06_16 . scb
            
      
! initialize constants
      pfa=0.01*pfa0         !cm to m
      rs=rs*0.001           !mm to m
      rw=rw*0.001           !mm to m
      tw=tw*0.001           !mm to m
      rgt=rgt*0.001         !mm to m
      rs2=rs*rs
      rw2=rw*rw
      rg=rw-tw
      rg2=rg*rg

      powfa=powfa0*1e6       !Mw to w

      if(rect) then   ! 2014_12_29 . scb added if statement
        xx=1
        do ia=1,nxa
          if(xx.lt.nthx(ia)) xx=nthx(ia)
        enddo
        yy=1
        do ja=1,nya
          if(yy.lt.nthy(ja)) yy=nthy(ja)
        enddo
        nperfa=xx*yy
      
        chanvf=1./nperfa         !channel volume fraction per assembly
        acf=(pfa**2-PI*(npint*rw2+ngt*rgt**2))*chanvf    !coolant flow area
      else
! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
        chanvf = 1.d0
        acf=(pfa**2*0.86602540378444-PI*(npint*rw2+ngt*rgt**2))*chanvf  !txk: hex assy area correction
      endif
! added end

      afp=npint*PI*rs2*chanvf                          !total fuel pellet area
      xi=2*pi*(npint*rw+ngt*rgt)*chanvf                !wetted perimeter
      zeta=npint*2*pi*rw*chanvf                        !heated perimeter
      cfrperfa=cfrperfa*chanvf                        ! coolant flow rate
      zetap=zeta/acf                                  !heated perimeter density
      deq=4*acf/xi                                    !equivalent diameter
      hac=coreheight*0.01                        !active core height in meter
      powlin=powfa/hac*chanvf                 !linear power density
      fracdf=1-fracdc                         !fuel heat deposit fraction            
      
! calculate relative channel volume that equals to radial area of channel.
      k=kfbeg
      do l=1,nxy
        la=ltola(l)
        iat=iassytyp(la)
        if(.not.iffuela(iat)) cycle
        lchan=ltochan(l)
        chanvol(lchan)=chanvol(lchan)+volnode(l,k)
      enddo
      do lchan=1,nchan
        chanvol(lchan)=chanvol(lchan)*rhmesh(ZDIR,1,k)
      enddo

! axial nodalization
      nzthp1=nzth+1                           !number of junctions
      do kth=1,nzth
         hzth(kth)=0
         do k=junb(kth-1)+1,junb(kth)
            ka=ktoka(k)
            hzth(kth)=hzth(kth)+hmesh(ZDIR,1,k)
            ktokth(k)=kth
         enddo
         hzth(kth)=hzth(kth)*0.01             !cm to m
      enddo
      
! radial nodalization (in fuel pin)
      delr=rs/float(nr)
      delr2=delr*delr
      tw2=tw*tw
      delrw=0.5*tw
      delrw2=delrw*delrw
      tworm=tw/(rg+0.5*tw)
      do i=1,nrp1
         r(i)=delr*(i-1)
      enddo
      r(nrp2)=rg
      r(nrp3)=rg+delrw
      r(nrp4)=rw
      kgap=hgap*delr
      kgap2=hgap*tw*rs/rg
      kgap4=hgap*tw*(4-tw/rg)*rs/rg
      
! assign inlet condition to all nodes
      tdopin=sqrt(tin+CKELVIN) !inlet doppler temperature
      din=fdens(tin) !inlet densit
      do l=1,nchan
        do k=1,nzth
          do ir=1,nrp1
            tfuel(ir,k,l)=tin
            tdopl(k,l)=tdopin
          enddo
          do ir=nrp2,nrp4
            tfuel(ir,k,l)=tin
          enddo
          tcool(k,l)=tin
          dcool(k,l)=din
        enddo
      enddo

! assign reflector t/h condition
      do kth=1,nzth
        dcool(kth,nchanp1)=din
        tcool(kth,nchanp1)=tin
        tdopl(kth,nchanp1)=tdopin
      enddo      
      
! initialize vloume and junction variables
      hin=fenthal(tin)
      rhoin=fdens(tin)
      rhouin=cfrperfa/acf
      rhohuin=rhouin*hin
      uin=rhouin/rhoin
      do l=1,nchan
         k=0
         rhou(k,l)=rhouin
         rhohu(k,l)=rhohuin
         u(k,l)=uin
         ud(k,l)=uin
         do k=1,nzth
            hcool(k,l)=hin
            rhou(k,l)=rhouin
            rhohu(k,l)=rhohuin
            u(k,l)=uin
            ud(k,l)=uin
         enddo
      enddo
      
! initialze power with cosine shape
! added in ARTOS ver. 0.2 . 2012_07_03 by SCB
      kfsth=ktokth(kfbeg)
      kfeth=ktokth(kfend)
! added end      
      if(.not.isflatpower) then
         rsaving=10.
         buckl=pi/(coreheight+2*rsaving)
         p0=hac/(cos(buckl*rsaving)-cos(buckl*(coreheight+rsaving)))
         cosi=cos(buckl*rsaving)
! added in ARTOS ver. 0.2 . 2012_07_03 by SCB
!         kfsth=ktokth(kfs)
!         kfeth=ktokth(kfe)
! added end      
         do kth=1,kfsth
            relp(kth,1)=0
         enddo
         do kth=kfeth,nzth
            relp(kth,1)=0
         enddo
         z=rsaving
         do k=kfsth,kfeth
            z=z+hzth(k)*100
            cosip1=cos(buckl*z)
            relpk=p0*(cosi-cosip1)/hzth(k)
            if(relpk.lt.0) relpk=0
            cosi=cosip1
            ptot=ptot+relpk*hzth(k)
            relp(k,1)=relpk
         enddo
         do k=1,nzth
            do lth=2,nchan
               relp(k,lth)=relp(k,1)
            enddo
         enddo
      else
        call updrelpow(FALSE)
      endif
      
!! initialize mask index
!      if(ifmslb) then               ! MSLB
!        do lth=0,nchan+1            ! MSLB
!          do kth=1,nzth             ! MSLB   
!            itabindx(1,kth,lth)=1   ! MSLB  
!            itabindx(2,kth,lth)=1   ! MSLB 
!          enddo                     ! MSLB 
!        enddo                       ! MSLB 
!      endif                         ! MSLB 
	
      return
    end subroutine
      