! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_04 by SCB	
    subroutine updxsecf(fdbk)
! calculated nodal cross sections at a given t/h condition
      use param
      use timer
      use allocs
      use trinx_cntl
      use tranxsec, only : xschifd
      use sfam_cntl, only : nodalupd, cmfdupd

      include 'global.h'
      include 'times.h'
      include 'xsec.h'
      include 'geom.h'
      include 'files.h'
      include 'srchppm.h'
      include 'ffdm.h'
      include 'thfdbk.inc'
      include 'thfuel.inc'
      include 'thgeom.inc'
      include 'thop.inc'
!      include 'mslb.inc' ! MSLB
      include 'thlink.inc'
      include 'thexpan.inc' ! THERMAL EXPANSION
      include 'trancntl.inc'

      logical, intent(in) :: fdbk   
      logical :: ifskip,  first_fuel, first_cool, iftran   
      character*6 :: compname
      logical,save :: first=TRUE
      real,pointer :: tfueld(:,:)

      iftran=isteptr.ne.0
      if(thexpan) then
        write(mesg,'(a20,2f20.3)')     'Avg. Temp : ', tfuelavg, tcoolavg
        call message(false,false,mesg)

        call updexpan(0)
        call updgeom
      endif

      if(first) then
        call dmalloc(tfueld,nzth,nchan)
        do k=kfbeg,kfend
          kth=ktokth(k)
          ka=ktoka(k)
          do l=1,nxy
            la=ltola(l)
            iat=iassytyp(la)
            icomp=icompnumz(ka,iat)
            if(iffuelc(icomp)) then
              lth=ltochan(l)
              tfueld(kth,lth)=tfuel(nrp5,kth,lth)
            endif
          enddo
        enddo	
      endif

      ! update xsec
      iskip=0
      kthd=0
      first_cool=TRUE
      do k=1,nz
        ifskip=FALSE
        ka=ktoka(k)

        kth=ktokth(k)
        if(kth.eq.kthd) ifskip=TRUE
        kthd=kth

        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          icomp=icompnumz(ka,iat)
          compname=ictocn(icomp)
          lchan=ltochan(l) 
          itype=icomptyp(icomp)
          do m=1,ng
            if(nodalupd) xstfn(m,l,k)=xstf(m,l,k)
            if(cmfdupd) xstfc(m,l,k)=xstf(m,l,k)
          enddo

          ! fuel assembly
          if(lchan.le.nchan) then
            if(ilink.lt.0) then
              tfueln=tfuel(nrp5,kth,lchan)+CKELVIN
              tcooln=tcool(kth,lchan)+CKELVIN
            else
              tfueln=tdopl(kth,lchan)*tdopl(kth,lchan)
              tcooln=tcool(kth,lchan)
            endif
            if(itype.eq.COREM) then
              if(.not.ifskip) then			
                dt=abs(tfuel(nrp5,kth,lchan)-tfueld(kth,lchan))
                if(dt.lt.epstemp) then
                  iskip=iskip+1
                  if(.not.first) cycle
                else
                  tfueld(kth,lchan)=tfuel(nrp5,kth,lchan)
                endif

                call updexpan(itype)											
                call trinx(FALSE,icomp,compname,itype,rdum) ! TRINX
                call nodexs(l,k,icomp)
              else
                km1=k-1
                do m=1,ng
                  xschif(m,l,k)=xschif(m,l,km1) 
                  xsdf(m,l,k)=xsdf(m,l,km1)
                  xsdf2(m,l,k)=9.d0/7.d0*xsdf(m,l,k)   ! 2013_07_15 . scb
                  xsaf(m,l,k)=xsaf(m,l,km1)
                  xstf(m,l,k)=xstf(m,l,km1)
                  xsnff(m,l,k)=xsnff(m,l,km1)
                  xsff(m,l,k)=xsff(m,l,km1)    ! 2014_07_31 . scb
                  xsnuf(m,l,k)=xsnff(m,l,k)/(xsff(m,l,k)+1.e-20)   ! 2014_10_24 . scb                  
                  !xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb                  
                  xskpf(m,l,k)=xskpf(m,l,km1)
                  xssf(sigsms(m,icomp):sigsme(m,icomp),m,l,k) = xssf(sigsms(m,icomp):sigsme(m,icomp),m,l,km1)

                  xbeta(m,l,k,XDIR)=xsdf(m,l,k)*rhmesh(XDIR,l,k)
                  xbeta(m,l,k,YDIR)=xsdf(m,l,k)*rhmesh(YDIR,l,k)
                  xbeta(m,l,k,ZDIR)=xsdf(m,l,k)*rhmesh(ZDIR,l,k)
! 2013_07_15 . scb            
	                xbeta2(m,l,k,XDIR)=xsdf2(m,l,k)*rhmesh(XDIR,l,k)
                  xbeta2(m,l,k,YDIR)=xsdf2(m,l,k)*rhmesh(YDIR,l,k)
	                xbeta2(m,l,k,ZDIR)=xsdf2(m,l,k)*rhmesh(ZDIR,l,k)
! added end            
                enddo ! m
              endif ! .not. ifskip
            else
              if(itype.eq.BACK) itype=AREFL
              call updexpan(itype)											
              call trinx(FALSE,icomp,compname,itype,rdum) ! TRINX
              call nodexs(l,k,icomp)
            endif
          ! control assembly
          elseif(lchan.eq.nchan+2) then
            if(itype.eq.CONTM) then
              irodtyp1=abs(irodtyp(la)) 
              rodfrac1=rodfrac(k,irodtyp1)
              call updexpan(itype)											
              call trinx(FALSE,icomp,compname,itype,rodfrac1) ! TRINX
              call nodexs(l,k,icomp)
            else
              if(itype.eq.BACK) itype=RREFL
              call updexpan(itype)											
              call trinx(FALSE,icomp,compname,itype,rdum) ! TRINX
              call nodexs(l,k,icomp)
            endif
          ! radial reflector
          else 
            if(first_cool) then
              call updexpan(RREFL)	
              call trinx(FALSE,icomp,compname,RREFL,rdum) ! TRINX
              first_cool=FALSE
            endif				
            if(itype.ne.CONTM) call nodexs(l,k,icomp)
          endif ! iffuelc 
        enddo ! l
      enddo ! k

      write(mesg,'(i10)') iskip
      call message(true,true,mesg)

      first=FALSE
      
      return
    end subroutine
!=======================================================================================!
    subroutine nodexs(l,k,icomp)
      use param
      use trinx_cntl
      use tranxsec, only : xschifd, xsbetak => betak, xslmbdk => lmbdk

      include 'global.h'
      include 'xsec.h'
      include 'geom.h'

      integer,intent(in) :: l,k,icomp
      logical :: iftran

      iftran=isteptr.ne.0
      do m=1,ng
        xschif(m,l,k)=sigchi(m,icomp)    
        xsdf(m,l,k)=sigd(m,icomp)
        xsdf2(m,l,k)=9.d0/7.d0*sigd(m,icomp)   ! 2013_07_15 . scb
        xsaf(m,l,k)=siga(m,icomp)
        xstf(m,l,k)=sigt(m,icomp)
        xsnff(m,l,k)=signf(m,icomp)
        xsff(m,l,k)=sigf(m,icomp)   ! 2014_07_31 . scb
        xsnuf(m,l,k)=xsnff(m,l,k)/(xsff(m,l,k)+1.e-20)   ! 2014_10_24 . scb        
        !xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb        
        if(iftran) then
          xsnff(m,l,k)=reigv0*xsnff(m,l,k) 
          xsnuf(m,l,k)=xsnff(m,l,k)/(xsff(m,l,k)+1.e-20)   ! 2014_10_24 . scb
          !xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb          
        endif
        xskpf(m,l,k)=sigkf(m,icomp)
        xssf(sigsms(m,icomp):sigsme(m,icomp),m,l,k) = sigsm(sigsms(m,icomp):sigsme(m,icomp),m,icomp)

        xbeta(m,l,k,XDIR)=xsdf(m,l,k)*rhmesh(XDIR,l,k)
        xbeta(m,l,k,YDIR)=xsdf(m,l,k)*rhmesh(YDIR,l,k)
        xbeta(m,l,k,ZDIR)=xsdf(m,l,k)*rhmesh(ZDIR,l,k)
! 2013_07_15 . scb            
	      xbeta2(m,l,k,XDIR)=xsdf2(m,l,k)*rhmesh(XDIR,l,k)
        xbeta2(m,l,k,YDIR)=xsdf2(m,l,k)*rhmesh(YDIR,l,k)
	      xbeta2(m,l,k,ZDIR)=xsdf2(m,l,k)*rhmesh(ZDIR,l,k)
! added end            
      enddo ! m

      return
    end subroutine
