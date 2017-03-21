! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
    subroutine updxsec_h(iffdbk)
! calculated nodal cross sections at a given t/h condition
      use param
      use timer
!
      include 'global.h'
      include 'times.h'
      include 'xsec.h'
      include 'geom.h'
      include 'files.h'
      include 'srchppm.h'
      include 'thfdbk.inc'
      include 'thfuel.inc'
      include 'thcntl.inc'
      include 'thgeom.inc'
      include 'thop.inc'
!      include 'mslb.inc' ! MSLB
      include 'trancntl.inc'
      include 'xesm.h'  ! 2013_09_27 . scb      
!
      logical, intent(in) :: iffdbk
!
      logical, save :: first=TRUE
!      
      real :: del(NUMOFDXS-1)
      equivalence(del(1), delppm)
      equivalence(del(2), deltm)
      equivalence(del(3), deldm)
      equivalence(del(4), deltf)

      call timeron() 
      
      sumdkp=0
! update xsec
      do k=1,nz
        ka=ktoka(k)
        if(first) then
          do l=1,nxy
            iat=iassytyp(l)
            icomp=icompnumz(ka,iat)
            xsmax(l,k)=sigsmax(icomp)   ! 2013_07_19 . scb
            do m=1,ng         
              xschif(m,l,k)=sigchi(m,icomp)
              xsadf(:,m,l,k)=sigadf(:,m,icomp)
              xscdf(:,m,l,k)=sigmdf(:,m,icomp)              
              
              xssfs(m,l,k)=sigsms(m,icomp)
              xssfe(m,l,k)=sigsme(m,icomp)
            enddo ! m
            
! 2013_09_27 . scb
            if(ixesmopt.ne.0) then
              gami(l,k)=gammafp(1,icomp)
              gamxe(l,k)=gammafp(2,icomp)
              gampm(l,k)=gammafp(3,icomp)
            endif              
! added end                          
          enddo ! la
        endif ! first

        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          icomp=icompnumz(ka,iat)
          irodtyp1=abs(irodtyp(la)) ! STUCK ROD

          do m=1,ng
            if(usesigtr) then
              xstrf1=sigtr(m,icomp)
              xstrf(m,l,k)=xstrf1
            else
              xsdf(m,l,k)=sigd(m,icomp)
              xsdf2(m,l,k)=9.d0/7.d0*sigd(m,icomp)   ! 2013_07_15 . scb
            endif
            
            if(usesiga) then
              xsaf(m,l,k)=siga(m,icomp)
            else
              xstf(m,l,k)=sigt(m,icomp)   ! 2014_03_20 . scb, ksw
            endif
            
            xsnff(m,l,k)=signf(m,icomp)
            xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb            
            xskpf(m,l,k)=sigkf(m,icomp)
            xssf(sigsms(m,icomp):sigsme(m,icomp),m,l,k) &
            = sigsm(sigsms(m,icomp):sigsme(m,icomp),m,icomp)

! control rod            
            if(irodtyp1 .ne. 0) then
              if(rodfrac(k,irodtyp1) .ne. 0) then
                if(usesigtr) then
                  xstrf1=xstrf1+rodfrac(k,irodtyp1)*dsigd_tr(DROD,m,icomp)
                  xstrf(m,l,k)=xstrf1
                  xsdf(m,l,k)=1./(3.*xstrf1)                        
                else
                  xsdf(m,l,k)=xsdf(m,l,k)+rodfrac(k,irodtyp1)*dsigd_tr(DROD,m,icomp)
                endif
                xsdf2(m,l,k)=9.d0/7.d0*xsdf(m,l,k)   ! 2013_07_15 . scb
                
                if(usesiga) then
                  xsaf(m,l,k)=xsaf(m,l,k)+rodfrac(k,irodtyp1)*dsigt_a(DROD,m,icomp)
                else
                  xstf(m,l,k)=xstf(m,l,k)+rodfrac(k,irodtyp1)*dsigt_a(DROD,m,icomp)
                endif
                
                xsnff(m,l,k)=xsnff(m,l,k)+rodfrac(k,irodtyp1)*dsignf(DROD,m,icomp)
                xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb                
                xskpf(m,l,k)=xskpf(m,l,k)+rodfrac(k,irodtyp1)*dsigkf(DROD,m,icomp)
                do ms=sigsms(m,icomp),sigsme(m,icomp)
                  xssf(ms,m,l,k)=xssf(ms,m,l,k)+dsigsm(ms,DROD,m,icomp)*rodfrac(k,irodtyp1)
                enddo
              endif !rodfrac(k,irodtyp1) .ne. 0
            endif ! irodtyp1 .ne. 0 

! density of boron, temperature of moderator
! density of moderator, temperature of fuel
            delppm=ppm-basecond(DPPM)
            ice=DPPM
            if(iffdbk) then
              kth=ktokth(k)
              lchan=ltochan(l)
! 2014_04_11 . scb              
              !if(iffuela(iat)) then
              !  deltm=tcool(kth,lchan)-basecond(DTM)
              !  deldm=dcool(kth,lchan)-basecond(DDM)
              !  deltf=tdopl(kth,lchan)-basecond(DTF)
              !else
              !  deltm=tin-basecond(DTM)
              !  deldm=din-basecond(DDM)
              !  deltf=tdopin-basecond(DTF)
              !endif
              deltm=tcool(kth,lchan)-basecond(DTM)
              deldm=dcool(kth,lchan)-basecond(DDM)
              deltf=tdopl(kth,lchan)-basecond(DTF)   
! test ! if there is no problem , below lines should be removed..              
#ifndef DLL
              ierr=0
              if(.not.iffuela(iat)) then
                if(lchan.ne.nchan+1) ierr=1
                if(tin.ne.tcool(kth,lchan)) ierr=1
                if(din.ne.dcool(kth,lchan)) ierr=1
                if(tdopin.ne.tdopl(kth,lchan)) ierr=1
              endif
              if(ierr.eq.1) stop 'there is error in inlet T/H condition'
#endif
! added end              
              ice=DTF
            endif

            do icond=DPPM,ice ! NUMOFDXS-1
              if(usesigtr) then
                xstrf1=xstrf1+dsigd_tr(icond,m,icomp)*del(icond)
                xstrf(m,l,k)=xstrf1
                xsdf(m,l,k)=1./(3.*xstrf1)
              else
                xsdf(m,l,k)=xsdf(m,l,k)+dsigd_tr(icond,m,icomp)*del(icond)
              endif
              xsdf2(m,l,k)=9.d0/7.d0*xsdf(m,l,k)   ! 2013_07_15 . scb

              if(usesiga) then
                xsaf(m,l,k)=xsaf(m,l,k)+dsigt_a(icond,m,icomp)*del(icond)
              else
                xstf(m,l,k)=xstf(m,l,k)+dsigt_a(icond,m,icomp)*del(icond)
              endif

              xsnff(m,l,k)=xsnff(m,l,k)+dsignf(icond,m,icomp)*del(icond)
              xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb              
              xskpf(m,l,k)=xskpf(m,l,k)+dsigkf(icond,m,icomp)*del(icond)
              sumdkp=sumdkp+del(icond)
              do ms=sigsms(m,icomp),sigsme(m,icomp)
               xssf(ms,m,l,k)=xssf(ms,m,l,k)+dsigsm(ms,icond,m,icomp)*del(icond)
              enddo
            enddo ! icond
          
	          xbeta(m,l,k,XDIR)=xsdf(m,l,k)*rhmesh(XDIR,l,k)
            xbeta(m,l,k,YDIR)=xsdf(m,l,k)*rhmesh(YDIR,l,k)
	          xbeta(m,l,k,ZDIR)=xsdf(m,l,k)*rhmesh(ZDIR,l,k)
            
! 2013_07_15 . scb            
	          xbeta2(m,l,k,XDIR)=xsdf2(m,l,k)*rhmesh(XDIR,l,k)
            xbeta2(m,l,k,YDIR)=xsdf2(m,l,k)*rhmesh(YDIR,l,k)
	          xbeta2(m,l,k,ZDIR)=xsdf2(m,l,k)*rhmesh(ZDIR,l,k)
! added end                        
	        enddo !m
        enddo !l
      enddo !k
      
      r4o3=4./3.
      r5o3=5./3.
      do k=1,nz
        do l=1,nxy
          do m=1,ng
            summ=0
            do md=1,ng
              if(m .ge. xssfs(md,l,k) .and. m.le.xssfe(md,l,k)) then
                summ=summ+xssf(m,md,l,k)
              endif
            enddo
            if(usesiga) then
              xstf(m,l,k)=xsaf(m,l,k)+summ
            else
              xsaf(m,l,k)=xstf(m,l,k)-summ
            endif
          enddo
        enddo
      enddo
      
      call timeroff(txsec)
      
      first=FALSE
      
      print '(a,6f16.4)', &
      'SUM OF XSEC', sum(xsdf),sum(xsaf),sum(xstf),sum(xsnff),sum(xskpf)

      return
    end subroutine
