! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_04 by SCB
    subroutine initxsecf
! calculated nodal cross sections at a given t/h condition for LFR
      use param
      use timer
!
      include 'global.h'
      include 'times.h'
      include 'xsec.h'
      include 'geom.h'
      include 'files.h'
      include 'srchppm.h'
      include 'ffdm.h'
      include 'thfdbk.inc'
      include 'thfuel.inc'
      include 'thcntl.inc'
      include 'thgeom.inc'
      include 'thop.inc'
!      include 'mslb.inc' ! MSLB
      include 'thlink.inc'
!
      logical, save :: first=TRUE

! initialize xsec
      if(.not.first) return
      
      if(first) then
        do k=1,nz
          ka=ktoka(k)
          do l=1,nxy
            la=ltola(l)
            iat=iassytyp(la)
            icomp=icompnumz(ka,iat)   
            xsmax(l,k)=sigsmax(icomp)   ! 2013_07_19 . scb
            do m=1,ng         
              xsadf(:,m,l,k)=1.
              xscdf(:,m,l,k)=1.                     
              xssfs(m,l,k)=sigsms(m,icomp)
              xssfe(m,l,k)=sigsme(m,icomp)

              xschif(m,l,k)=sigchi(m,icomp)  
              xstrf(m,l,k)=sigtr(m,icomp)  
              xsdf(m,l,k)=sigd(m,icomp)
              xsdf2(m,l,k)=9./7.*sigd(m,icomp)

              xsaf(m,l,k)=siga(m,icomp)
              xstf(m,l,k)=sigt(m,icomp)
              xsnff(m,l,k)=signf(m,icomp)
              xsff(m,l,k)=sigf(m,icomp)   ! 2014_07_31 . scb
              xsnuf(m,l,k)=xsnff(m,l,k)/(xsff(m,l,k)+1.e-20)   ! 2014_10_24 . scb
              !xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb
              xskpf(m,l,k)=sigkf(m,icomp)
              xss2nf(m,l,k)=sigsm(m,m,icomp)
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
          enddo ! l
        enddo !k
      endif
      first=FALSE
 
      return
    end subroutine
