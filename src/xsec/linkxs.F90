! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_04 by SCB
    subroutine linkxs(icomp,compname,comp) 

      use param
      use BasicParam
      use MaterialMod

      include 'global.h'
      include 'files.h'
      include 'cards.h'      
      include 'xsec.h'
      !include 'ntbytes.h'

      integer, intent(in) :: icomp
      character*6, intent(in) :: compname
      type(Composition), intent(in) :: comp

      do m=1,ng
        sigtr(m,icomp) = comp%xs_tr(m)
        siga(m,icomp)  = comp%xs_a(m)
        sigf(m,icomp)  = comp%xs_f(m)    ! 2014_10_24 . scb
        signf(m,icomp) = comp%xs_nf(m)
        sigkf(m,icomp) = comp%xs_kf(m)
        sigchi(m,icomp)= comp%chip(m)
        sigchid(m,icomp)=comp%chid(m)

        sigd(m,icomp) = 1._8/(3._8*sigtr(m,icomp))
        if(signf(m,icomp).ne.0) iffuelc(icomp)=TRUE
      enddo

      do mp=1,nprec
        betak(mp,icomp)=comp%beta(mp)
        lmbdk(mp,icomp)=comp%dct(mp)
      enddo

      !if(iffuelc(icomp)) then
      sum=0.d0
      sumd=0.d0
      do m=1,ng
        sum=sum+sigchi(m,icomp)
        sumd=sumd+sigchid(m,icomp)
      enddo
        
      if(sum.lt.1.e-30) sum=1.e-30   ! 2013_10_10 . scb 
      if(sumd.lt.1.e-30) sumd=1.e-30   ! 2013_10_10 . scb  
        
      fnorm=1.d0/sum
      fnormd=1.d0/sumd
      do m=1,ng
        sigchi(m,icomp)=sigchi(m,icomp)*fnorm
        sigchid(m,icomp)=sigchid(m,icomp)*fnormd
      enddo
      !endif

! scattering matrix
      smsq=0
      do m = 1, ng
        do md = 1, ng
          smsq(m,md) = comp%xs_s(m,md) 
        enddo

        sum=0
        do md=1,ng
          sum=sum+smsq(m,md)
        enddo

        sigs(m,icomp)=sum
        sigt(m,icomp)=siga(m,icomp)+sum
        if(siga(m,icomp).lt.0.) then
          siga(m,icomp)=0
          sigt(m,icomp)=sum
        endif
      enddo
!
      do m=1,ng
        ib=m
        ie=m
        do ms=1,ng
          if(smsq(ms,m).ne.0) then
            ib=min(ib,ms)
            exit
          endif
        enddo
        do ms=ng,1,-1
          if(smsq(ms,m).ne.0) then
            ie=max(ie,ms)
            exit
          endif
        enddo

        ! upscattering bound
        if(ie.gt.m) maxupscatgr=min(maxupscatgr,m)
        !
        sigsms(m,icomp)=ib
        sigsme(m,icomp)=ie
        sigsm(ib:ie,m,icomp)=smsq(ib:ie,m)

        ! remove self scattering from the scattering matrix
        sigt(m,icomp)=sigt(m,icomp)-sigsm(m,m,icomp)
        sigs(m,icomp)=sigs(m,icomp)-sigsm(m,m,icomp)
        sigsm(m,m,icomp)=0        
      enddo 

      sigsmax(icomp)=maxupscatgr   ! 2013_07_19 . scb      
      
      return          
    end subroutine 

