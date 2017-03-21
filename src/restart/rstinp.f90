    subroutine rstinp
!
      use param
      use sfam,   only : eigv, phif, psi, psid
      use sfam,   only : eigv0  ! 2014_09_04 . scb
      use tran,   only : deltm, deltmd, psit, psitd, prec, betap
      use cmfdmg, only : dhatrf, dhatzf, ccrf, cczf
      use cmfd2g, only : am, af, ccr, ccz
      use bdf,    only : phifbdf, phifbdf2, deltmarray  ! 2013_05_20 . scb
      use bdf,    only : flagbdf    ! 2015_07_07 . scb
      use senm2n, only : phicff  ! 2013_05_20 . scb
      use trinx_cntl, only : isteptr  ! 2013_05_20 . scb
      use allocs    ! 2013_05_22 . scb
            
      use tpen,       only : fcnto, fcntzo, hflx, pflx, aflx, xmom, ymom, hflxf, zmom, fhflx
      use tpen_sp3,   only : cnto2, cntzo2, hflx2, pflx2, aflx2, xmom2, ymom2, hflxf2, zmom1, zmom2, fhflx2
      use cmfdhex2g,  only : dhat, dhatz, dsum, betaphis, betaphisz, betaphis2, betaphisz2
      use geomhex,    only : wtdhat
      
      use sp3senm,    only : phi,phishp,lkg2   ! 2014_12_23 . scb      
            
      include 'global.h'
      include 'files.h'
      include 'trancntl.inc'
      include 'xsec.h'
      include 'thgeom.inc'
      include 'cntl.h'
      include 'srchppm.h'
      include 'geom.h'
      include 'thcntl.inc'
      include 'thcool.inc'
      include 'thfdbk.inc'
      include 'thfuel.inc'
      include 'pow.h'
      include 'thop.inc'
      include 'thexpan.inc'
      include 'geomh.h'
      include 'thlink.inc' ! 2014_04_10 . scb
      include 'xesm.h'   ! 2014_09_04 . scb
!      
      character*65 amesg65
! 2013_05_22 . scb
      real,pointer :: facscb(:,:,:)
      real,pointer :: facscbp(:,:,:)  ! 2014_03_02 . scb
      integer,pointer :: nassyp(:)    ! 2014_03_02 . scb
      common / facscbsp3 / facscb
      common / facscbsp32 / facscbp, nassyp
      
      real,pointer,dimension(:,:) :: psrc,pcoef,psrctem,pflxt,xsdmwt(:,:,:)
      common / rsthexsp3 / psrc,pcoef,psrctem,pflxt,xsdmwt
      
      if(hex) then
        call dmalloc(facscb,ng,nxy,nz)
        call dmalloc(facscbp,ng,ncorn,nz) 
        call dmalloc(nassyp,ncorn)
        if(ifhexsp3) then
          call dmalloc(psrc,2,ncorn)
          call dmalloc(pcoef,2,ncorn)
          call dmalloc(psrctem,2,ncorn)
          call dmalloc(pflxt,2,ncorn)
          call dmalloc(xsdmwt,2,3,nassy)
        endif
      endif      
! added end
!      
      !if(hex) then
      !  call rstinphex
      !  return
      !endif
!--- Read time-independent restart variables.
! decay heat...
      read(irstin) omalphatot
! eigenvalue...
      read(irstin) eigv
      read(irstin) eigv0   ! 2014_09_04 . scb
! precursor...
      read(irstin) nprec
200   continue

!--- Read time-dependent restart variables.
      read(irstin,end=999,err=999) irstadv,rsttime,deltm,isteptr
! boron
      read(irstin) ppm
! rod step
      do ib=1,nrodtyp
        read(irstin) rodstep(ib)
      enddo
      
      if(fdbk) then
  ! coolant
        do l=1,nchan
          do k=1,nzth
            read(irstin) hcool(k,l), qeff(k,l)
          enddo
        enddo
        do l=1,nchan
          do k=0,nzth
            read(irstin) rhou(k,l), rhohu(k,l), u(k,l), ud(k,l)
          enddo
        enddo
  ! decay heat
        do k=1,nz
          do l=1,nxy
            do m=1,ndecgrp
              read(irstin) precconc(m,l,k)
            enddo
          enddo
        enddo
        do m=1,ndecgrp
          read(irstin) ezetdelt(m)
          read(irstin) omexpprod(m)
        enddo
  ! T/H
        do l=0,nchan+1
          do k=0,nzth
            read(irstin) tdopl(k,l)
            read(irstin) tcool(k,l)
            read(irstin) dcool(k,l)
          enddo
        enddo
        do k=1,nz
          do l=1,nxy
            read(irstin) ppml(l,k)
          enddo
        enddo

        do l=1,nchan
          do k=1,nzth
            do m=1,nrp5
              read(irstin) tfuel(m,k,l)
            enddo
            read(irstin) htcoef(k,l), qvol(k,l)
          enddo
        enddo      
        
! power.h variables
        do l=1,nchan
          do k=1,nzth
            read(irstin) relp(k,l)
          enddo
        enddo        
! T/H property          
        read(irstin) tfmax, tmmax, tfuelavg, tcoolavg, dcoolavg    ! 2014_09_05 . scb added dcoolavg
      endif      
! power.h variables
      read(irstin) plevel, plevel0, plevel00, pleveld  ! 2014_09_15 . scb added plevel00
! flux
      do k=1,nz
        do l=0,nxy
          do m=1,ng
            read(irstin) phif(m,l,k)
          enddo
        enddo
      enddo
!
      do k=1,nz
        do l=1,nxy
          read(irstin) psi(l,k),psid(l,k)
        enddo
      enddo
! trcntl.h variables
      read(irstin) deltmd, trip, scram, tripbeg, tscrmbeg
      do ib=1,nrodtyp
        read(irstin) pscrmbeg(ib)
      enddo
! xsec.h variables
      do k=1,nz
        do l=1,nxy
          do m=1,ng
            read(irstin) xsdf(m,l,k), xsaf(m,l,k), xstf(m,l,k), xsnff(m,l,k), xskpf(m,l,k), xsff(m,l,k)
            do ms=1,ng
              read(irstin) xssf(ms,m,l,k)
            enddo              
          enddo
        enddo
      enddo
! 2013_05_20 . scb      
      do icomp=1,ncomp
        do m=1,ng
          read(irstin) signf(m,icomp)
          do icond=DPPM,NUMOFDXS
            read(irstin) dsignf(icond,m,icomp)
          enddo 
        enddo
      enddo      
! added end      
! neutronics

      !do i=1,bdforder
      !do i=1,bdforder0   ! 2014_09_04 . scb
      if(flagbdf) then   ! 2015_07_07 . scb      
        do i=1,5   ! 2014_12_29 . scb
          do k=1,nz
            do l=1,nxy
              do m=1,ng
                read(irstin) phifbdf(m,l,k,i)
              enddo              
            enddo
          enddo
        enddo        
      endif
      
      do k=1,nz
        do l=1,nxy
          do m=1,nprec
            read(irstin) prec(m,l,k)
          enddo
          read(irstin) psit(l,k),psitd(l,k)
          read(irstin) betap(l,k)
        enddo
      enddo
      if(flagbdf) read(irstin) deltmarray  ! 2013_05_20 . scb
      ! 2015_07_07 . scb added flagbdf        
      
      if(ifsp3) then
        do k=1,nz
          do l=1,nxy
            do m=1,ng
              !do i=1,bdforder
              !do i=1,bdforder0   ! 2014_09_04 . scb
              do i=1,5   ! 2014_12_29 . scb
                read(irstin) phifbdf2(m,l,k,i)
              enddo              
            enddo
          enddo
        enddo
      endif
      
      
      if(rect) then
! dhat...
        do k=1,nz
          do l=1,nsurf
            do m=1,ng
              read(irstin) dhatrf(m,l,k)
            enddo
          enddo
        enddo
        
        do k=1,nzp1
          do l=1,nxy
            do m=1,ng
              read(irstin) dhatzf(m,l,k)
            enddo
          enddo
        enddo
      
        !read(irstin) phicff  ! 2013_05_20 . scb
        if(ifsp3) then   ! 2014_12_23 . scb
          read(irstin) phishp
          read(irstin) phi
          read(irstin) lkg2
        else
          read(irstin) phicff  ! 2013_05_20 . scb
        endif
! linear system
        do k=1,nz
          do l=1,nxy
            do m=1,ng2
              read(irstin) am(m,l,k), af(m,l,k), ccr(1:4,m,l,k), ccz(1:2,m,l,k)
            enddo
          enddo
        enddo
        
        if(ng.gt.2) then
          do k=1,nz
            do l=1,nxy
              do m=1,ng
                read(irstin) ccrf(1:4,m,l,k), cczf(1:2,m,l,k)
              enddo
            enddo
          enddo
        endif        
      else
        read(irstin) dhat
        read(irstin) dhatz
        read(irstin) wtdhat
        read(irstin) dsum
        read(irstin) af
        read(irstin) betaphis
        read(irstin) betaphisz    
        
        if(ng.ne.ng2) then
          read(irstin) fcnto
          read(irstin) fcntzo
        endif                 
        
        read(irstin) facscb
        read(irstin) facscbp
        read(irstin) nassyp
        if(ifhexsp3) then
          read(irstin) hflx2
          read(irstin) pflx2
          read(irstin) betaphis2
          read(irstin) betaphisz2  
! 2013_05_22 . scb    
          read(irstin) psrc
          read(irstin) pcoef
          read(irstin) psrctem
          read(irstin) pflxt
          read(irstin) xsdmwt    
          read(irstin) cnto2
          read(irstin) cntzo2  
! added end             
          read(irstin) aflx2
          read(irstin) xmom2
          read(irstin) ymom2
          read(irstin) hflxf2
          read(irstin) zmom1
          read(irstin) zmom2
          read(irstin) fhflx2
        else          
          read(irstin) hflx
          read(irstin) pflx
          if(ng.ne.ng2) then
            read(irstin) aflx
            read(irstin) xmom
            read(irstin) ymom
            read(irstin) hflxf
            read(irstin) zmom
            read(irstin) fhflx
          endif       
        endif        
      endif      
  
! 2014_09_04 . scb       
      read(irstin)  flxlevel
      !flxlevel = 1.d0  ! 2014_09_14 . scb for dbg
      !print *, flxlevel  ! 2014_09_05 . scb for dbg
      read(irstin)  xeave, smave
      do k=1,nz
        do l=1,nxy
          read(irstin)  rni(l,k),rnxe(l,k),rnpm(l,k),rnsm(l,k)
          read(irstin)  rnip(l,k),rnxep(l,k),rnpmp(l,k),rnsmp(l,k)
        enddo
      enddo

!! 2014_09_16 . scb for dbg
!        print *, rnip(18,2),rnxep(18,2),rnpmp(18,2),rnsmp(18,2)
!        pause
! !added end      
! added end
            
      if (irstadv.lt.irstbeg) goto 200
      if (irstadv.eq.irstbeg) then
        call message(TRUE,TRUE,'Restart File Successfully Processed.')
        write(amesg65,'(11x,a,i5,a,f15.3)') 'Advancement Unit: ',irstadv,', Time: ',rsttime
        call message(FALSE,TRUE,amesg65)
        time=rsttime
        goto 500
      endif
      
 999  continue

      print *
      print *, " *** RESTART TIME REQUESTED DOES NOT EXIST. "
      print *, " ***   Requested block: ",irstbeg
      print *, " ***   Last block read: ",irstadv
      print *
!
      write(io8,*)
      write(io8,*)" *** RESTART TIME REQUESTED DOES NOT EXIST. "
      write(io8,*)" ***   Requested block: ",irstbeg
      write(io8,*)" ***   Last block read: ",irstadv
      write(io8,*)
      
      stop  ! 2013_05_16 . scb
      
500   continue
      
      return
    end

