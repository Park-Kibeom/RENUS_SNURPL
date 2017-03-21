    subroutine rstedit(ivaropt)
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
      
      use tpen,       only : fcnto, fcntzo, hflx, pflx, aflx, xmom, ymom, hflxf, zmom, fhflx
      use tpen_sp3,   only : cnto2, cntzo2, hflx2, pflx2, aflx2, xmom2, ymom2, hflxf2, zmom1, zmom2, fhflx2
      use cmfdhex2g,  only : dhat, dhatz, dsum, betaphis, betaphisz, betaphis2, betaphisz2
      use geomhex,    only : wtdhat, ifhexsp3
      
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
      include 'thlink.inc' ! 2014_04_10 . scb
      include 'xesm.h'   ! 2014_09_04 . scb
      
! Write restart information to file.
      character*80 amesg80
!
      integer ivaropt
!      
! 2013_05_22 . scb
      real,pointer :: facscb(:,:,:)
      real,pointer :: facscbp(:,:,:)  ! 2014_03_02 . scb
      integer,pointer :: nassyp(:)    ! 2014_03_02 . scb
      common / facscbsp3 / facscb
      common / facscbsp32 / facscbp, nassyp
      
      real,pointer,dimension(:,:) :: psrc,pcoef,psrctem,pflxt,xsdmwt(:,:,:)
      common / rsthexsp3 / psrc,pcoef,psrctem,pflxt,xsdmwt
! added end      
!
!--- Write time-independent restart variables.
      if (ivaropt.eq.0) then
! decay heat...
        write(irstout) omalphatot
! eigenvalue...
        write(irstout) eigv
        write(irstout) eigv0   ! 2014_09_04 . scb
! precursor...
        write(irstout) nprec
      else
!--- Write time-dependent restart variables.
        write(amesg80,810) irstadv,time
        call message(FALSE,TRUE,amesg80) 
810     format(11x,"Restart advancement: ",i4," written at time: ",f10.3," s")
! transient variables: irstadv, time
        write(irstout) irstadv, time, deltm, isteptr
! boron
        write(irstout) ppm
! rod step
        do ib=1,nrodtyp
          write(irstout) rodstep(ib)
        enddo
        if(fdbk) then
  ! coolant
          do l=1,nchan
            do k=1,nzth
              write(irstout) hcool(k,l), qeff(k,l)
            enddo
          enddo
          do l=1,nchan
            do k=0,nzth
              write(irstout) rhou(k,l), rhohu(k,l), u(k,l), ud(k,l)
            enddo
          enddo
  ! decay heat
          do k=1,nz
            do l=1,nxy
              do m=1,ndecgrp
                write(irstout) precconc(m,l,k)
              enddo
            enddo
          enddo
          do m=1,ndecgrp
            write(irstout) ezetdelt(m)
            write(irstout) omexpprod(m)
          enddo
  ! T/H
          do l=0,nchan+1
            do k=0,nzth
              write(irstout) tdopl(k,l)
              write(irstout) tcool(k,l)
              write(irstout) dcool(k,l)
            enddo
          enddo
          do k=1,nz
            do l=1,nxy
              write(irstout) ppml(l,k)
            enddo
          enddo

          do l=1,nchan
            do k=1,nzth
              do m=1,nrp5
                write(irstout) tfuel(m,k,l)
              enddo
              write(irstout) htcoef(k,l), qvol(k,l)
            enddo
          enddo
! power.h variables
          do l=1,nchan
            do k=1,nzth
              write(irstout) relp(k,l)
            enddo
          enddo        
! T/H property          
          write(irstout) tfmax, tmmax, tfuelavg, tcoolavg, dcoolavg    ! 2014_09_05 . scb added dcoolavg
        endif        
! power.h variables
        write(irstout) plevel, plevel0, plevel00, pleveld   ! 2014_09_15 . scb added plevel00
! flux
        do k=1,nz
          do l=0,nxy
            do m=1,ng
              write(irstout) phif(m,l,k)
            enddo
          enddo
        enddo
!
        do k=1,nz
          do l=1,nxy
            write(irstout) psi(l,k),psid(l,k)
          enddo
        enddo
! trcntl.h variables
        write(irstout) deltmd, trip, scram, tripbeg, tscrmbeg
        do ib=1,nrodtyp
          write(irstout) pscrmbeg(ib)
        enddo
! xsec.h variables
        do k=1,nz
          do l=1,nxy
            do m=1,ng
              write(irstout) xsdf(m,l,k), xsaf(m,l,k), xstf(m,l,k), xsnff(m,l,k), xskpf(m,l,k), xsff(m,l,k)
              do ms=1,ng
                write(irstout) xssf(ms,m,l,k)
              enddo              
            enddo
          enddo
        enddo
! 2013_05_20 . scb        
        do icomp=1,ncomp
          do m=1,ng
            write(irstout) signf(m,icomp)
            do icond=DPPM,NUMOFDXS
              write(irstout) dsignf(icond,m,icomp)
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
                  write(irstout) phifbdf(m,l,k,i)
                enddo              
              enddo
            enddo
          enddo        
        endif
        
        do k=1,nz
          do l=1,nxy
            do m=1,nprec
              write(irstout) prec(m,l,k)
            enddo
            write(irstout) psit(l,k),psitd(l,k)
            write(irstout) betap(l,k)
          enddo
        enddo
        if(flagbdf) write(irstout) deltmarray  ! 2013_05_20 . scb      
        ! 2015_07_07 . scb added flagbdf        
        
        !if(ifsp3) then
        if(ifsp3 .and. flagbdf) then   ! 2015_07_07 . scb
          do k=1,nz
            do l=1,nxy
              do m=1,ng
                !do i=1,bdforder
                !do i=1,bdforder0   ! 2014_09_04 . scb
                do i=1,5   ! 2014_12_29 . scb
                  write(irstout) phifbdf2(m,l,k,i)
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
                write(irstout) dhatrf(m,l,k)
              enddo
            enddo
          enddo
        
          do k=1,nzp1
            do l=1,nxy
              do m=1,ng
                write(irstout) dhatzf(m,l,k)
              enddo 
            enddo
          enddo
        
          !write(irstout) phicff  ! 2013_05_20 . scb
          if(ifsp3) then   ! 2014_12_23 . scb
            write(irstout) phishp
            write(irstout) phi
            write(irstout) lkg2
          else
            write(irstout) phicff  ! 2013_05_20 . scb
          endif
          
! linear system
          do k=1,nz
            do l=1,nxy
              do m=1,ng2
                write(irstout) am(m,l,k), af(m,l,k), ccr(1:4,m,l,k), ccz(1:2,m,l,k)
              enddo
            enddo
          enddo
        
          if(ng.gt.2) then
            do k=1,nz
              do l=1,nxy
                do m=1,ng
                  write(irstout) ccrf(1:4,m,l,k), cczf(1:2,m,l,k)
                enddo
              enddo
            enddo
          endif      
        else
          write(irstout) dhat
          write(irstout) dhatz
          write(irstout) wtdhat
          write(irstout) dsum
          write(irstout) af
          write(irstout) betaphis
          write(irstout) betaphisz  
        
          if(ng.ne.ng2) then
            write(irstout) fcnto
            write(irstout) fcntzo
          endif                 
        
          write(irstout) facscb  
          write(irstout) facscbp
          write(irstout) nassyp
          if(ifhexsp3) then
            write(irstout) hflx2
            write(irstout) pflx2
            write(irstout) betaphis2
            write(irstout) betaphisz2  
! 2013_05_22 . scb            
            write(irstout) psrc
            write(irstout) pcoef
            write(irstout) psrctem
            write(irstout) pflxt
            write(irstout) xsdmwt  
            write(irstout) cnto2
            write(irstout) cntzo2
! added end            
            write(irstout) aflx2
            write(irstout) xmom2
            write(irstout) ymom2
            write(irstout) hflxf2
            write(irstout) zmom1
            write(irstout) zmom2
            write(irstout) fhflx2   
          else          
            write(irstout) hflx
            write(irstout) pflx  
            if(ng.ne.ng2) then
              write(irstout) aflx
              write(irstout) xmom
              write(irstout) ymom
              write(irstout) hflxf
              write(irstout) zmom
              write(irstout) fhflx
            endif       
          endif     
        endif
  
! 2014_09_04 . scb        
        write(irstout)   flxlevel
        write(irstout)   xeave, smave
        do k=1,nz
          do l=1,nxy
            write(irstout)  rni(l,k),rnxe(l,k),rnpm(l,k),rnsm(l,k)
            write(irstout)  rnip(l,k),rnxep(l,k),rnpmp(l,k),rnsmp(l,k)
          enddo
        enddo
!! 2014_09_16 . scb for dbg
!        print *, rnip(18,2),rnxep(18,2),rnpmp(18,2),rnsmp(18,2)
!        pause
! !added end
      endif      

      return
    end

    
      subroutine rstedithex(ivaropt)
!!
!! Write restart information to file for hexagonal geometry.
!#include <param.h>
!      
!#include <geom.h>
!#include <geomh.h>
!#include <thgeom.h>
!#include <adj.h>
!#include <cntrod.h>
!#include <coolth.h>
!#include <defsfc.h>
!#include <lscoefh.h>
!#include <defpnt.h>
!#include <deffg.h>
!#include <decayht.h>
!#include <fbvar.h>
!#include <fuelth.h>
!#include <kin.h>
!#include <lscoef.h>
!#include <partrod.h>
!#include <power.h>
!#include <solvec.h>
!#include <ssvar.h>
!#include <xesm.h>
!#include <xsec.h>
!!
!#include <implicit.h>
!      include 'cntl.h'
!      include 'itrcntl.h'
!      include 'trcntl.h'
!      include 'files.h'
!      character*80 amesg80
!!
!      integer ivaropt
!! ======================================================================
!! Steady-State Initialization Restart Data
!! ======================================================================
!      if (ssinit) then    ! temp
!        write(amesg80,800)irstadv,sstime
!        call message(FALSE,TRUE,amesg80) 
! 800    format(11x,"Restart advancement: ",i4                               &
!     &,      " written at time: ",f10.3," s")
!        write(irstout) irstadv,sstime
!! defsfc.h variables
!        write(irstout) dhat
!        write(irstout) dhatz
!        write(irstout) wtdhat
!        write(irstout) dsum
!        write(irstout) betaphis
!        write(irstout) betaphisz
!! itrcntl.h variables
!        write(irstout) eigv
!! partrod.h variables
!        write(irstout) fweight
!! solvec.h variables
!        write(irstout) phi
!        write(irstout) psi
!        write(irstout) psid
!        write(irstout) psidd
!! xsec.h variables
!        write(irstout) xsnf
!        write(irstout) xstd
!! defpnt.h variables
!        write(irstout) pflx
!! ======================================================================
!! Transient Restart Data
!! ======================================================================
!      else
!!--- Write time-independent restart variables.
!        if (ivaropt.eq.0) then
!! adj.h variables
!          write(irstout) phia
!          write(irstout) sumv0
!! cntl.h variables
!          write(irstout) flxlevel
!! decayht.h variables
!          write(irstout) omalphatot
!! itrcntl.h variables
!          write(irstout) eigv, psil10
!! kin.h variables
!          write(irstout) nprec
!! ssvar.h variables
!          write(irstout) tdopl0
!          write(irstout) tcool0
!          write(irstout) dcool0
!          write(irstout) ppml0
!          write(irstout) crbdens0
!          write(irstout) fweight0
!! ssvar.h variables
!          write(irstout) xsxea0
!          write(irstout) xssma0
!          write(irstout) rnxe0
!          write(irstout) rnsm0
!!--- Write time-dependent restart variables.
!        else
!          write(amesg80,810)irstadv,time
!          call message(FALSE,TRUE,amesg80) 
! 810      format(11x,"Restart advancement: ",i4                         &
!     &,        " written at time: ",f10.3," s")
!! trcntl.h variables: irstadv, time
!          write(irstout) irstadv, time, cetadelt
!! cntl.h variables
!          if (.not.extth) then
!            write(irstout) ppm
!          endif
!! cntrod.h variables
!          write(irstout) crbpos
!! coolth.h variables
!          if (.not.extth) then
!            write(irstout) hcool
!            write(irstout) qeff
!            write(irstout) rhou
!            write(irstout) rhohu
!            write(irstout) u
!            write(irstout) ud
!          endif
!! defsfc.h variables
!          write(irstout) dhat
!          write(irstout) dhatz
!          write(irstout) wtdhat
!          write(irstout) dsum
!          write(irstout) betaphis
!          write(irstout) betaphisz
!! decayht.h variables
!          write(irstout) precconc
!          write(irstout) ezetdelt
!          write(irstout) omexpprod
!! fbvar.h variables
!          if (.not.extth) then
!            write(irstout) tdopl
!            write(irstout) tcool
!            write(irstout) dcool
!          endif
!! fuelth.h variables
!          if (.not.extth) then
!            write(irstout) tfuel
!            write(irstout) htcoef
!            write(irstout) qvol
!          endif
!! kin.h variables
!          write(irstout) prec
!          write(irstout) psit
!          write(irstout) psitd
!! lscoef.h variables
!          write(irstout) af
!! partrod.h variables
!          write(irstout) fweight
!! power.h variables
!          write(irstout) relp
!          write(irstout) plevel, plevel0, pleveld, rtdblrd
!! solvec.h variables
!!
!! store aphi for theta method
!          call axb(phi,aphi)
!!
!          write(irstout) phi
!          write(irstout) phid
!          write(irstout) aphi
!          write(irstout) phif
!          write(irstout) aphif
!          write(irstout) psi
!          write(irstout) psid
!          write(irstout) psidd
!! trcntl.h variables
!          write(irstout) delt, deltmd, trip, scram, tripbeg, tscrmbeg
!          write(irstout) pscrmbeg
!! xesm.h variables
!          write(irstout) rnxe
!          write(irstout) rni
!          write(irstout) rnsm
!          write(irstout) rnpm
!! xsec.h variables
!          write(irstout) xsnf
!          write(irstout) xstd
!          write(irstout) rvdelt
!          write(irstout) betap
!! lscoefh.h variables
!          write(irstout) hflx
!! defpnt.h variables
!          write(irstout) pflx
!          if(multigroup) then
!! deffg.h variables
!            write(irstout) aflx
!            write(irstout) xmom
!            write(irstout) ymom
!            write(irstout) hflxf
!            write(irstout) phiadjf
!            write(irstout) zmom1
!            write(irstout) zmom2
!            write(irstout) fhflx
!            write(irstout) fcnto
!            write(irstout) fcntzo
!! solvec
!            write(irstout) phif
!          endif
!        endif
!      endif
!      return
      end
