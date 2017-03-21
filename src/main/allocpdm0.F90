    subroutine allocpdm0
!
      use param
      use allocs
      use xsec,   only : xstfsp3   ! 2014_10_07 . scb
      use Mod_FixedSource
      use Masterxsl
!
      include 'global.h'
      include 'itrcntl.h'
      include 'xsec.h'
      include 'geom.h'
      include 'geomh.h'    ! added in ARTOS ver. 0.2 ( for hexagonal ). 2012_07_05 by SCB
      include 'ffdm.h'
      !include 'ntbytes.h'
      include 'editout.h'
      include 'ff.h'
      include 'thcntl.inc'
      include 'perturb.inc'
!      include 'mslb.inc' ! MSLB
      include 'trancntl.inc' ! added in ARTOS ver. 0.2 ( for decay heat ). 2012_07_05 by SCB
      include 'thexpan.inc'  ! added in ARTOS ver. 0.2 ( for thermal expansion ). 2012_07_05 by SCB
      include 'xesm.h'   ! 2013_09_26 . scb

!
      nrdir2=nrdir*2
      nzap1=nza+1
!      nsurf=1       ! commented in ARTOS ver. 0.2 . 2012_07_05 by SCB
!
! geom varaibles
      call dmalloc(hxa,nxa)
      call dmalloc(hya,nya)
      call dmalloc(hza,nza)
      
! added in ARTOS ver. 0.2 . 2012_07_05 by SCB
      call dmalloc(hxa0,nxa)
      call dmalloc(hya0,nya)
      call dmalloc(hza0,nza)
! added end
!
      call dmalloc(iassytyp,nxya)
      call dmalloc(icompnumz,nza,nassytyp)
      call dmalloc(iffuela,nassytyp)
      call dmalloc(iffuella,nxya)
      call dmalloc(nxsa,nya)
      call dmalloc(nxea,nya)
      call dmalloc(nrowxa,nya)
      call dmalloc(nysa,nxa)
      call dmalloc(nyea,nxa)
      call dmalloc(nrowya,nxa)
      call dmalloc0(nodela,0,nxa+1,0,nya+1)
      call dmalloc(itoia,nx)
      call dmalloc(jtoja,ny)
      call dmalloc0(latoia,0,nxya)
      call dmalloc0(latoja,0,nxya)
      call dmalloc(neibz,2,nz)
      call dmalloc(nmeshx,nxa)
      call dmalloc(nmeshy,nya)
      call dmalloc(nmeshz,nza)
      nmeshx=1
      nmeshy=1
      nmeshz=1
      
      call dmalloc(nxsfa,nya) ! fuel region 
      call dmalloc(nxefa,nya) ! fuel region 
      call dmalloc(nxsf,ny) ! fuel region 
      call dmalloc(nxef,ny) ! fuel region 
      !nxefa=-1;nxef=-1   ! 2014_12_18 . scb commented

!control rod
      call dmalloc(irodtyp,nxya)
      call dmalloc(rodfullpos,nrodtyp)
      call dmalloc(rodstep,nrodtyp)
      call dmalloc(rodstep0,nrodtyp)   ! 2014_12_22 . scb
      call dmalloc(rodstepd,nrodtyp)      ! 2012_09_28 . scb
      call dmalloc(rodstepsize,nrodtyp)
      call dmalloc(rodstepsize0,nrodtyp)   ! added in ARTOS ver. 0.2 ( scram ). 2012_07_05 by SCB
      !call dmalloc(rodfrac,nz,nrodtyp)
      call dmalloc0(rodfrac,1,nz,0,nrodtyp)   ! 2014_08_21 . scb
      call dmalloc(pscrmbeg,nrodtyp)      ! 2012_09_28 . scb
!
! added in ARTOS ver. 0.2 ( scram, stuck rod ). 2012_07_05 by SCB
      call dmalloc(srfrac,nz,nrodtyp)  ! SCRAM
      call dmalloc(srstep,nrodtyp)     ! SCRAM
      call dmalloc(lstroda,nrodtyp)    ! STUCK ROD
      call dmalloc(lstrodb,nrodtyp)    ! STUCK ROD
      call dmalloc0(lcrbptr,0,nrodtyp) ! STUCK ROD
      call dmalloc(rodtola,nxya)       ! STUCK ROD
! added end
!	
      call dmalloc0(ktoka,0,nz+1)
      call dmalloc(kfmb,nza)
      call dmalloc0(kfme,0,nza)
      call dmalloc(heleva,nza)
      call dmalloc(helev,nz)
!
      if(.not.associated(latol))  allocate(latol(nxya))   ! 2014_12_22 . scb added if statement
!
      call dmalloc(albr,ng)
      call dmalloc(albrl,ng)
      call dmalloc(albrr,ng)
      call dmalloc(albxl,ng)
      call dmalloc(albxr,ng)
      call dmalloc(albyl,ng)
      call dmalloc(albyr,ng)
      call dmalloc(albzb,ng)
      call dmalloc(albzt,ng)
!
      call dmalloc(nxs,ny)
      call dmalloc(nxe,ny)
      call dmalloc(nrowx,ny)
      call dmalloc(nys,nx)
      call dmalloc(nye,nx) 
      call dmalloc(nrowy,nx)
! added in ARTOS ver. 0.2 ( hex ). 2012_07_05 by SCB
!      call dmalloc0(nodel,0,nx+1,0,ny+1)
      if(rect) then
        call dmalloc0(nodel,0,nx+1,0,ny+1)
      else
        call dmalloc0(nodel,-3,nx+4,-1,ny+2)
      endif
! added end
!
! xsec variables
      call dmalloc(smsq,ng,ng)
      call dmalloc(sigt,ng,ncomp)
      call dmalloc(siga,ng,ncomp)
      call dmalloc(sigd,ng,ncomp)
      call dmalloc(sigtr,ng,ncomp)
      call dmalloc(signf,ng,ncomp)
      call dmalloc(sigkf,ng,ncomp)
      call dmalloc(sigf,ng,ncomp)   ! 2014_07_31 . scb
      call dmalloc(signu,ng,ncomp)  ! 2014_07_31 . scb
      call dmalloc(sigt_p3,ng,ncomp)  ! 2014_11_28 . scb
      call dmalloc(sigs,ng,ncomp)
      call dmalloc(sigchi,ng,ncomp)
      call dmalloc(sigphi,ng,ncomp)
      call dmalloc(sigkinf,ncomp)
! 2014_05_13 . scb      
      !call dmalloc(sigadf,nrdir2,ng,ncomp)
      !call dmalloc(sigcdf,nrdir2,ng,ncomp)
      !call dmalloc(sigmdf,nrdir2,ng,ncomp)
      if(rect) then
        call dmalloc(sigadf,nrdir2,ng,ncomp)
        call dmalloc(sigcdf,nrdir2,ng,ncomp)
        call dmalloc(sigmdf,nrdir2,ng,ncomp)  
      else
        call dmalloc(sigadf,6,ng,ncomp)
        call dmalloc(sigcdf,6,ng,ncomp)
        call dmalloc(sigmdf,6,ng,ncomp)  
      endif      
! added end      

      call dmalloc(iffuelc,ncomp)
! added in ARTOS ver. 0.2 ( control assembly ). 2012_07_05 by SCB
      call dmalloc(ifcontc,ncomp) ! CONTROL ASSEMBLY
      call dmalloc(icomptyp,ncomp)
! added end
      do icomp=1,ncomp
        do m=1,ng
          sigadf(:,m,icomp)=1.d0
          sigcdf(:,m,icomp)=1.d0
          sigmdf(:,m,icomp)=1.d0
        enddo
      enddo
!    - scattering matrix
      call dmalloc(sigsm,ng,ng,ncomp)
      call dmalloc(sigsms,ng,ncomp)
      call dmalloc(sigsme,ng,ncomp)
      call dmalloc(sigsmax,ncomp)   ! 2013_07_19 . scb
!
! dxsec variables
      call dmalloc(dsigt_a,NUMOFDXS,ng,ncomp)
      call dmalloc(dsigd_tr,NUMOFDXS,ng,ncomp)
      call dmalloc(dsignf,NUMOFDXS,ng,ncomp)
      call dmalloc(dsigkf,NUMOFDXS,ng,ncomp)
      call dmalloc(dsigs,NUMOFDXS,ng,ncomp)
      call dmalloc(dsigchi,NUMOFDXS,ng,ncomp)
      call dmalloc(dsigf,NUMOFDXS,ng,ncomp)   ! 2014_07_31 . scb
      call dmalloc(dsigt_p3,NUMOFDXS,ng,ncomp)   ! 2014_11_28 . scb
      call dmalloc(dsigsm,ng,NUMOFDXS,ng,ncomp)     

! origin dxsec variables 
      Allocate(dsigt_a0(NUMOFDXS,ng,ncomp))
      Allocate(dsigd_tr0(NUMOFDXS,ng,ncomp))
      Allocate(dsignf0(NUMOFDXS,ng,ncomp))
      Allocate(dsigkf0(NUMOFDXS,ng,ncomp))
      Allocate(dsigs0(NUMOFDXS,ng,ncomp))
      Allocate(dsigchi0(NUMOFDXS,ng,ncomp))
      Allocate(dsigf0(NUMOFDXS,ng,ncomp))
      Allocate(dsigt_p30(NUMOFDXS,ng,ncomp))
      Allocate(dsigsm0(ng,NUMOFDXS,ng,ncomp))
! transient xsec
! 2014_04_10 . scb
      call dmalloc(lmbdk0,nprec)
      call dmalloc(betak0,nprec)
! added end
      call dmalloc(lmbdk,nprec,ncomp)
      call dmalloc(betak,nprec,ncomp)
! 2014_04_10 . scb for precursor default
      lmbdk0(1)=0.0128
      lmbdk0(2)=0.0318
      lmbdk0(3)=0.119
      lmbdk0(4)=0.3181
      lmbdk0(5)=1.4027
      lmbdk0(6)=3.9286
      betak0(1)=0.0002584
      betak0(2)=0.00152
      betak0(3)=0.0013908
      betak0(4)=0.0030704
      betak0(5)=0.001102
      betak0(6)=0.0002584
      do icomp=1,ncomp
        lmbdk(:,icomp)=lmbdk0(:)
        betak(:,icomp)=betak0(:)
      enddo
! added end      

      call dmalloc(sigchid,ng,ncomp)
      call dmalloc(velof,ng,ncomp)

! ff.h
      call dmalloc(iffset,ncomp)            ! ff set number for each comp. 
      call dmalloc(ff,npin,npin,ng,nffset)  ! ff for unrodded set
      call dmalloc(ffr,npin,npin,ng,nffset) ! ff for rodded set

! perturb.inc
      call dmalloc(idpbank,npbank)
      call dmalloc(ntpbank,npbank)
      call dmalloc(tbank,maxntm,2,npbank)
      
      call allocth0

!FIXME remove
! scb : MOX bench 구현 안할 것. 삭제. 20120620
!#ifdef MOXBENCH
!      call allocmoxbch     
      if(ixsecver.eq.2) call allocmoxbch     ! 2015_07_31 . scb recover MOXBENCH
!#endif
! added in ARTOS ver. 0.2 ( hexagonal ). 2012_07_05 by SCB
      if(.not.rect) then
         call dmalloc(wtass,nxya)        !assembly weight
         call dmalloc(rhzbar2,2,nz)      !2/(hz(k)*(hz(k)+hz(k+1)))
         if(ng.ne.2) then
            call dmalloc(reflratf,ng)
            call dmalloc(reflratzbf,ng)
            call dmalloc(reflratztf,ng)
            call dmalloc(reflratzf,ng) 
            call dmalloc(alphazf,ng) 
         else
            reflratf=>reflrat 
            reflratzbf=>reflratzb
            reflratztf=>reflratzt
            reflratzf=>reflratz
         endif
         call alloctpen
      endif
! added end      

! 2013_09_27 . scb
! Xe/Sm      
      call dmalloc(gammafp,3,ncomp)         !comp-wise FP yields  
      call dmalloc(sigxea,ng,ncomp)       !Xe comp-wise micro xsec            
      call dmalloc(dsigxea,5,ng,ncomp)     !Xe feedback partial micro !bwrxsec 
      
      call dmalloc(dcontxea,ng,ncomp)       !Xe control rod partial micro       
      call dmalloc(sigsma,ng,ncomp)       !Sm comp-wise micro xsec            
      call dmalloc(dsigsma,5,ng,ncomp)     !Sm feedback partial micro !bwrxsec 
      call dmalloc(dcontsma,ng,ncomp)       !Sm control rod partial micro       
! added end 
! 2016_04_07~~2016_04_12 . pkb
      Allocate(iSrctyp(nxya))    
      Allocate(Srcdenza(nza,nSrctyp))
      Allocate(Srcdenz(nz,nSrctyp))
      Allocate(ExtSrcMap(nx,ny))
! Added End
      return
!
! allocate memory after node structure is determined
      entry allocpdm
!
! 
!    - node xsec
      call dmalloc(xstf,ng,nxy,nz)
! added in ARTOS ver. 0.2 ( hexagonal ). 2012_07_05 by SCB
      call dmalloc(xstrf,ng,nxy,nz)
      call dmalloc(xstfc,ng,nxy,nz)
      call dmalloc(xstfn,ng,nxy,nz)
! added end      
      !call dmalloc(xstrfdcsp,ng,nxy,nz)  ! 2012_09_28 . scb
      if(.not.associated(xsaf)) call dmalloc(xsaf,ng,nxy,nz)
      call dmalloc(xsdf,ng,nxy,nz)
      call dmalloc(xsdf2,ng,nxy,nz)   ! added in ARTOS ver. 0.2 . 2012_07_05 by SCB
      if(.not.associated(xsnff)) call dmalloc(xsnff,ng,nxy,nz)
      if(.not.associated(xsff)) call dmalloc(xsff,ng,nxy,nz)   ! 2013_10_02 . scb
      if(.not.associated(xsnuf)) call dmalloc(xsnuf,ng,nxy,nz)   ! 2014_10_24 . scb
      call dmalloc(xskpf,ng,nxy,nz)
      call dmalloc(xschif,ng,nxy,nz)
      call dmalloc(xbeta,ng,nxy,nz,ndirmax)
      call dmalloc(xbeta2,ng,nxy,nz,ndirmax)   ! 2013_07_15 . scb
      
! 2014_05_13 . scb      
      !call dmalloc(xsadf,nrdir2,ng,nxy,nz)
      !call dmalloc(xscdf,nrdir2,ng,nxy,nz)
      if(rect) then
        call dmalloc(xsadf,nrdir2,ng,nxy,nz)
        call dmalloc(xscdf,nrdir2,ng,nxy,nz)
      else
        call dmalloc(xsadf,6,ng,nxy,nz)
        call dmalloc(xscdf,6,ng,nxy,nz)
      endif     
! added end      
      call dmalloc(xssf,ng,ng,nxy,nz)
      call dmalloc(xssfs,ng,nxy,nz)
      call dmalloc(xssfe,ng,nxy,nz)
      call dmalloc(xss2nf,ng,nxy,nz)   ! added in ARTOS ver. 0.2 . 2012_07_05 by SCB
      !call dmalloc(xsmax,nxy,nz)   ! 2013_07_19 . scb
      if(.not.associated(xsmax)) call dmalloc(xsmax,nxy,nz)   ! 2013_07_19 . scb
! added in ARTOS  2014_02_17. pkb
      call dmalloc(xstf0,ng,nxy,nz)
      call dmalloc(xsaf0,ng,nxy,nz)
      call dmalloc(xsdf0,ng,nxy,nz)
      call dmalloc(xsnff0,ng,nxy,nz)
      call dmalloc(xskpf0,ng,nxy,nz)
      call dmalloc(xssf0,ng,ng,nxy,nz)
      call dmalloc(xsff0,ng,nxy,nz)
! added end
      
!
! ffdm varaibles
      if(.not.associated(phif)) call dmalloc0(phif,1,ng,0,nxy,0,nz)
      if(.not.associated(phifp))call dmalloc0(phifp,1,ng,0,nxy,0,nz) ! neutron flux from previous time step, 06.02.16.alc
      if(.not.associated(phifp_))call dmalloc0(phifp_,1,ng,0,nxy,0,nz) ! neutron flux from previous time step, 06.02.16.alc (it is necessary for calculate error between two iterations)
      call dmalloc0(phi,1,ng,0,nxy,0,nz)   ! added in ARTOS ver. 0.2 . 2012_07_05 by SCB
      call dmalloc(psif,nxy,nz)
! 
!editout.h 
      call dmalloc(powa,nxya,nza)
      call dmalloc(powat,nxya,nza)
      call dmalloc(psia,nxya,nza)
      call dmalloc(psiat,nxya,nza)
      call dmalloc(phia,ng2,nxya,nza)
!
!      nbytes=(nbytesf*NBF+nbytesi*NBI)/(2**10)   ! added in ARTOS ver. 0.2 . 2012_08_08 by SCB
      !nbytes=(nbytesf*8+nbytesi*4)/(2**10)   ! added in ARTOS ver. 0.2 . 2012_08_08 by SCB
      nbytes=(nbytesf*4)/(2**10)   ! 2014_07_21 . scb
      write(mesg,'(a,i8,a)') 'Allocated Memory:',nbytes,' KBytes'
      !write(mesg,'(a,f8.0,a)') 'Allocated Memory:',nbytes,' MBytes'
      call message(TRUE,TRUE,mesg)

!nodal.h
      call dmalloc(jnet,2,ng,nxy,nz,ndirmax)
      call dmalloc(phisfc,2,ng,nxy,nz,ndirmax)
! MSLB
!      if(ifmslb) then    
!         call dmalloc(velo,ng,nxy,nz)
!         call dmalloc(xstr,ng,nxy,nz)
!		do l=1,nxy
!			do k=1,nz
!				xssfs(1,l,k)=1
!				xssfe(1,l,k)=1
!				xssfs(2,l,k)=1
!				xssfe(2,l,k)=2							
!			enddo
!		enddo
!      endif
  
! 2013_09_26 . scb         
! Xe/Sm     
      call dmalloc(gami,nxy,nz)          !Iodine FP yield                    
      call dmalloc(gamxe,nxy,nz)          !Xenon FP yield                     
      call dmalloc(gampm,nxy,nz)          !Promethium FP yield                
      call dmalloc(rni,nxy,nz)          !Iodine Number Density              
      call dmalloc(rnxe,nxy,nz)          !Xenon Number Density               
      call dmalloc(rnpm,nxy,nz)          !Promethium Number Density          
      call dmalloc(rnsm,nxy,nz)          !Samarium Number Density            
      call dmalloc(rnip,nxy,nz)          !Previous Iodine Number Density     
      call dmalloc(rnxep,nxy,nz)          !Previous Xenon Number Density      
      call dmalloc(rnpmp,nxy,nz)          !Previous Promethium Number Density 
      call dmalloc(rnsmp,nxy,nz)          !Previous Samarium Number Density   

      call dmalloc(arate,nxy,nz)          ! absorbtion by xenon: phi * saxe
      call dmalloc(frate,nxy,nz)          ! fission rate: phi * Sf   
      call dmalloc(darate,nxy,nz)          ! derive of xenon absorbtion by time
      call dmalloc(dfrate,nxy,nz)          ! derive if fission rate by time

      call dmalloc(xsxea,ng2,nxy,nz)      !Xe node-wise micro xsec            
      call dmalloc(xssma,ng2,nxy,nz)      !Sm node-wise micro xsec       
      
      if(ng.ne.ng2) then
          call dmalloc(xsxeaf,ng,nxy,nz)        
          call dmalloc(xssmaf,ng,nxy,nz)
      else
          xsxeaf=>xsxea
          xssmaf=>xssma
      endif  
! added end

   
      if(ifsp3) call dmalloc(xstfsp3,ng,nxy,nz)   ! 2014_10_07 . scb    
                 
      return      
    end subroutine
