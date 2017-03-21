    subroutine runtr  
! transient calculation
      use param
      use allocs
      use trinx_cntl  ! added in ARTOS ver. 0.2 . 2012_07_10 by SCB. 
      use tran,     only : malloctran
      use sfam,     only : psisfam => psi, phisfam=>phi, reigv
      use sfam,     only : eigv0   ! 2014_09_05 . scb
      use bdf     ! 2012_09_05 . scb
      use tran,     only : deltm, deltmd
      use xsec,     only : xsnf,xsf   ! 2013_10_11 . scb
      use vtk,      only : flagvtk,ivtkfreq,ivtkstep  ! 2014_01_14 . scb
      USE MASTERXSL, ONLY : FLAGMASXSL     ! 2014_08_23 . SCB
      use linkdt,   only : mmidata   ! 2014_12_16 . scb
      use cmfd2g,     only : reigvs    ! 2014_12_22 . scb
      use dumping,  only : dumpfunc    !function for dump data into output file (with using of external function or by default)
      
      include 'global.h'
      include 'xsec.h'
      include 'pow.h'
      include 'files.h'
      include 'geom.h'
      include 'ffdm.h'
      include 'thcntl.inc'
      include 'thgeom.inc'
      include 'thfuel.inc'
      include 'thop.inc'
      include 'thfdbk.inc'
      include 'trancntl.inc'
      include 'pinpr.h'
! added in ARTOS ver. 0.2 . 2012_07_10 by SCB. 
!      include 'mslb.inc' ! MSLB
      include 'thlink.inc' ! FREK/MARS
      include 'thexpan.inc'
! added end      
      include 'cntl.h'  ! 2013_05_16 . scb for restart
      include 'xesm.h'  ! 2013_09_30 . scb
!#ifdef MOXBENCH
      include 'moxbch.inc'      ! 2015_07_31 . scb recover MOXBENCH
!#endif   
! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.    
      data dtsysa /0./ ! FREK/MARS
      save dtsysa      ! FREK/MARS
       
      logical :: flag2g=FALSE
      integer,save :: iskip=0
      integer :: count            ! added for xenon iterations .alc
      !integer,save :: iodbg=2001, iodbg2=2002  ! 2012_09_05 . scb
      !integer,save :: iorod=2003    ! 2014_09_22 . scb
      character*256 :: form, form2, form3*15            ! 2014_09_22 . scb
      logical,save :: first=TRUE
! added end      
      
      real,pointer,dimension(:,:,:) :: xschifd,xslmbdk,xsbetak,xsrvelof,prec
      character*80  :: chartime
      
      real :: peaktime=0. , peaklevel=0.      ! 2012_09_28 . scb      
      
! 2013_10_18 . scb      
      integer :: i_tunit=0   ! = 0 : sec
                             ! = 1 : min
                             ! = 2 : hr
                             ! = 3 : day
      real*8 :: timep     ! printing time                       
      
! added end                             

! 2014_09_27 . scb
      real :: timesrcdtdbg  ! 2014_09_27 . scb
      common / srcdtdbg / timesrcdtdbg
! added end      

#ifdef DLL_TASS
! 2014_12_16 . scb
      type(MMIDATA)  art2mmi2      ! mmi display data
      common / artosmmi / art2mmi2
! added end
#endif      
      logical :: firstdll   ! 2014_12_18 . scb      
      common / firstdata / firstdll   ! 2014_12_19 . scb
                 
      real,save :: tdbg = 0.d0  ! 2015_06_23 . scb
      integer,save :: itr=0   ! 2015_06_23 . scb
            
      flag2g=ng.eq.ng2
      nintot=0;noutbeg=0

! adjust fission xsec to make the core critical      
! 2013_10_28 . scb (below line is commented)
      !if(ilink.lt.0) then  ! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.     ! FREK/MARs  
        if(ixsecver.eq.2) then
!#ifdef MOXBENCH
          ! 2015_07_31 . scb recover MOXBENCH
          do ic=1,ncomp !dmm
            do m=1,ng
              sigbchnf(:,m,ic)=sigbchnf(:,m,ic)*reigv
            enddo
          enddo
          ! added end
!#else
          !call terminate('ERROR : Compile RENUS again with a preprocessor "MOXBENCH"')
!#endif
        else          
          do icomp=1,ncomp !dmm
            do m=1,ng !mg
              signf(m,icomp)=reigv*signf(m,icomp)
              signu(m,icomp)=reigv*signu(m,icomp)   ! 2014_07_31 . scb
              do icond=DPPM,NUMOFDXS
                dsignf(icond,m,icomp)=reigv*dsignf(icond,m,icomp)
              enddo
            enddo
          enddo
        endif
! 2013_10_28 . scb 
      !endif  ! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.     ! FREK/MARs  
      
      !reigvs=1.d0  ! 2014_12_22 . scb
      
! 2014_12_19 . scb      
      if(.not.first) then
        if(.not.firstdll) then
          firstdll=.true.
          !go to 1219
        else
          go to 1000  ! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.     ! FREK/MARs        
        endif
      endif
! added end      
      
! 2012_10_17 . scb
      do icomp=1,ncomp
        if(sigchid(1,icomp).gt.1.e-10) cycle
        
        sigchid(:,icomp)=sigchi(:,icomp)
      enddo
      
! 2014_01_22 . scb for normalizing delayed source fraction..      
      do icomp=1,ncomp
        sum=0.d0
        sumd=0.d0 ! 2014_01_22 . scb
        do m=1,ng
          sum=sum+sigchi(m,icomp)
          sumd=sumd+sigchid(m,icomp) ! 2014_01_22 . scb
        enddo
        if(sum.lt.1.e-30) sum=1.e-30   ! 2013_10_10 . scb 
        if(sumd.lt.1.e-30) sumd=1.e-30   ! 2014_01_22 . scb
        fnorm=1/sum
        fnormd=1/sumd   ! 2014_01_22 . scb
        do m=1,ng
          sigchi(m,icomp)=sigchi(m,icomp)*fnorm
          sigchid(m,icomp)=sigchid(m,icomp)*fnormd   ! 2014_01_22 . scb
        enddo  
      enddo  
! added end
      
! syncronize data
      !call dmalloc(xslmbdk,nprec,nxy,nz)
      !call dmalloc(xsbetak,nprec,nxy,nz)
      !!call dmalloc(prec,nprec,nxy,nz)  ! 2013_05_18 . scb
      !
      !call dmalloc(xschifd,ng,nxy,nz)
      !call dmalloc(xsrvelof,ng,nxy,nz)
      !if(.not.associated(xslmbdk))  allocate(xslmbdk(nprec,nxy,nz))
      !if(.not.associated(xsbetak))  allocate(xsbetak(nprec,nxy,nz))
      !if(.not.associated(xschifd))  allocate(xschifd(nprec,nxy,nz))
      !if(.not.associated(xsrvelof))  allocate(xsrvelof(nprec,nxy,nz))
      !if(first) then   ! 2014_12_18 . scb
      allocate(xslmbdk(nprec,nxy,nz))
      allocate(xsbetak(nprec,nxy,nz))
      allocate(xschifd(nprec,nxy,nz))
      allocate(xsrvelof(nprec,nxy,nz))
      xslmbdk=0.d0
      xsbetak=0.d0
      xschifd=0.d0
      xsrvelof=0.d0
      !endif      
        
      do k=1,nz
        ka=ktoka(k)
        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          icomp=icompnumz(ka,iat)
          psif(l,k)=0
          do m=1, ng
            xsnff(m,l,k)=reigv*xsnff(m,l,k) 
            !xsff(m,l,k)=xsnff(m,l,k)/2.3    ! 2013_10_02 . scb
            xschifd(m,l,k)=sigchid(m,icomp)
            xsrvelof(m,l,k)=1/velof(m,icomp)
            psif(l,k)=psif(l,k)+xsnff(m,l,k)*phif(m,l,k)
          enddo
          
! 2013_10_11 . scb
          if(ng.ne.ng2) then
            do m=1, ng2
              xsnf(m,l,k)=reigv*xsnf(m,l,k) 
              !xsf(m,l,k)=xsnf(m,l,k)/2.3    ! 2013_10_02 . scb
            enddo
          endif          
! added end
          psif(l,k)=psif(l,k)*volnode(l,k)
          do mp=1, nprec
            xsbetak(mp,l,k)=betak(mp,icomp)
            xslmbdk(mp,l,k)=lmbdk(mp,icomp)
          enddo
        enddo
      enddo
! added end
      
      !reigv0=reigv  ! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.
      reigv0=1.d0/eigv0  ! 2014_09_15 . scb
      
      if(fdbk) then     ! 2015_07_10 . scb  for dbg
          kgap=thetaf*hgap*delr
          kgapb=thetafb*hgap*delr
          tmp1=hgap*tw*(4-tw/rg)*rs/rg
          kgap2=thetaf*hgap*tw*rs/rg
          kgap4=thetaf*tmp1
          kgap4b=thetafb*tmp1
      endif         
        
      do k=1,nz
        do l=1,nxy
          psisfam(l,k) = psif(l,k)
        enddo
      enddo

! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.
!      if(ifmslb) call initmslb ! MSLB
      if(decayht) call initdecay(deltm0) ! DECAY HEAT
! added end      
           
      call malloctran(ng,nxy,nz,nprec,xslmbdk,xsbetak,xschifd,xsrvelof,prec)
      !call inittran()  ! 2013_05_20 . scb
      
! 2014_08_23 . scb tr initialization
      IF(flagmasxsl)  CALL INITTRAN_MASXS(REIGV0)

      !call inittran(rstrt)  ! 2013_05_20 . scb
      call inittran(rstrt,deltm0)  ! 2013_05_20 . scb
      
! 2012_09_05 . scb
! 2013_05_18 . scb . moved to subroutine init
      !if(flagbdf) then
      !  call initbdf(unit_tm, bdforder, ng, nxy, nz)
      !  if(autodt) then
      !    call initautodt(deltm0,bdforder,ng,nxy,nz)
      !    
      !    write(filename(2), '("out/",a,"_delt.plt")') trim(caseid)
      !    open(iodbg2, file=trim(filename(2)),status='unknown')
      !  endif
      !endif
! added end  
! added end
      
      deltmd=deltm0
      deltm=deltm0     ! 2012_09_05 . scb
! 2013_05_20 . scb      
      !time = 0.0
      !isteptr=0
! 2014_04_28 . scb     
      !if(.not.rstrt) then
      !  time = 0.d0
      !  isteptr=0
      !endif 
      if(isteptr.eq.0) time=0.d0
! added end      
! added end            
      first=.false.  ! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.      
      
! 2012_09_05 . scb . move the output file generation routine from the time loop to here...
! 2013_12_16 . scb added : #ifndef
      
1219  continue  ! 2014_12_19 . scb
      
! 2014_12_29 . scb      
      if(.not.rstrt) then
        isteptr = 0
        time=0.d0   ! 2014_12_19 . scb for dbg
      endif      
! added end      

      if(.not. autodt) flagout(5)=.false.   ! 2014_10_15 . scb
#ifndef CVF
      res = makedirqq('out/')
! 2014_10_15 . scb      
      !if(rect) then
      !  write(filename(1), '("out/",a,"_"i2.2,"box.plt")') trim(caseid),(nx/nxa)**2
      !else
      !  write(filename(1), '("out/",a,".plt")') trim(caseid)
      !endif
      if(flagout(2))  write(filename(1), '("out/",a,".plt")') trim(caseid)
! added end      
      if(flagout(5))  write(filename(2), '("out/",a,".dt")') trim(caseid)
      if(flagout(4))  write(filename(3), '("out/",a,".rod")') trim(caseid)   ! 2014_09_22 . scb
#else
! 2014_10_15 . scb
      !if(rect) then
      !  write(filename(1), '(a,"_"i2.2,"box.plt")') trim(caseid),(nx/nxa)**2
      !else
      !  write(filename(1), '(a,".plt")') trim(caseid)   ! 2014_08_08 . scb
      !endif
      if(flagout(2))  write(filename(1), '(a,".plt")') trim(caseid)   ! 2014_08_08 . scb
! added end      
      if(flagout(5))  write(filename(2), '(a,".dt")') trim(caseid)  ! 2014_08_08 . scb
      if(flagout(4))  write(filename(3), '(a,".rod")') trim(caseid)   ! 2014_09_22 . scb
#endif      
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
! added end
      if(flagout(2))  open(io1, file=trim(filename(1)),status='unknown')
      if(flagout(5))  open(io2, file=trim(filename(2)),status='unknown')
      if(flagout(4))  open(io3, file=trim(filename(3)),status='unknown')   ! 2014_09_22 . scb
! 2013_10_02 . scb      
      !write(iodbg,'(a10,1p,a15,0p,3a15)') 'Time(sec)','Power(%)','TFMAX(K)','TFAVG(K)','TMAVG(K)'
      !write(iodbg,'(f10.5,1p,e15.4,0p,3f15.2,1p,2e15.4,0p,3i4)') time,plevel*100,tfmax+CKELVIN,tfuelavg,tcoolavg
! 2013_10_18 . scb      
      !write(iodbg,'(a10,1p,a15,0p,4a15)') 'Time(sec)','Power(%)','TFMAX(K)','TFAVG(K)','TMAVG(K)','Avg ND of Xe'
      !write(iodbg,'(f10.5,1p,e15.5,0p,3f15.2,1p,e15.4,0p,3i4)') time,plevel*100,tfmax+CKELVIN,tfuelavg,tcoolavg,xeave      
      !if(deltm0.lt.1) then
      i_tunit=0
      !write(iodbg,'(a15,1p,a15,0p,6a15)') 'Time(sec)','Power(%)','TFMAX(C)','TFAVG(K)','TMAVG(K)','DMAVG(kg/m3)','Avg ND of Xe','Avg ND of Sm'
      !if(flagout(2))  write(io1,'(a15,1p,a15,0p,7a15)') 'Time(sec)','Power(%)','TFMAX(C)','TFAVG(C)','TMMAX(C)','TMAVG(C)','DMAVG(kg/m3)','Avg ND of Xe','Avg ND of Sm'
      if(flagout(2))  write(io1,'(a15,1p,a15,0p,8a15)') 'Time(sec)','Power(%)','TFMAX(C)','TFAVG(C)','TMMAX(C)','TMAVG(C)','DMAVG(kg/m3)','Avg ND of Xe','Avg ND of Sm','Axial Offset'
        
      !elseif(deltm0.lt.60) then
      !  i_tunit=1
      !  write(iodbg,'(a15,1p,a15,0p,5a15)') 'Time(min)','Power(%)','TFMAX(K)','TFAVG(K)','TMAVG(K)','Avg ND of Xe','Avg ND of Sm'
      !elseif(deltm0.lt.3600) then
      !  i_tunit=2
      !  write(iodbg,'(a15,1p,a15,0p,5a15)') 'Time(hr)','Power(%)','TFMAX(K)','TFAVG(K)','TMAVG(K)','Avg ND of Xe','Avg ND of Sm'
      !elseif(deltm0.lt.86400) then        
      !  i_tunit=3
      !  write(iodbg,'(a15,1p,a15,0p,5a15)') 'Time(day)','Power(%)','TFMAX(K)','TFAVG(K)','TMAVG(K)','Avg ND of Xe','Avg ND of Sm'
      !endif             
        
      !write(iodbg,'(f15.5,1p,e15.5,0p,4f15.4,1p,2e15.4)') time,plevel*100,tfmax+CKELVIN,tfuelavg,tcoolavg,dcoolavg,xeave,smave 
      !if(flagout(2))  write(io1,'(f15.5,1p,e15.5E3,0p,5f15.4,1p,2e15.4)') time,plevel*100,tfmax,tfuelavg-CKELVIN,tmmax,tcoolavg,dcoolavg,xeave,smave   ! 2014_10_06 . scb
      if(flagout(2))  write(io1,'(f15.5,1p,e15.5E3,0p,5f15.4,1p,2e15.4,E15.4)') time,plevel*100,tfmax,tfuelavg-CKELVIN,tmmax,tcoolavg,dcoolavg,xeave,smave,ao   ! 2015_07_07 . scb
! added end      
! added end      

! 2014_09_22 . scb      
      if(flagout(4))  then
        write(form,'("(a15,",i2,"a15)")')  nrodtyp
        form2='      Time(sec)'
        do irod=1,nrodtyp
          write(form3,'(" CR_type_",i2.2,"(cm)")') irod
          form2=trim(form2)//form3
        enddo           
        write(form,'("(a",i3.3,")")')  (nrodtyp+1)*15
        write(io3,form)   trim(form2)      

! 2014_09_22 . scb
        write(form,'("(",i2,"f15.4)")') nrodtyp+1
        write(io3,form) 0.d0,rodstep(1:nrodtyp)
! added end          
      endif
! added end      
              
1000  continue   ! added in ARTOS ver. 0.2 . 2012_07_10 by SCB. ! FREK/MARS    

      if(rsted)  call rstedit(0)  ! 2013_05_16 . scb for restart
#ifdef DLL_TASS   ! 2014_04_10 . scb
      if(rsted) then
        call rstedit(1)
        call restart(2)
        return
      endif      
#endif                                 

      do while(TRUE)
        if(ilink.ge.0) goto 410   ! 2014_11_06 . scb
        
! 2013_05_16 . scb for restart
        if(irstfreq.eq.0) then
          goto 410
        elseif(rsted .and. mod(itimeadv,irstfreq).eq.0) then
          call rstedit(1)
          irstadv=irstadv+1
        endif
        
410     continue
        itimeadv=itimeadv+1
! added end        

        isteptr=isteptr+1
! time expansion        
        if(automain)  goto 500   ! 2012_09_27 . scb
        
        if(time.lt.tswitch) then
          deltm=deltm0
        elseif(time.le.tswitch2*(1-1e-7)) then
          deltm=deltm0*texpand
! 2012_08_22 . scb          
!        else
!          deltm=deltm0*texpand2
        elseif(time.le.tswitch3*(1-1e-7)) then
          deltm=deltm0*texpand2
        elseif(time.le.tswitch4*(1-1e-7)) then
          deltm=deltm0*texpand3            
        else
          deltm=deltm0*texpand4
! added end
        endif

500     continue  ! 2012_09_05 . scb
! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.
#ifdef DLL
        if(ilink.ge.0) then
          dtsysa=dtsysa+dtsys
          if(dtsysa.gt.deltm*0.95 .or. rsted) then
            deltm=dtsysa
            dtsysa=0.0
            iconv=1
            write(mesg,'(a)') '========================================================================================='
            call message(true,true,mesg)
            write(mesg,'(a10,i10,a10,f15.6)')  'skip :: ', iskip, 'time :: ', time+deltm
            call message(FALSE,true,mesg) 
          else
            isteptr=isteptr-1
            iskip=iskip+1
            iconv=0
            return
          endif
        else
          write(mesg,'(a)') '========================================================================================='
          call message(true,true,mesg)
        endif
#endif        
! added end
                
        !if(flagbdf)  call updbdf(isteptr, bdforder, autodt, deltm, ng, nxy, nz)     ! 2012_09_05 . scb
        
        !itr=itr+1
        !t0=dclock()
        if(flagbdf)  call updbdf(isteptr, bdforder, autodt, deltm, ng, nxy, nz,ifsp3)     ! 2013_08_12 . scb   
        !t1=dclock()
        !tdbg = tdbg + t1 - t0
        !write(623,*) itr, tdbg    ! 2015_06_23 . scb
        
        time=time+deltm
        timesrcdtdbg=time  ! 2014_09_27 . scb
        
! 2013_10_18 . scb
        if(i_tunit.eq.0) timep=time
        if(i_tunit.eq.1) timep=time/60
        if(i_tunit.eq.2) timep=time/3600
        if(i_tunit.eq.3) timep=time/86400

! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.
!        if(fdbk) call drivetftr(deltm)

        if(fdbk .and. ilink.lt.0) call drivetftr(deltm)

! 2015_07_10 . scb         for dbg
!        if(fdbk .and. ilink.lt.0) then
        !  if(deltm0.lt.30) then
        !    !call drivetftr(deltm)
!            call drivethss
        !  else
        !    call drivethss
        !  endif
!        endif        
! added end          
! added end

! adjust pertubation
! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.
!        call perturb(time)
        call perturb(deltm) ! SCRAM
! added end

        ! additional iterations between neutrons and poisons equations
        ! 2016_06_02.alc
          if(flagxesm) call calxesm(true)  ! 2013_09_30 . scb

          call runtrtm(flag2g,deltm)

        
        if(automain) then
          call detcont(ng,bdforder,deltm,errtol,safefac)
          
          if(cont>1.0.and.imaxiter.ge.maxiter.and.deltm.ne.maxdeltm) then
            if(cont>maxbound(ibound))  cont=maxbound(ibound)
            
            deltmsave=cont*deltm
                        
            if(ibound.lt.3) then
              ibound=ibound+1
            else
              maxiter=maxiter-1
            endif
            
            deltmupd=TRUE
            deltmsave=aint(deltmsave/unit_tm)*unit_tm
            if(deltmsave.lt.unit_tm)   deltmsave=unit_tm
            if(deltmsave.gt.maxdeltm)  deltmsave=maxdeltm
            imaxiter=0
          elseif(errtol<error) then
            imaxiter=0
            ibound=1
            if(maxiter.le.bdforder)  maxiter=bdforder
            if(retry) then
              time=time-deltm
              deltm=0.8*deltm
              deltm=aint(deltm/unit_tm)*unit_tm
              if(deltm.lt.unit_tm) deltm=unit_tm
              deltmsave=deltm
              updneed=FALSE
              if(scram) then
                do id=1,nrodtyp              
                  rodstep(id) = rodstepd(id)
                enddo
              endif
              goto 500
            else
              if(cont<0.8) then
                deltmsave=0.8*deltm
              else
                deltmsave=cont*deltm
              endif
              deltmupd=TRUE
              deltmsave=aint(deltmsave/unit_tm)*unit_tm
              if(deltmsave.lt.unit_tm) deltmsave=unit_tm
            endif
          endif
        endif
        
        call updprec   ! 2012_09_26 . scb
        
        call updrelpow(TRUE)
        if(decayht) call upddecay(deltm) ! added in ARTOS ver. 0.2 . 2012_07_10 by SCB. ! DECAY HEAT
        
        call signal(deltm)   ! 2012_09_28 . scb
        
! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.
!        if(fdbk) call drivethtr(deltm)
        if(fdbk .and. ilink.lt.0) call drivethtr(deltm)   ! 2015_07_10 . scb
!        if(fdbk .and. ilink.lt.0) call drivethss

!        plevel = 1.D0

!        call drivethss

! added end
        if(flagxesm) call updxesm   ! 2014_07_29 . scb

        fq =0
        fr =0
        lamax = 0;
        ipmax = 0;
        jpmax = 0;
        if(pinpower .and. mod(time, 0.1).eq.0.0) then       
            call driveppr(TRUE)
            fq = 0
            do k=1,nz
              do la=1,nxya
                if(pppeak(la,k) .gt. fq) then
                    fq = pppeak(la,k)
                endif
              enddo
            enddo
            
            fr = 0
            ipmax = 0
            jpmax = 0
            lamax = 0
            do la=1,nxya
              do jp=1,npin
                do ip=1,npin
                  if(powvalr(ip,jp,la) .gt. fr) then
                    fr = powvalr(ip,jp,la)
                    ipmax = ip
                    jpmax = jp
                    lamax = la
                  endif
                enddo
              enddo
            enddo          
        endif

! 2012_09_26 . scb        
        if(flagout(5))  write(io2,'(2f20.8)') time,deltm
        
! added in ARTOS ver. 0.2 . 2012_07_10 by SCB.
!        write(mesg,'(a,1p,e15.4,0p,a,f6.3)') 'Power : ', plevel*100, '%   Time : ', time
        !write(mesg,'(a20,1p,e15.4,0p,a20,f10.3)') 'Power : ', plevel*100, '%           Time : ', time
        write(mesg,'(a20,1p,e15.4,0p,a20,f15.3)') 'Power : ', plevel*100, '%           Time : ', time
! added end

        call message(true,true,mesg)
        
! 2012_10_12 . scb        
!        call editout  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
        !if(flagouttr) then    ! 2015_06_22 . scb
        if(rect) then
          !call editout  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
          call editout( .true. , time )   ! 2014_10_06 . scb
! 2014_01_14 . scb          
          if(flagvtk .and. ivtkfreq.ne.0) then
            ivtkstep=ivtkstep+1
            if(mod(ivtkstep,ivtkfreq).eq.0) then
              ivtkstep=0
              call makevtk(false)
            endif   
          endif
! added end          
        else
          !call editout_hex
          call editout_hex(.true. , time)   ! 2014_10_06 . scb
        endif
        !endif        
! added end        

! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
!        write(iodbg,'(f8.5,1p,e12.4,0p,2f12.2,1p,2e12.4,0p,3i4)')  &
!              time,plevel*100,tfmax,toutavg,fq/plevel,fr,lamax,ipmax,jpmax
 !       write(iodbg,'(f10.5,1p,e15.4,0p,3f15.2,1p,2e15.4,0p,3i4)') time,plevel*100,tfmax+CKELVIN,tfuelavg,tcoolavg
 !!    +,              fq/plevel,fr,lamax,ipmax,jpmax
        !write(iodbg,'(f10.5,1p,e15.5,0p,3f15.2,1p,e15.4,0p,3i4)') time,plevel*100,tfmax+CKELVIN,tfuelavg,tcoolavg,xeave     ! 2013_10_04 . scb
        !write(iodbg,'(f15.5,1p,e15.5,0p,4f15.4,1p,2e15.4)') time,plevel*100,tfmax+CKELVIN,tfuelavg,tcoolavg,dcoolavg,xeave,smave 
        !if(flagout(2))  write(io1,'(f15.5,1p,e15.5E3,0p,5f15.4,1p,2e15.4)') time,plevel*100,tfmax,tfuelavg-CKELVIN,tmmax,tcoolavg,dcoolavg,xeave,smave     ! 2014_10_06 . scb
        if(flagout(2))  write(io1,'(f15.5,1p,e15.5E3,0p,5f15.4,1p,2e15.4,e15.4)') time,plevel*100,tfmax,tfuelavg-CKELVIN,tmmax,tcoolavg,dcoolavg,xeave,smave,ao     ! 2015_07_07 . scb

! 2014_09_22 . scb
        if(flagout(4)) then
          write(form,'("(",i2,"f15.4)")') nrodtyp+1
          write(io3,form) time,rodstep(1:nrodtyp)
        endif        
! added end              
 
        if(ilink.le.0 .and. plevel.gt.peaklevel) then
          peaklevel=plevel
          peaktime=time
          !peaktime=timep  ! 2013_10_18 . scb
        endif
! added end
        
        deltmd=deltm
        
! 2012_09_26 . scb
        if((tend-time).le.1.2*deltmsave) then
          deltmsave=tend-time
          deltm=deltmsave
        endif
! added end        
        
! 2013_05_16 . scb for restart
        !if(ilink.ge.0) return ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.  ! FREK/MARS
        !
        !if(time .ge. tend*(1-1e-7)) exit 
        
        if(ilink.ge.0) then         
! 2014_12_16 . scb         
#ifdef DLL_TASS
          art2mmi2%powcore = plevel*100
          !art2mmi2%toutcore = int(tcoolavg)
#endif         
! added end
          if(rsted) then
            call rstedit(0)   ! 2014_04_28 . scb
            call rstedit(1)
            call restart(2)
          endif
          
          return ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.  ! FREK/MARS
        endif
                
        if(time .ge. tend*(1-1e-7)) then
          if(.not.flagouttr)  call editout( .false. , time )  ! 2015_07_07 . scb
          if(rsted) call rstedit(1)          
          exit 
        endif 

        call dumpfunc() ! save state arrays into output file, 2016_06_17 . alc

! added end        
      enddo 

! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      if(ilink.le.0) then
        print '(/a12,f15.6)','  Delt   :: ',deltm
        print '(a12,f15.6)', '  Ptime  :: ',peaktime
        print '(a12,f15.6)', '  Plevel :: ',peaklevel

! 2012_11_21 . scb        
        write(io8,'(/a12,f15.6)'),'  Delt   :: ',deltm
        write(io8,'(a12,f15.6)'), '  Ptime  :: ',peaktime
        write(io8,'(a12,f15.6)'), '  Plevel :: ',peaklevel
! added end        
      endif
! added end

      if(rsted .or. rstrt) call restart(2)  ! 2013_05_16 . scb
             
      if(rect .and. flagvtk)  call makevtk(true)  ! 2014_01_14 . scb   
	
      return
    end subroutine
