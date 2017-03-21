! added in ARTOS ver. 0.2. 2012_07_24 by SCB 
!#ifdef DLL_MARS
!    Subroutine master(iflag,nthdim,fbdata,bankpos,deltat,iok,p3dmaster,dnbmaster)     									
!      !MS$ ATTRIBUTES DLLEXPORT :: master 
!      use msflib            ! file handling constructs
!      use linkdt            ! derived types for Neutronics/T-H linking
!      use wordsize
#ifdef DLL_TASS
		!Subroutine hectos(iflag,nthdim,fbdata,bankpos,deltat,rstflag,rstfname,art2mmi,iok)
		Subroutine hectos(iflag,nthdim,fbdata,bankpos,deltat,rstflag,rstfname,art2mmi,iok)
      !MS$ ATTRIBUTES DLLEXPORT :: HECTOS
      !use msflib            ! file handling constructs
      use wordsize
		  use linkdt            ! derived types for ARTOS/T-H linking
! added end      
#else
    program ARTOS 
#endif
      use param
      use timer
      use allocs
      use geom,      only : setgeom
      use xsec,      only : setxsec
      use sfam,      only : mallocsfam, initsfam,eigv
      use sfam,      only : setcntl   ! 2012_11_08 . scb
      use sfam,      only : eigv0  ! 2014_09_04 . scb
      use itrinfo,   only : tnodal,tcmfd,tcmfd2g,tother,nnodal,ncmfd2g,ncmfdmg
! TPEN
      use geomhex,   only : setgeomhex
      use tpen,      only : malloctpen, inittpen
      use cmfdhex2g, only : mallochex2g
      use hexfdm,    only : malloc_hexfdm
      use hexfdm3,   only : malloc_hexfdm3
      use tpen_sp3,      only : malloctpen_sp3, inittpen_sp3
      use cmfdhex2g_sp3, only : mallochex2g_sp3
      use vtk,        only : flagvtk   ! 2014_01_14 . scb
      !use tran,       only : deltm0    ! 2014_09_03 . scb
      
      use readN2A           ! 2016.9.5. jjh
      
      use MASTERXSL   ! 2014_05_22 . pkb
      USE DEPLMOD     ! 2014_10_06 . PKB

      include 'global.h'
      include 'itrcntl.h'
      include 'ffdm.h'
      include 'xsec.h'
      include 'geom.h'
      include 'nodal.h'
      include 'pinpr.h'
      include 'times.h'
      include 'cntl.h'
      include 'files.h'  ! 2014_04_10 . scb
      include 'ff.h'
      include 'srchppm.h'
      include 'thcntl.inc'
      include 'thlink.inc' ! FREK/MARS
      include 'thexpan.inc'
! hexnodal
      include 'geomh.h'
      include 'geomhfc.h'
      include 'defhex.h'
      include 'defsfc.h'
      include 'defpnt.h'
      include 'deffg.h'
      include 'lscoefh.h'
      include 'hexfdm.h'
      include 'trancntl.inc'   ! 2014_08_29 . scb
!
      include 'pow.h'  ! 2014_12_17 . scb
      
      integer :: ncmfdmgtot=0,ncmfd2gtot=0,ninmgtot=0,nin2gtot=0
      logical,save :: first=TRUE
      logical,save :: first_test = .true.  ! 2014_12_16 . scb for dbg
      
      real,save :: eigvsave=1.d0   ! 2014_09_04 . scb
      
      logical :: firstdll = .true.   ! 2014_12_18 . scb
      logical,save :: firstdll2 = .true.  ! 2014_12_19 . scb for dbg
      
      common / firstdata / firstdll   ! 2014_12_19 . scb

#ifdef DLL
      integer :: iflag  ! =0 for initialization
                         ! =1 for steady-state advancement
                         ! =2 for transient advancement (normal)
                         ! =3 for transient advancement after trip
                         ! =4 for wrapping up ARTOS calculation
                         ! <0 for standalone ARTOS execution 
      type(THDIM) nthdim  ! array dimension data for T/H-side variables
      type(FBVAR) fbdata(0:*) ! T/H feedback data
                              ! * index 0 is for core average data
      !real ::  bankpos(*)    ! bank positions in cm withdrawn 
      real(NBF) ::  bankpos(*)    ! bank positions in cm withdrawn       
#endif
#ifdef DLL_MARS  
      real ::  p3dmaster(*)  ! 3D normalized power distribution 
      real ::  dnbmaster(*)  ! 3D normalized power distribution 
#endif
#ifdef DLL_TASS
      logical  rstflag(2)         ! restart flag
                                  ! rstflag(1) = true for read
                                  ! rstflag(2) = true for write
      character(*)  rstfname(2)   ! restart file name
          
      type(MMIDATA)  art2mmi      ! mmi display data
      type(MMIDATA)  art2mmi2       ! 2014_12_16 . scb
      common / artosmmi / art2mmi2  ! 2014_12_16 . scb      
#endif      
#ifdef DLL
      real(4) :: deltat     ! time step size, s
      integer :: iok       ! return condition flag
                            ! =0 for neutronic-T/H inconsistency in input data
                            ! =1 for OK sign to proceed
                            ! =2 for SS convergence achieved in ARTOS
                            ! =3 for convergence problem in ARTOS

      logical ifin
      save timemas
      data timemas/0.0/
      integer,save :: iskip=0
      integer,parameter :: imod=1
      real,save :: dtsum=0.d0    ! 2014_08_29 . scb
      
      integer :: iflagskip  ! 2013_09_17 . scb
      integer :: ic, icomp
      
      integer,save :: istepdll=0
      
      logical,save :: firstout=.true.  ! 2014_04_15 . scb
      
      call cpu_time(tstart)

      ilink=iflag
      dtsys=deltat
      
#ifdef DLL_TASS   ! 2014_04_10 . scb
      rstrt=rstflag(1)
      rsted=rstflag(2)
      if(rstrt)  filename(9)=trim(rstfname(1))//'.rst'
      if(rsted)  filename(10)=trim(rstfname(2))//'.rst'
      
      istepdll=istepdll+1
      !write(24,*) 'istep in DLL : ',istepdll, 'ilink : ',ilink, 'rstrt : ',rstflag(1), 'rsted : ',rstflag(2)
      !write(22,*) 'ilink : ',ilink, 'Bank Position : ', bankpos(1:5)
#endif      

! 2014_08_29 . scb
      !if(ilink.eq.1) then
      !  iskip=iskip+1
      !  iflagskip=mod(iskip,imod)  ! 2013_09_17 . scb
      !  if(iflagskip.eq.1) then
      !    write(mesg,'(a18,i10)')       'skip :: ', iskip
      !    call message(FALSE,true,mesg) 
      !    return  ! 2014_04_10 . scb
      !  endif
      !  !if(iflagskip.ne.1) return  ! 2014_04_10 . scb commented
      !endif
      !if(firstdll) then
      !  write(mesg,'(a28)')       'ARTOS initialization'
      !  goto 829
      !endif
      
      if(ilink.eq.1) then
        iskip=iskip+1
        dtsum=dtsum+dtsys
        if(dtsum.gt.deltm0) then
          iskip=0
          dtsum=0.d0
          write(mesg,'(a28)')       'ARTOS calculation'
          call message(FALSE,true,mesg)          
        else
          write(mesg,'(a18,i10)')       'ARTOS skip : ', iskip
          call message(FALSE,true,mesg)          
          return
        endif        
      endif
! added end      

829   continue
! incoming data processing
      ifin=TRUE
      if(ilink.gt.0) then
        call exfbdata(ifin,nthdim,fbdata,bankpos,iok)
        ifin=FALSE
      endif
      
      if(ilink.eq.0) then
        nmars=0
        sstime=0.d0
        trtime=0.d0
      endif
   
      nmars=nmars+1

      select case(ilink)
        case(-1) ! standalone FREK execution
          goto 1000 
        case(0) ! initialization 
          nssstep=0
          ntrstep=0
          goto 1000
        case(1) ! steady-state advancement
          sstime=sstime+dtsys
          nssstep=nssstep+1
          goto 2000
        case(2,3) ! transient advancement (normal) or after trip
!#ifdef DLL_TASS
!          if(rstrt) then
!            call restart(0)
!            nmars=0
!            trtime=0.d0
!            
!            goto 1000
!          endif          
!#endif                     
#ifdef DLL_TASS   ! 2014_04_10 . scb
! 2014_09_11 . scb
          !if(rsted) then
          !  call restart(1)
          !  if(trtime.lt.1.e-30) then
          !    sstime=sstime+dtsys
          !    nssstep=nssstep+1
          !    goto 2000                         
          !  else
          !    trtime=trtime+dtsys
          !    ntrstep=ntrstep+1
          !    !goto 3000            
          !    goto 3500  ! 2014_04_10 . scb
          !  endif
          !endif
          if(rsted) then
            call restart(1)
            if(trtime.lt.1.e-30) then
              sstime=sstime+dtsys
              nssstep=nssstep+1
              goto 2000                         
            endif
          endif
          
          trtime=trtime+dtsys
          ntrstep=ntrstep+1
          goto 3500 
! added end           
#endif          
          trtime=trtime+dtsys
          ntrstep=ntrstep+1
          goto 3000
        case(4) ! wrapping up FREK calculation
! 2014_04_15 . scb          
          !if(firstout)  call editout_hex
          if(firstout)  call editout_hex(.true. , 0. )   ! 2014_10_06 . scb
          firstout=.false.   
! added end          
          goto 4000
      end select
#endif
 
1000 continue

      call timeron()   ! ttotal

! initialize
      call timeron()   ! tinit

      call preproc 
      
      if(firstdll) then  ! 2014_12_18 . scb
        call readinput
        call init			
      else
        plevel=plevel00    ! 2014_12_22 . scb
        plevel0=plevel00   ! 2014_12_22 . scb
        rodstep=rodstep0   ! 2014_12_22 . scb
        call updrodfrac   ! 2014_12_22 . scb
      endif          ! 2014_12_18 . scb
        
      if(ilink.eq.0) then
        call exfbdata(ifin,nthdim,fbdata,bankpos,iok) 
! 2014_04_23 . scb        
        call updrodfrac
        call updxsec(TRUE,false)
! added end        
        ifin=FALSE 
      endif

      if(firstdll) then  ! 2014_12_18 . scb
        call setcntl(ninmax,ifcmfd2g)  ! 2012_11_08 . scb
      
        call setgeom( ng,nx,ny,nz,nxy,ndir,symopt,isymang,                       &
                      isymloc,nxs,nxe,nys,nye,nxsf,nxef,                          &
                      kfbeg,kfend,jfbeg,jfend,nodel,                              &
                      hmesh,volnode,volcore,volfuel,                              &
                      (/albxl(1),albxr(1),albyl(1),albyr(1),albzb(1),albzt(1)/),  &
                      ltola,ktoka,iassytyp,iffuela,                               &
                      ifrecsp3,                                                   &   ! 2013_07_05 . scb
                      nxa,nya,nza,rhmesh,ltoi,ltoj)        ! 2013_07_15 . scb

        if(.not.rect) then 
          call setgeomhex( ng,nz,nxy,nassy,nxpnt,ncorn,nxsfc,                    &
                            nsurf,hf2f,ndivhs,kfbeg,kfend,                        &                   
                            volnode,hz,albr(1),albzb(1),albzt(1),                 &
                            wtass,neignd,neigjin,neigpt,neigz,                    &
                            neigsfc,neigsnd,neigsndz,                             &
                            codpnt,pbdv,wtdhat,neigsfcz,ineigcond,                &
                            ipntr,ilubnd,iastopnt,ifhexfdm,ifhexsp3,              &
                            ndiv,neigtri,neigtria)
          if(ifhexfdm) then
            if(ifhexsp3) then
              call malloc_hexfdm3
            else
              call malloc_hexfdm
            endif
          else
            if(ifhexsp3) then
              call malloctpen_sp3
              call mallochex2g_sp3
            endif
            call malloctpen
            call mallochex2g
          endif
        endif

        call setxsec(ng,mgb(2),nxy,nz,xstf,xsaf,xsdf,xsnff,      &
                     xskpf,xschif,xsff,xbeta,xsadf,xscdf,xssf,xssfs,   &
                     xssfe,xss2nf,xsdf2,xstrf,xbeta2,xsmax)   
        ! 2013_07_15 . scb :: xbeta2 added
        ! 2013_07_19 . scb :: xsmax added
        ! 2013_10_02 . scb :: xsff added

        call dmalloc(curil,ndirmax,nxy,nz,ng)
        call dmalloc(curir,ndirmax,nxy,nz,ng)
        call dmalloc(curol,ndirmax,nxy,nz,ng)
        call dmalloc(curor,ndirmax,nxy,nz,ng)

        call mallocsfam(ng,nxy,nz,phif,phisfc,jnet,curol,curor,curil,curir)  

  ! 2013_07_16 . scb      
        if(ifrecsp3) then
          call malloccmfd(ng,nxy,nz,nzp1,nsurf)
          call initsp3senm(ng)
        endif
  ! added end      

        if(.not.ifhexfdm) then
          call initsfam(TRUE)
          if(.not.rect) then   
            if(ifhexsp3) then
              call inittpen_sp3
            else
              call inittpen
            endif 
          endif
        endif
      endif    ! 2014_12_18 . scb       

      call timeroff(tinit)

2000  continue      

! 2014_11_01 . pkb for depletion

! added end

! 2013_05_16 . scb for restart     
      if(rstrt) then
        call rstinp
        call restart(2)    ! 2014_12_02 . scb
#ifdef DLL_TASS  ! 2014_04_10 . scb
        call exfbdata(TRUE,nthdim,fbdata,bankpos,iok)     
        eigvsave=eigv
        !eigv0 = eigv   ! 2014_09_04 . scb
        art2mmi2%powcore = plevel*100  ! 2014_12_17 . scb
        
        goto 5000
#else
        goto 3500
#endif  
      endif

! standard code
      call timeron()    ! 2015_06_22 . scb . tsteady
      call runss(ncmfdmgtot,ninmgtot,ncmfd2gtot,nin2gtot)
      call timeroff(tsteady)
                    
      eigvsave=eigv 
      eigv0 = eigv   ! 2014_09_04 . scb
#ifdef DLL      
      fbdata(0)%dm=1.d0   ! 2014_04_11 . scb      
#endif      
      
      if(ilink.ge.0) then
        call exfbdata(ifin,nthdim,fbdata,bankpos,iok)
        !call editout_hex   ! 2014_04_16 . scb
        iok=1 !iconv
        write(io8,*) "iok & iconv :", iok, iconv
        if(nmars.eq.1) then
          ifin=FALSE
        endif
#ifdef DLL_TASS   ! 2014_04_10 . scb
        art2mmi2%powcore = plevel*100  ! 2014_12_17 . scb

        if(rsted) goto 3000
#endif           
        ! 2014_04_14 . scb for dbg

        
        !write(1223, *) tinavg, toutavg, art2mmi2%tincore, art2mmi2%toutcore

        go to 5000   
      endif

      if(ifreact) call cal_beta

3000  continue

      !if(first) then  
      if(first .or. .not.firstdll) then   ! 2014_12_22 . scb
        call normalize  
        first=FALSE    ! 2017_03_08 . bys
! 2014_01_10 . scb. Move the output processing routine here...
        if(rect) then
          call editout(.false. , 0. )   ! 2014_10_06 . scb
          if(flagvtk) call makevtk(false)   ! 2014_01_14 . scb
      
          !stop 'vtk file is made successfully !!'
        else
#ifdef DLL          
          !if(firstout)  call editout_hex
          !if(firstout)  call editout_hex(.true. , 0. )   ! 2014_10_06 . scb
          !firstout=.false.   ! 2014_04_15 . scb
#else
          !call editout_hex
          call editout_hex(.true. , 0. )   ! 2014_10_06 . scb
#endif          
        endif
        
        if(printff) call writeff

! pin power calculation
        if(pinpower) call driveppr(FALSE)
! added end   
      else  ! depletion 
         ! print( ID
          call editout(.false. , 0. )   ! 2014_10_06 . scb
        
      endif

#ifdef DLL_TASS   ! 2014_04_10 . scb call transient only for restart
      if(ilink.ge.0) then
        if(trtime.lt.1.e-30 .and. rsted) then
          call runtr
          goto 5000
        endif        
      endif      
#endif
      If(ifdepl) then  
        DO ID=1,Ndepl
          If(FirstDepl) then
            Call Allocdepl
            Call NuclideChainDefinition
            Call Mas2Dep
          Endif
          If(depltyp.eq.0) then
            If(firstDepl) then
              Call Edit_Depl(id,inpbu(id),eigv)
            Else
              Call Edit_Depl(id,inpbu(id-1),eigv)
            Endif
          Elseif(depltyp.eq.1) then
            call Edit_Depl(id,burnn(1,1),eigv)
          Endif
          Call Normalization_depl
          Call DriveDepletion(id)
          FirstDepl = .FALSE.
        END DO
        If(depltyp.eq.0) then
          Call Edit_Depl(ndepl+1,inpbu(id-1),eigv)
        Elseif(depltyp.eq.1) then
          call Edit_Depl(ndepl+1,burnn(1,1),eigv)
        Endif
      End if
      
! transient calculation            
3500  continue         
    
      call timeron()   ! 2015_06_22 . scb . ttransient
      
      !if(ifdepl) call depletion()
      if(qtransient) call qtrans() !quasi-stationary calculations with xenon dynamic
      if(transient) call runtr

      call timeroff(ttransient)

      !if(ifdepl) call Depletion
      
#ifdef DLL      
      fbdata(0)%dm=-1.d0   ! 2014_04_11 . scb 
#endif      

      if(ilink.ge.0 .and. iconv.eq.0) go to 5000
      
      if(ilink.ge.0) then
        iok=1
        call exfbdata(ifin,nthdim,fbdata,bankpos,iok)
        !call editout_hex   ! 2014_04_16 . scb
        ifin=FALSE
        if(ilink.eq.2 .or. ilink.eq.3) then
          go to 5000
        elseif(ilink==4) then
          write(906,*) "Time for ARTOS",timemas
        endif
      endif
    
      call timeroff(ttotal) 

4000  continue     
        
! write results to file	
      !call cpu_time(tend)

      write(mesg, '("    ")')
      call message(FALSE,true,mesg) 
      write(mesg, '("============================================    Result    ============================================")')
      call message(FALSE,true,mesg) 
      !write(mesg,'(a18,f20.5)')       'K-EFF           : ', eigvsave
      write(mesg,'(a18,f20.5)')       'K-EFF           : ', eigv0   ! 2014_09_04 . scb
      call message(FALSE,true,mesg)
      if(srchppm .or. fdbk) then
        write(mesg,'(a18,f10.4)')       'BORON           : ', ppm
        call message(FALSE,true,mesg)
      endif

      negpnt=0
      do m=1, ng
        do k=1,nz
          do l=1,nxy
            if(phif(m,l,k) .lt. 0) negpnt=negpnt+1
          enddo
        enddo
      enddo

      write(mesg,'(a18,i8,a3,i8)') 'Negative Flux   : ',negpnt,'/',ng*nxy*nz
      call message(FALSE,true,mesg)

      write(mesg, '("    ")')
      write(mesg, '("======================================    Performance Summary    =====================================")')
      call message(FALSE,true,mesg) 
      write(mesg,'(a18,f10.3)')       'Initial    Time : ',tinit
      call message(FALSE,true,mesg)
      write(mesg,'(a18,f10.3)')       'XSEC       Time : ',txsec
      call message(FALSE,true,mesg)      
      write(mesg,'(a18,f10.3,a9,i4)') 'Nodal-Calc Time : ',tnodal, '##: ',nnodal
      call message(FALSE,true,mesg)

! 2012_09_27 . scb      
!      write(mesg,'(a18,f10.3,a9,i4,a9,i4)')   'CMFD-Calc  Time : ',tcmfd, '2G: ',ncmfd2g,'MG: ',ncmfdmg
!      call message(FALSE,true,mesg)
!      write(mesg,'(a18,f10.3,a9,i4,a9,i4)')   '        2G Time : ',tcmfd2g
!      call message(FALSE,true,mesg)	
      write(mesg,'(a18,f10.3,a9,i4,a9,i4)')   'CMFD MG    Time : ',tcmfd, 'MG: ',ncmfdmg
      call message(FALSE,true,mesg)
      write(mesg,'(a18,f10.3,a9,i4,a9,i4)')   'CMFD 2G    Time : ',tcmfd2g, '2G: ',ncmfd2g
      call message(FALSE,true,mesg)	
! added end      

      write(mesg,'(a18,f10.3)')       'T/H        Time : ',tth
      call message(FALSE,true,mesg)
      write(mesg,'(a18,f10.3)')       'PPR        Time : ',tppr
      call message(FALSE,true,mesg)      
      write(mesg,'(a18,f10.3)')       'OTHER      Time : ',tother
      call message(FALSE,true,mesg)      
      write(mesg, '("    ")')
      call message(FALSE,true,mesg)   
      write(mesg,'(a18,f10.3)')       'Steady     Time : ',tsteady
      call message(FALSE,true,mesg)
      write(mesg,'(a18,f10.3)')       'Transient  Time : ',ttransient
      call message(FALSE,true,mesg)      
      write(mesg, '("    ")')
      call message(FALSE,true,mesg)   
      write(mesg,'(a18,f10.3)')       'Total-CPU  Time : ',ttotal
      call message(FALSE,true,mesg)

      !if(ifrecsp3) call dbgsp3
      
! 2014_12_16 . scb
      close(io8)
      close(io9)
      close(io4)
      close(io1)
      close(io2)
      close(io3)
      close(io5)
       
      call dealloc
      
#ifdef DLL_TASS
      firstdll=.false.   ! 2014_12_18 . scb
#endif    

!#define debug_dll_rst
#ifdef debug_dll_rst
      if(firstdll .and. firstdll2) then
        firstdll = .false.
        firstdll2 = .false.
        goto 1000
      endif
#endif      
      
      
! added end      
      
5000  continue
! 2014_12_17 . scb      
#ifdef DLL_TASS
      !write(1223, *) tinavg, toutavg, art2mmi2%tincore, art2mmi2%toutcore
      if(trtime.lt.1.e-30)  call editout_hex(.false. , 0.)   ! 2014_12_22 . scb
      !write(1223, *) tinavg, toutavg, art2mmi2%tincore, art2mmi2%toutcore
      art2mmi = art2mmi2   
#ifdef mmi_print      
      write(1222,'(es12.5e2,6i10)') art2mmi%powcore,art2mmi%tincore,art2mmi%toutcore,art2mmi%peaklpd,art2mmi%ao,art2mmi%crbpos(1:2)     
#endif      
#endif    
! added end

    end
