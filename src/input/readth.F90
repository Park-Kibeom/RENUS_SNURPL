    subroutine readth
! read T/H data
      use param

      include 'global.h'
      include 'geom.h'
      include 'cards.h'
      include 'files.h'
      include 'thfuel.inc'
      include 'thgeom.inc'
      include 'frzfdbk.inc'
      include 'thcntl.inc'
      include 'thop.inc'
! 
! set defaults
! t/h mesh same as neutronic mesh
      do ia=1,nxa 
        nthx(ia)=nmeshx(ia)
      enddo
      do ja=1,nya
        nthy(ja)=nmeshy(ja)
      enddo
      nzth=nz
      do k=1,nz
!        junb(i)=k     ! 2012_08_22 . scb
        junb(k)=k     ! 2012_08_22 . scb
      enddo

      rhof=10282.                   !density of fuel kg/m^3
      arcpfuel(0)=162.3*rhof
      arcpfuel(1)=0.3038*rhof
      arcpfuel(2)=-2.391e-4*rhof
      arcpfuel(3)=6.404e-8*rhof
      rhoc=6600.                    !density of Zr
      arcpclad(0)=252.54*rhoc
      arcpclad(1)=0.11474*rhoc
      arcpclad(2)=0
      arcpclad(3)=0 
!       
      ndataf=0
      iffile=FALSE
      indev=io5
      
      if(probe.eq.DOT)  probe=''   ! 2014_12_17 . scb
!
  100 continue
      do while (probe.ne.DOT)
        read(indev,'(a512)',end=1000) oneline
        write(io8,'(a)') trim(oneline)
        if(probe.eq.BANG .or. oneline.eq.BLANK .or.ifnumeric(oneline)) cycle
        if(probe.eq.DOT .or. probe.eq.SLASH) exit
        if(probe.ne.BLANK) then
          backspace(indev)
          backspace(io8)
          go to 2000
        endif
        read(oneline,*) cardname
        call toupper(cardname)
        if(cardname.eq.'FILE') then
          indev=io5+100
          call openlf(indev,oneline)
          iffile=TRUE
          go to 100
        endif
        ndataf=nfields(oneline)-1
        select case(cardname) 
! added in ARTOS ver. 0.2 ( for LFR ). 2012_07_03 by SCB
          case('FUEL_COMP')
            read(oneline,*) cardname,wzr,wpu
! added end
          case('N_PINGT')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            read(oneline,*) cardname,(idum(i),i=1,ndataf)
            if(idum(1).ne.0) npint=idum(1)
            if(idum(2).ne.0) ngt=idum(2)
          case('FA_POWPIT')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            if(fdum(1).ne.0.) powfa0=fdum(1)
            if(fdum(2).ne.0.) pfa0=fdum(2)
          case('PIN_DIM')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            if(fdum(1).ne.0.) rs=fdum(1)
            if(fdum(2).ne.0.) rw=fdum(2)
            if(fdum(3).ne.0.) tw=fdum(3)
            if(fdum(4).ne.0.) rgt=fdum(4)
          case('FLOW_COND')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            if(fdum(1).ne.0.) tin=fdum(1)
            if(fdum(2).ne.0.) cfrperfa=fdum(2)
          case('GAMMA_FRAC')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            if(ndataf.ge.1) then
              if(fdum(1).ge.0 .and. fdum(1).lt.1) then
                fracdc=fdum(1)
              else
                mesg='Coolant Gamma Fraction < 0.0 or > 1.0'
                call terminate(mesg)
              endif
            endif
            if(ndataf.ge.2) then
              if(fdum(2).ge.0 .and. fdum(2).lt.1) then
                fracdvb=fdum(2)
              else
                mesg='Vessel Bypass Gamma Fraction < 0.0 or > 1.0'
                call terminate(mesg)
              endif
            endif
            if(ndataf.eq.3) then
              if(fdum(3).ge.0 .and. fdum(3).lt.1) then
                fracdwr=fdum(3)
              else
                mesg='Water Rod Gamma Fraction < 0.0 or > 1.0'
                call terminate(mesg)
              endif
            endif
          case('HGAP')
            read(oneline,*) cardname,hgap          
          case('N_RING')
            read(oneline,*) cardname,nr
          case('THMESH_X')
            read(oneline,*) cardname,(nthx(i),i=1,nxa)
          case('THMESH_Y')
            read(oneline,*) cardname,(nthy(i),i=1,nya)
          case('THMESH_Z')
            nzth=ndataf
            junb(0)=0
            read(oneline,*) cardname,(junb(i),i=1,nzth)          
          case('FREEZE_TF')
            frozentf=0
            read(oneline,*) cardname,freezetf
            if(ndataf.eq.2 .and.freezetf) then
              read(oneline,*) cardname,freezetf,frozentf
              frozentf=sqrt(frozentf+273.15)
            endif          
          case('FREEZE_DM')
            frozendm=0
            read(oneline,*) cardname,freezedm
            if(ndataf.eq.2 .and.freezedm) then
              read(oneline,*) cardname,freezedm,frozendm
              frozendm=frozendm*1000
            endif
          case('WRITE_FBV')
            read(oneline,*) cardname,writefbv
            if(ndataf.gt.1) call getfn(oneline,3,filename(17))
            io17=17
! open binary output FB variable file
            if(writefbv) call openfile(io17,FALSE,TRUE,filename(17))
          case('READ_DOPL')
            read(oneline,*) cardname,readdopl
            if(readdopl) then
              freezetf=TRUE
              call getfn(oneline,3,filename(21))
              io21=21
! open binary input FB variable file
!             call openfile(io21,TRUE,TRUE,filename(21))
            else
              freezetf=FALSE
            endif
          case('READ_TMDM')
            read(oneline,*) cardname,readtmdm
            if(readtmdm) then
              freezedm=TRUE
              call getfn(oneline,3,filename(22))
              io22=22
! open binary feedback variable file
!             call openfile(io22,TRUE,TRUE,filename(22))
            else
              freezedm=FALSE
            endif
          case('EFF_DOPLT')
            read(oneline,*) cardname,ifeffdopl
            if(ifeffdopl .and. ndataf.ne.1) then
              read(oneline,*) cardname,ifeffdopl,wfsurf
              wfcl=1-wfsurf
            endif
          case('KCOND_FUEL')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            do i=1,ndataf
              akfuel(i-1)=fdum(i)
            enddo
            do i=ndataf,4
              akfuel(i)=0
            enddo
          case('RHOCP_FUEL')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            do i=1,ndataf
              arcpfuel(i-1)=fdum(i)
            enddo
            do i=ndataf,3
              arcpfuel(i)=0
            enddo
          case('KCOND_CLAD')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            do i=1,ndataf
              akclad(i-1)=fdum(i)
            enddo
            do i=ndataf,3
              akclad(i)=0
            enddo
          case('RHOCP_CLAD')
            read(oneline,*) cardname,(fdum(i),i=1,ndataf)
            do i=1,ndataf
              arcpclad(i-1)=fdum(i)
            enddo
            do i=ndataf,3
              arcpclad(i)=0
            enddo
          case('UNIF_TH')
             read(oneline,*) cardname,(fdum(i),i=1,ndataf)
             unifdm=fdum(1)*1000
             uniftf=fdum(2)+CKELVIN
             uniftm=fdum(3)                      
          case('K_RATIO')   ! 2015_08_03 . scb
             read(oneline,*) cardname,(kratio(i),i=1,nassytyp)  ! 2015_08_03 . scb
          case default
            call terminate(trim(cardname)//' Card Not Allowed')
        end select
      enddo
      
 1000 continue

! reset the input device after done with local file
      if(iffile) then
        close(indev)
        indev=io5
        iffile=FALSE
! return to the next card in the input file
        go to 100
      endif
!
2000  continue
      
      nrp1=nr+1
      nrp2=nr+2
      nrp3=nr+3
      nrp4=nr+4
      nrp5=nr+5
      
      return
      
    end subroutine