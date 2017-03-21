    !subroutine editout
    subroutine editout(iftran,time2)    ! 2014_10_06 . scb 
!
      use param
      use sfam,     only : psi
!
      include 'global.h'
      include 'files.h'
      include 'itrcntl.h'
      include 'xsec.h'
      include 'geom.h'
      include 'ffdm.h'
      include 'times.h'
      include 'editout.h'
      include 'pow.h'
      include 'ff.h'
      include 'pinpr.h'
      include 'thgeom.inc'   ! added in ARTOS ver. 0.2 . 2012_07_05 by SCB     
      include 'xesm.h'    ! 2014_01_10 . scb
      include 'trancntl.inc'   ! 2015_07_07 . scb
      
      integer :: nkf   ! 2012_08_22 . scb
      
! 2014_10_06 . scb
      logical :: iftran
      real    :: time2
      
      real    :: fr=0.d0, fz=0.d0, fxyz=0.d0
! added end

      logical,save :: firstout3=.true.   ! 2015_08_21 . scb
      
! assembly wise flux and power
      do ka=1,nza
        do la=1,nxya
          psia(la,ka)=0
          powa(la,ka)=0
          
          psia1=0
          powa1=0
          do m2=1,ng2
            phia1=0
            do k=kfmb(ka),kfme(ka)
              do li=1,latol(la)%nfm
                l=latol(la)%fm(li)
                volf=volnode(l,k)
                do m=mgb(m2),mge(m2)
                  phia1=phia1+phif(m,l,k)*volf
                  powa1=powa1+phif(m,l,k)*volf*xskpf(m,l,k)
                enddo
!                psia1=psia1+psif(l,k)
                psia1=psia1+psi(l,k)  ! 2012_09_27 . scb
              enddo
            enddo 
            !phia(m2,la,ka)=phia1/volnodea(la,ka)
            phia(m2,la,ka)=phia1/volnodea(la,ka)*flxlevel  ! 2014_01_10 . scb
          enddo !of m2
          !psia(la,ka)=psia1/volnodea(la,ka)
          !powa(la,ka)=powa1/volnodea(la,ka)
          psia(la,ka)=psia1/volnodea(la,ka)*flxlevel  ! 2014_01_10 . scb
          powa(la,ka)=powa1/volnodea(la,ka)*flxlevel  ! 2014_01_10 . scb
          !if(powa(la,ka) .gt. fxyz)  fxyz=powa(la,ka)    ! 2014_10_06 . scb
        enddo
      enddo   

! 2015_07_07 . scb for axial offset calculation
!============== calculate total power in each axial node=========
      do ka=kfbega,kfenda
        powat(1,ka)=0.d0
        do la=1,nxya
          if(.not.iffuella(la)) cycle
          vol=volnodea(la,ka)
          powat(1,ka)=powat(1,ka)+powa(la,ka)*vol
        enddo
      enddo
!============== calculate half plane and overapped ratio... ====
      sumheight=0.d0
      halfheight=coreheight*0.5
      kfbega1=kfbega
      kfenda2=kfenda
      powbot=0.d0
      powtop=0.d0  
      ao=0.d0
      if(kfbega.ne.kfenda) then
        do ka=kfbega,kfenda
          sumheight=sumheight+hza(ka)
          dif=halfheight-sumheight
        
          if(dif .lt. hza(ka+1)) then
            kfenda1=ka;   kfbega2=ka+2
            ratio=dif/hza(ka+1)
            exit                    
          endif
        enddo
!============ calculate total power in bottom and top... =======
        do ka=kfbega1,kfenda1
          powbot=powbot+powat(1,ka)
        enddo
        do ka=kfbega2,kfenda2
          powtop=powtop+powat(1,ka)
        enddo
        ka=kfenda1+1
        powbot=powbot+powat(1,ka)*ratio
        powtop=powtop+powat(1,ka)*(1.d0-ratio)
      
        ao=(powtop-powbot)/(powbot+powtop)
      endif
! added end


! 2014_10_06 . scb
      if(.not. iftran .and. flagout(3) .and. firstout3) then
#ifndef CVF
        write(filename(4), '("out/",a,".peak")') trim(caseid)  ! 2014_10_06 . scb
#else      
        write(filename(4), '(a,".peak")') trim(caseid)   ! 2014_10_06 . scb
#endif      
        open(io4,file=trim(filename(4)),status='unknown')
        write(io4,'(4a12)') 'Time(sec)','Fz','Fr','Fxyz'
        firstout3 = .false.
      endif
      
      inquire(file=filename(4),exist=filexist)
      if(.not. flagout(3))  filexist=.false.
      
      fr=0.d0
      fz=0.d0
      fxyz=0.d0
! added end

      if(iftran .and. .not. flagouttr) return   ! 2015_08_21 . scb
      
      write(io8,'(/a)') 'Core Average Axial Power Distribution'
      totpow=0.       ! 2012_08_22 . scb
      do ka=kfbega,kfenda
        powat(1,ka)=0
        volpl=0
        do la=1,nxya
          if(.not.iffuella(la)) cycle !if(.not.iffuela(iassytyp(la))) cycle
          vol=volnodea(la,ka)
          powat(1,ka)=powat(1,ka)+powa(la,ka)*vol
          volpl=volpl+vol
        enddo
        powat(1,ka)=powat(1,ka)/volpl
! 2012_08_22 . scb
        totpow=totpow+powat(1,ka)
      enddo
      nkf=kfenda-kfbega+1
      do ka=kfbega,kfenda
        powat(1,ka)=powat(1,ka)/totpow*nkf
! added end        
        write(io8,'(i3,1p,e15.4)') ka,powat(1,ka)
        if(powat(1,ka) .gt. fz)  fz=powat(1,ka)   ! 2014_10_06 . scb
      enddo
!
      write(io8,'(/a)') 'Assembly Power Distribution'
      totpow=0.
      do la=1,nxya
        powat(la,1)=0
        volfa=0
        do ka=kfbega,kfenda
          vol=volnodea(la,ka)
          powat(la,1)=powat(la,1)+powa(la,ka)*vol
          volfa=volfa+vol
        enddo
        powat(la,1)=powat(la,1)/volfa
        totpow=totpow+powat(la,1)
        !totpow=totpow+powat(la,1)/volfa
      enddo
      do la=1,nxya
! 2012_08_22 . SCB     
        powat(la,1)=powat(la,1)/totpow*nfuela
!        powat(la,1)=powat(la,1)/totpow*nchan

        if(powat(la,1) .gt. fr)  fr=powat(la,1)   ! 2014_10_06 . scb
! added end
      enddo
      call radmap1(powat(1,1),1,1,nxya,io8)
      
! 2014_10_06 . scb
      powa1=0.
      l=0
      volpl=0.
      do ka=kfbega,kfenda
        do la=1,nxya
          if(.not.iffuella(la)) cycle 
          l=l+1
          powa1=powa1 + powa(la,ka)*volnodea(la,ka)
          volpl=volpl + volnodea(la,ka)
        enddo
      enddo
      powa=powa/powa1*volpl
      do ka=kfbega,kfenda
        do la=1,nxya
          if(powa(la,ka) .gt. fxyz)  fxyz=powa(la,ka)    ! 2014_10_06 . scb
        enddo
      enddo

      if(filexist)  write(io4,'(4f12.4)')  time2,fz,fr,fxyz    ! 2014_10_06 . scb           
! added end

      write(io8,'(/a)') 'Core Average Axial Fission Source Distribution'
      do ka=kfbega,kfenda
        psiat(1,ka)=0
        volpl=0
        do la=1,nxya
          if(.not.iffuella(la)) cycle !if(.not.iffuela(iassytyp(la))) cycle
          vol=volnodea(la,ka)
          psiat(1,ka)=psiat(1,ka)+psia(la,ka)*vol
          volpl=volpl+vol
        enddo
        psiat(1,ka)=psiat(1,ka)/volpl
!        write(io8,'(i3,f10.6)') ka,psiat(1,ka)
        write(io8,'(i3,1p,e15.4)') ka,psiat(1,ka)  ! 2012_09_27 . scb
      enddo
!
      write(io8,'(/a)') 'Assembly Fission Source Distribution'
      do la=1,nxya
        psiat(la,1)=0
        volfa=0
        do ka=kfbega,kfenda
          vol=volnodea(la,ka)
          psiat(la,1)=psiat(la,1)+psia(la,ka)*vol
          volfa=volfa+vol
        enddo
        psiat(la,1)=psiat(la,1)/volfa
      enddo
!      call radmap(psiat(1,1),1,1,nxya,io8)
      call radmap1(psiat(1,1),1,1,nxya,io8)  ! 2012_09_27 . scb
!
      write(io8,'(/a)') 'Assemblywise Fission Source Distribution'
      do ka=kfbega,kfenda
        write(io8,'(a,i3)') 'Plane',ka
!        call radmap(psia(1,ka),1,1,nxya,io8)
        call radmap1(psia(1,ka),1,1,nxya,io8)  ! 2012_09_27 . scb
      enddo
!!
      do m=1,ng2
        write(io8,'(/a,i3)') 'Assemblywise Flux Distribution for 2 Group Condensing - Group ',m
        do ka=kfbega,kfenda
          write(io8,'(a,i3)') 'Plane',ka
          call radmap1(phia(1,1,ka),2,m,nxya,io8)
        enddo 
      enddo

      if(flagout(1))  write(io9,'(/a)') 'Mesh Fission Source Distribution'
      psi=psi*flxlevel
      do k=kfbeg,kfend
        if(flagout(1))  write(io9,'(a,i3)') 'Plane',k
!        call radmap(psif(1,k),1,1,nxy,io9)
        !call radmap(psi(1,k),1,1,nxy,io9)  ! 2012_09_27 . scb
        if(flagout(1))  call radmap(psi(1,k),1,1,nxy,io9)  ! 2014_01_10 . scb
      enddo
      psi=psi/flxlevel
      
      return
    end subroutine
!=======================================================================================!
    subroutine radmap(fx,n1,i1,nrad,iodev)
!
      use param
!
      include 'global.h'
      include 'files.h'
      include 'geom.h'
!
      character(4) maptype
      real :: fx(n1,nrad)
!
      if(nrad.eq.nxy) then
        do j=1,ny
          if(nxs(j).ne.1) then
            write(mesg,'(a,i3,a,i3,a)') '(',nxs(j)-1,'(10x),1p,',nrowx(j),'e10.3)'
          else
            write(mesg,'(a,i3,a)') '(1p,',nrowx(j),'e10.3)'
          endif
!          write(iodev,mesg) (fx(i1,nodel(i,j)),i=nxs(j),nxe(j))  ! 2015.03.28 pkb - Add a Comment
        enddo
      elseif(nrad.eq.nxya) then
        do ja=1,nya
          if(nxsa(ja).ne.1) then
            write(mesg,'(a,i3,a,i3,a)') '(',nxsa(ja)-1,'(8x),',nrowxa(ja),'f8.4)'
          else
            write(mesg,'(a,i3,a)') '(',nrowxa(ja),'f8.4)'
          endif
          write(iodev,mesg) (fx(i1,nodela(ia,ja)),ia=nxsa(ja),nxea(ja))
        enddo
      endif
!
      return
    end subroutine
!=======================================================================================!
    subroutine radmap1(fx,n1,i1,nrad,iodev)
!
      use param
!
      include 'global.h'
      include 'files.h'
      include 'geom.h'
!
      character(4) maptype
      real :: fx(n1,nrad)
!
      if(nrad.eq.nxy) then
        do j=1,ny
          if(nxs(j).ne.1) then
            write(mesg,'(a,i3,a,i3,a)') '(',nxs(j)-1,'(12x),1p,',nrowx(j),'e12.5)'
          else
            write(mesg,'(a,i3,a)') '(1p,',nrowx(j),'e12.5)'
          endif
          write(iodev,mesg) (fx(i1,nodel(i,j)),i=nxs(j),nxe(j))           ! 2016.9.21. jjh
        enddo
      elseif(nrad.eq.nxya) then
        do ja=1,nya
          if(nxsa(ja).ne.1) then
            write(mesg,'(a,i3,a,i3,a)') '(',nxsa(ja)-1,'(12x),1p,',nrowxa(ja),'e12.5)'
          else
            write(mesg,'(a,i3,a)') '(1p,',nrowxa(ja),'e12.5)'
          endif
          write(iodev,mesg) (fx(i1,nodela(ia,ja)),ia=nxsa(ja),nxea(ja))
        enddo
      endif
      return
    end subroutine
