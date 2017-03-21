! added in ARTOS ver. 0.2 (Hexagonal). 2012_07_06 by SCB   
    !subroutine editout_hex
    !subroutine editout_hex(iftran,time)    ! 2014_10_06 . scb 
    subroutine editout_hex(ifprint,time)    ! 2014_12_22 . scb  added ifprint
!
      use param
      use allocs
      use fluxyoyo, only : expphi
      use sfam,     only : phi, phif, fphi, psi
      use linkdt,   only : mmidata   ! 2014_12_16 . scb

!
      include 'global.h'    
      include 'cards.h'
      include 'files.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'defhex.h'
      include 'defpnt.h'
      include 'lscoefh.h'
      include 'deffg.h'
      include 'xsec.h'
      include 'thlink.inc'   ! 2014_09_01 . scb
!
      real,pointer :: psia(:,:), powa(:,:), phia(:,:,:), powu(:), powl(:)
      !real,pointer :: xsfa(:,:
      logical,save :: first=TRUE
      real    :: fr=0.d0, fz=0.d0, fxyz=0.d0   ! 2014_10_06 . scb
      real    :: peaklpd   ! 2014_12_22 . scb
      
      common/edithexvar/psia, powa, phia, powu, powl
      
#ifdef DLL_TASS
! 2014_12_16 . scb
      real :: ao, pow_top, pow_bot
      integer :: kfhalf
      
      type(MMIDATA)  art2mmi2      ! mmi display data
      common / artosmmi / art2mmi2
! added end
#endif      

! 2014_10_06 . scb
      !if(.not. iftran .and. flagout(3)) then
      if(first .and. flagout(3)) then
#ifndef CVF
        write(filename(4), '("out/",a,".peak")') trim(caseid)  ! 2014_10_06 . scb
#else      
        write(filename(4), '(a,".peak")') trim(caseid)   ! 2014_10_06 . scb
#endif      
        open(io4,file=trim(filename(4)),status='unknown')
        write(io4,'(4a12)') 'Time(sec)','Fz','Fr','Fxyz'
      endif
      
      inquire(file=filename(4),exist=filexist)    
      !
      !print *, 'filename(4)', filename(4)
      !print *, 'flagout', flagout(3)
      !print *, 'filexist', filexist
      !
      !pause
      
      if(.not. flagout(3))  filexist=.false.
      
      fr=0.d0
      fz=0.d0
      fxyz=0.d0
! added end

      if(first) then
        call dmalloc(psia,nxya,nza)
        call dmalloc0(powa,0,nxya,0,nza)
        call dmalloc0(phia,1,ng,0,nxya,0,nza)
        call dmalloc0(powu,0,nxya)
        call dmalloc0(powl,0,nxya)
        first=FALSE
      endif
!     
! 2014_07_28 . scb      
      powa=0.d0
      powu=0.d0
      powl=0.d0
! added end

      do ka=1,nza
	    do la=1,nxya
          l=la
          iat=iassytyp(la)
          psia(la,ka)=0
          powa(la,ka)=0
          phia(:,la,ka)=0   ! 2014_10_24 . scb

          psia1=0
          powa1=0
          !phia1=0
          do m2=1,ng2			
            do k=kfmb(ka),kfme(ka)
              volf=volnode(l,k)
              do m=mgb(m2),mge(m2)
                !phia(m,l,k)=phia(m,l,k) + phif(m,l,k)*volf
                phia(m,l,ka)=phia(m,l,ka) + phif(m,l,k)*volf   ! 2014_11_05 . scb
                powa1=powa1 + phif(m,l,k)*volf*xskpf(m,l,k)
              enddo
              psia1=psia1+psi(l,k)	
            enddo			
          enddo !of m2
          phia(:,la,ka)=phia(:,la,ka)/volnodea(la,ka)
          psia(la,ka)=psia1/volnodea(la,ka)
          powa(la,ka)=powa1/volnodea(la,ka)
          if(powa(la,ka) .gt. fxyz)  fxyz=powa(la,ka)    ! 2014_10_06 . scb
        enddo
      enddo
	
      if(ifprint) then   ! 2014_12_22 . scb added ifprint statement
        write(io8,*)
        write(io8,*)
        write(io8,*)
        write(io8,'(a100)') "===========================================  Output Edit   ==========================================="
        write(io8,*)        
				 
! assembly power
        write(io8,*)
        write(io8,*) "Assembly Power Distribution"
        write(io8,*)
      endif      

      powa(:,0)=0.d0   ! 2013_10_14 . scb 
      do ka=kfbega,kfenda
        do la=1,nxya
          iat=iassytyp(la)
          if(.not.iffuela(iat)) cycle

          powa(la,0)=powa(la,0)+powa(la,ka)*volnodea(la,ka)
        enddo
      enddo 

      avgpow=sum(powa(:,0))/volfuel
      do la=1,nxya
        powa(la,0)=powa(la,0)/avgpow/sum(volnodea(la,kfbega:kfenda))
        if(powa(la,0) .gt. fr)  fr=powa(la,0)   ! 2014_10_06 . scb
      enddo

      if(ifprint)  call radformhex(io8,powa(0:nxya,0),8,4,1)   ! 2014_12_22 . scb added ifprint
    
! 2014_12_16 . scb         
#ifdef DLL_TASS
      if(nxya.eq.295) then
        do la=1,nxya
          art2mmi2%powrad(la) = powa(la,0)
        enddo      
      endif      
      
      peaklpd = power/real(nxya)*fr/coreheight
      art2mmi2%peaklpd = int(peaklpd)
      !write(1223,*)  power, nxya, fr, coreheight   ! 2014_12_22 . scb for dbg
      !write(1223,*)  peaklpd, art2mmi2%peaklpd   ! 2014_12_22 . scb for dbg
#endif         
! added end      

!#define VVER1000
!#ifdef VVER1000
!! upper core assembly power
!      write(io8,*)
!      write(io8,*) "Upper Core Assembly Power Distribution"
!      write(io8,*)
!
!      do ka=6,10
!        do la=1,nxya
!          iat=iassytyp(la)
!          if(.not.iffuela(iat)) cycle
!
!          powu(la)=powu(la)+powa(la,ka)*volnodea(la,ka)
!        enddo
!      enddo 
!      do la=1,nxya
!        powu(la)=powu(la)/avgpow/sum(volnodea(la,6:10))
!      enddo
!      call radformhex(io8,powu(0:nxya),8,4,1)
!
!! lower core assembly power
!      write(io8,*)
!      write(io8,*) "Lower Core Assembly Power Distribution"
!      write(io8,*)
!
!      do ka=1,5
!        do la=1,nxya
!          iat=iassytyp(la)
!          if(.not.iffuela(iat)) cycle
!
!          powl(la)=powl(la)+powa(la,ka)*volnodea(la,ka)
!        enddo
!      enddo 
!      do la=1,nxya
!        powl(la)=powl(la)/avgpow/sum(volnodea(la,1:5))
!      enddo
!      call radformhex(io8,powl(0:nxya),8,4,1)
!#endif

! axial power
      if(ifprint) then   ! 2014_12_22 . scb added ifprint statement
        write(io8,*)
        write(io8,*) "Axial Power Distribution"
        write(io8,*)
      endif      
		
      do ka=kfbega,kfenda
        volpl=0
        do la=1,nxya
          iat=iassytyp(la)
          if(.not.iffuela(iat)) cycle

          vol=volnodea(la,ka)
          powa(0,ka)=powa(0,ka)+powa(la,ka)*vol
          volpl=volpl+vol
        enddo
      enddo
      avgpow=sum(powa(0,kfbega:kfenda))/coreheight
      do ka=kfbega,kfenda
        powa(0,ka)=powa(0,ka)/avgpow/hza(ka)
        if(ifprint)  write(io8,'(i7,f10.4)') ka,powa(0,ka)   ! 2014_12_22 . scb added ifprint statement
        if(powa(0,ka) .gt. fz)  fz=powa(0,ka)   ! 2014_10_06 . scb
      enddo
      
! 2014_12_17 . scb         
#ifdef DLL_TASS
      k=0
      do ka=kfbega,kfenda
        k=k+1
        art2mmi2%powax(k) = powa(0,ka)
      enddo      
        
      ao=0.d0
      pow_top = 0.d0
      pow_bot = 0.d0
      
      kfhalf=(kfenda - kfbega + 1)/2 + kfbega -1
            
      do ka=kfbega,kfhalf
        pow_bot = pow_bot + powa(0,ka)*volnodea(1,ka)
      enddo
      
      kfhalf=kfhalf+1
      
      do ka=kfhalf,kfenda
        pow_top = pow_top + powa(0,ka)*volnodea(1,ka)
      enddo
        
      ao=(pow_top - pow_bot)/(pow_top + pow_bot) * 100
      art2mmi2%ao = int(ao)
      
      !write(1223,*)  art2mmi2%ao   ! 2014_12_22 . scb for dbg      
#endif         
! added end          
      
! planar edits
      do ka=kfbega,kfenda
        do la=1,nxya
          vol=volnodea(la,ka)
          powa(la,ka)=powa(la,ka)*vol
        enddo
      enddo

      avgpow=sum(powa(1:nxya,1:nza))/volfuel
      do ka=kfbega,kfenda
        do la=1,nxya
          vol=volnodea(la,ka)
          powa(la,ka)=powa(la,ka)/avgpow/vol
        enddo
      enddo

      if(ifprint) then   ! 2014_12_22 . scb added ifprint statement
        do ka=kfbega,kfenda
          write(io8,*)
          write(io8,*) "Planar Power Distribution at Plane",ka
          write(io8,*)
          call radformhex(io8,powa(0:nxya,ka),8,4,1)
        enddo
      endif      

! 2014_09_01 . scb commented 
      !do m=1,2
      !  write(io8,*)
      !  write(io8,*) "Group Flux Distribution at Plane",m
      !  write(io8,*)
      !  call radformhex(io8,phif(m,0:nxya,1),8,4,1)
      !enddo     

! planar edits 2
      if(ifprint)  call radformhex2(io8,powa(1:nxya,1:nza),8,4,1,fxyz)   ! 2014_12_22 . scb added ifprint statement
      ! 2014_10_06 . scb add fxyz

      if(filexist .and. ifprint)  write(io4,'(4f12.4)')  time,fz,fr,fxyz    ! 2014_10_06 . scb                
      
! 2014_10_24 . scb
      if(ifprint .and. flagout(6)) call radformhex3(1,io8,phia(1:ng,1:nxya,1:nza),8,4,1)  ! 2014_12_22 . scb added ifprint statement
      if(nz.ne.nza) then
        if(flagout(7) .or. flagout(8)) then
          print *, "Each axial node should not be split to print fission XS and nu !"
        endif
      !else      
      elseif(ifprint) then   ! 2014_12_22 . scb
        if(flagout(7)) call radformhex3(2,io8,xsff(1:ng,1:nxya,1:nza),8,4,1)
        if(flagout(6)) call radformhex3(3,io8,xsnuf(1:ng,1:nxya,1:nza),8,4,1)
      endif
! added end      
      
      

      return		
      
      
      
      
      
      
      
      
      
      
! planar edits
      !write(io8,*) "Flux Distribution"
      !do ka=kfbega,kfenda
      !  write(io8,*)
      !  write(io8,*) "Planar Flux Distribution at Plane",ka
      !  write(io8,*)
      !  call radformhex_flux(io8,phia(0:nxya,ka),10,6,1)
      !enddo

      
      do ka=1,nza
        write(io8,*)
        write(io8,*) "Planar Flux Distribution at Plane",ka
        do m=1,ng
          write(io8,*)
          write(io8,*) "Group",m
          write(io8,*) 
          call radformhex_flux(io8,phif(m,0:nxya,ka),10,6,1)
        enddo
      enddo


! planar edits 2
      !call radformhex2_flux(io8,phia(1:nxya,1:nza),8,0,1)

      return
    end subroutine
!=======================================================================================!
    subroutine radformhex(idev,assdat,nlong,nprec,ife)
!
      use param
      use allocs
!
      include 'global.h'      
      include 'cards.h'
      include 'files.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'defhex.h'
      include 'defpnt.h'
      include 'lscoefh.h'
      include 'deffg.h'
!
      character*66 if1,if2,if3
      dimension assdat(0:nxya)
      character*256 form
      integer :: mx, my  ! 2012_10_15 . scb
      
      if(ife.eq.1) then
!         write(form,'("(5x,100i",i2,")")') nlong/2
         write(form,'("(12x,100i",i2,")")') nlong/2  ! 2012_10_15 . scb
      else
         write(form,'("(3x,100i",i2,")")') nlong/2
      endif
      write(idev,form) (ia,ia=icoreff,icordxe,2)
      write(idev,*)
      do iy=icoreyfs,icoreyfe,2
         nblnk=icorexfs(iy)-icoreff
         nblnk=nblnk*nlong/4+1
         if(ife.eq.1) then
           write(form,'("(i7,",i3,"x,100f",i2,".",i2,")")') nblnk,nlong,nprec
         else
           write(form,'("(i7,",i3,"x,1p,100e",i2,".",i2,")")') nblnk,nlong,nprec
         endif
         write(idev,form) iy,(assdat(ltola(iaass(ix,iy))),ix=icorexfs(iy),icorexfe(iy),4)
      enddo
      
! 2012_10_15 . scb
      peak=0.
      do iy=icoreyfs,icoreyfe,2
        do ix=icorexfs(iy),icorexfe(iy),4
          if(assdat(ltola(iaass(ix,iy))) .gt. peak) then
            mx=ix
            my=iy
            peak=assdat(ltola(iaass(ix,iy)))
          endif
        enddo
      enddo
! added end      
      
      write(idev,*)
      write(idev,*)"       Maximum Pos. Maximum Value "
      if(ife.eq.1) then
         write(form,'("(9x,a,i3,a,i3,a,f",i2,".",i2,")")') nlong,nprec
      else
         write(form,'("(9x,a,i3,a,i3,a,1p,e",i2,".",i2,")")') nlong,nprec
      endif
! 2012_10_15 . scb      
!      write(idev,form) '(',i1,',',i2,' ) ',peak
      write(idev,form) '(',mx,',',my,' ) ',peak
! added end      
      write(idev,*)
      return
    end subroutine
!=======================================================================================!
    subroutine radformhex_flux(idev,assdat,nlong,nprec,ife)
!
      use param
      use allocs
!
      include 'global.h'  
      include 'cards.h'
      include 'files.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'defhex.h'
      include 'defpnt.h'
      include 'lscoefh.h'
      include 'deffg.h'
!
      character*66 if1,if2,if3
      dimension assdat(0:nxya)
      character*256 form
      if(ife.eq.1) then
!         write(form,'("(5x,100i",i2,")")') nlong/2
         write(form,'("(12x,100i",i2,")")') nlong/2  ! 2012_10_15 . scb
      else
         write(form,'("(3x,100i",i2,")")') nlong/2
      endif
      write(idev,form) (ia,ia=icoreff,icordxe,2)
      write(idev,*)
      do iy=icoreys,icoreye,2
         nblnk=icorexs(iy)-icoreff
         nblnk=nblnk*nlong/4+1
         if(ife.eq.1) then
           write(form,'("(i7,",i3,"x,100f",i2,".",i2,")")') nblnk,nlong,nprec
         else
           write(form,'("(i7,",i3,"x,1p,100e",i2,".",i2,")")') nblnk,nlong,nprec
         endif
         write(idev,form) iy,(assdat(ltola(iaass(ix,iy))), ix=icorexs(iy),icorexe(iy),4)
      enddo
      
      write(idev,*)
      write(idev,*)"       Maximum Pos. Maximum Value "
      if(ife.eq.1) then
         write(form,'("(9x,a,i3,a,i3,a,f",i2,".",i2,")")') nlong,nprec
      else
         write(form,'("(9x,a,i3,a,i3,a,1p,e",i2,".",i2,")")') nlong,nprec
      endif
      write(idev,form) '(',i1,',',i2,' ) ',peak
      write(idev,*)
      return
    end subroutine
!=======================================================================================!
    subroutine radformhex2(idev,assdat,nlong,nprec,ife,peak)
!
      use param
      use allocs
      include 'global.h'    
      include 'cards.h'
      include 'files.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'defhex.h'
      include 'defpnt.h'
      include 'lscoefh.h'
      include 'deffg.h'
!
      dimension assdat(1:nxya,1:nza)
      character*256 form
      logical :: first
      real    :: peak  ! 2014_10_06 . scb
      integer :: iapeak, kapeak

      write(idev,*) 
      write(idev,*) "Region Index"
      write(idev,*) 
!
      peak=0.d0  ! 2014_07_28 . scb
      
      if(ife.eq.1) then
!         write(form,'("(5x,100i",i2,")")') nlong/2
         write(form,'("(12x,100i",i2,")")') nlong/2  ! 2012_10_15 . scb
      else
        write(form,'("(3x,100i",i2,")")') nlong/2
      endif
      write(idev,form) (ia,ia=icoreff,icordxe,2)
      write(idev,*)
      ia=0
      itemp=0
      do iy=icoreyfs,icoreyfe,2
        nblnk=icorexfs(iy)-icoreff
        nblnk=nblnk*nlong/4+1
        first=TRUE
        write(form,'("(i7,",i3,"x,100i",i2,".",i2,")")') nblnk,nlong,nprec
        ia=ia+itemp+1
        itemp=(icorexfe(iy)-icorexfs(iy))/4
        write(idev,form) iy,(ix,ix=ia,ia+itemp)
      enddo
      write(idev,*)
!

      write(idev,*) 
      write(idev,*) "Power Distribution in All Region"
      write(idev,*) 
            
      
      write(idev,'(a7,50a15)') '      ', 'Bottom', ('',ka=kfbega+1,kfenda-1), 'Top'
      write(idev,'(a7,50i15)') 'Assem.', (ka,ka=kfbega,kfenda)
!
      ia=0
      do iy=icoreyfs,icoreyfe,2
        do ix=icorexfs(iy),icorexfe(iy),4
          ia=ia+1
          write(idev,'(i7,50(1pe15.6))') ia,(assdat(ltola(iaass(ix,iy)),ka),ka=kfbega,kfenda)
          !if(ia.eq.180) print *, ix, iy, iaass(ix,iy), ltola(iaass(ix,iy))  ! 2015_02_11 . scb for SC
! 2014_04_16 . scb          
          do ka=kfbega,kfenda
            if(assdat(ltola(iaass(ix,iy)),ka) .gt. peak) then
              peak=assdat(ltola(iaass(ix,iy)),ka)
              iapeak=ia
              kapeak=ka   ! 2014_07_25 . scb fixed bug
            endif
          enddo
! added end          
        enddo
      enddo
      
! 2014_04_21 . scb      
      write(idev,*)
      write(idev,*)"       Maximum Pos. Maximum Value "
      
      if(ife.eq.1) then
         write(form,'("(9x,a,i3,a,i3,a,f",i2,".",i2,")")') nlong,nprec
      else
         write(form,'("(9x,a,i3,a,i3,a,1p,e",i2,".",i2,")")') nlong,nprec
      endif
      
      write(idev,form) '(',iapeak,',',kapeak,' ) ',peak
      write(idev,*)      
! added end
!
!! 2014_10_24 . scb
!      if(flagout(6)) then
!        write(idev,*) 
!        write(idev,*) "Flux Distribution in All Region"
!        write(idev,*) 
!      
!        do m=1,ng
!          write(idev,*) "Group : ",m
!          write(idev,*) 
!      
!          write(idev,'(a7,50a15)') '      ', 'Bottom', ('',ka=kfbega+1,kfenda-1), 'Top'
!          write(idev,'(a7,50i15)') 'Assem.', (ka,ka=kfbega,kfenda)
!    !
!          ia=0
!          do iy=icoreyfs,icoreyfe,2
!            do ix=icorexfs(iy),icorexfe(iy),4
!              ia=ia+1
!              write(idev,'(i7,50(1pe15.6))') ia,(assdat(ltola(iaass(ix,iy)),ka),ka=kfbega,kfenda)  
!            enddo
!          enddo
!          write(idev,*) 
!        enddo
!      endif
!      
      

      return
    end subroutine
!=======================================================================================!

    subroutine radformhex2_flux(idev,assdat,nlong,nprec,ife)
!
      use param
      use allocs
      include 'global.h'   
      include 'cards.h'
      include 'files.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'defhex.h'
      include 'defpnt.h'
      include 'lscoefh.h'
      include 'deffg.h'
!
      dimension assdat(1:nxya,1:nza)
      character*256 form
      logical :: first

      write(idev,*) 
      write(idev,*) "Flux Distribution in All Region"
      write(idev,*) 
!
      if(ife.eq.1) then
!         write(form,'("(5x,100i",i2,")")') nlong/2
         write(form,'("(12x,100i",i2,")")') nlong/2  ! 2012_10_15 . scb
      else
        write(form,'("(3x,100i",i2,")")') nlong/2
      endif
      write(idev,form) (ia,ia=icoreff,icordxe,2)
      write(idev,*)
      ia=0
      itemp=0
      do iy=icoreys,icoreye,2
        nblnk=icorexs(iy)-icoreff
        nblnk=nblnk*nlong/4+1
        first=TRUE
        write(form,'("(i7,",i3,"x,100i",i2,".",i2,")")'),nblnk,nlong,nprec
        ia=ia+itemp+1
        itemp=(icorexe(iy)-icorexs(iy))/4
        write(idev,form) iy,(ix,ix=ia,ia+itemp)
      enddo
      write(idev,*)
!
      write(idev,'(a7,50a15)') '      ', 'Bottom',('',ka=kfbega+1,kfenda-1), 'Top'
      write(idev,'(a7,50i15)') 'Assem.', (ka,ka=kfbega,kfenda)
!
      ia=0
      do iy=icoreys,icoreye,2
        do ix=icorexs(iy),icorexe(iy),4
          ia=ia+1
          write(idev,'(i7,50(1pe15.6))') ia,(assdat(ltola(iaass(ix,iy)),ka),ka=kfbega,kfenda)
        enddo
      enddo

      return
    end subroutine
    
!=======================================================================================!
! 2014_10_24 . scb    
    subroutine radformhex3(itype,idev,assdat,nlong,nprec,ife)
!
      use param
      use allocs
      include 'global.h'    
      include 'cards.h'
      include 'files.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'defhex.h'
      include 'defpnt.h'
      include 'lscoefh.h'
      include 'deffg.h'
!
      dimension assdat(1:ng,1:nxya,1:nza)
      character*256 form
      logical :: first
      integer :: itype
      
      write(idev,*) 
      if(itype.eq.1) write(idev,*) "Flux Distribution in All Region"
      if(itype.eq.2) write(idev,*) "Fission XS in All Region"
      if(itype.eq.3) write(idev,*) "Nu value in All Region"
      write(idev,*) 
      
      do m=1,ng
        write(idev,*) "Group : ",m
        write(idev,*) 
      
        write(idev,'(a7,50a15)') '      ', 'Bottom', ('',ka=kfbega+1,kfenda-1), 'Top'
        write(idev,'(a7,50i15)') 'Assem.', (ka,ka=kfbega,kfenda)
  !
        ia=0
        do iy=icoreyfs,icoreyfe,2
          do ix=icorexfs(iy),icorexfe(iy),4
            ia=ia+1
            write(idev,'(i7,50(1pe15.6))') ia,(assdat(m,ltola(iaass(ix,iy)),ka),ka=kfbega,kfenda)  
          enddo
        enddo
        write(idev,*) 
      enddo
      
      

      return
    end subroutine
! added end    
!=======================================================================================!    