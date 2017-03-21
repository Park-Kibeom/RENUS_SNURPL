    subroutine homopin
    
      use allocs
      use param

      include 'global.h'
      include 'geom.h'
      include 'xsec.h'
      include 'nodal.h'
      include 'pinpr.h'
      include 'ff.h'
      include 'ffdm.h'
      include 'files.h'
!
      real, parameter                     :: wxn=2.0, wyn=2.0
      real :: wxa,wya, &         ! normalized total Width of Assembly: x- and y- 
              dxpin,dypin, &     ! wxa/npin, wya/npin
              lnx,lny, &         ! coordinate of the  target  node in one assembly
              txu,tyu, &         ! temporary x-upper, temporary y-upper
              epspin, &          ! 1.0E-5
              itnx,itny, &       ! coordinate of the related node 
              exx,exy, &         ! value 
              xl(2),xu(2),yl(2),yu(2), &
              phit(4),tarea(4)
!!FIXME : undefine MATLAB_DEBUG
!!define MATLAB_DEBUG
!#ifdef MATLAB_DEBUG
!     +,           phih2(ng2,npin,npin,nxya,nz)
!     +,           powh(npin,npin,nxya,nz)
!#endif
!
      epspin=1.0E-5
      do m=1,ng
        do k=1,nz
          ka=ktoka(k)

          do la=1,nxya
            ia=latoia(la)
            ja=latoja(la)
            wxa=nmeshx(ia)*wxn
            wya=nmeshy(ja)*wyn 
            dxpin=wxa/npin
            dypin=wya/npin
            yl(2)=-one ! y-lower
            lny=1    ! coordinate 
            do ipy=1,npin   ! pin y-directional              
              yl(1)=yl(2)
              tyu=yl(2)+dypin

              xl(2)=-one
              lnx=1
              do ipx=1,npin
                xl(1)=xl(2)
                txu=xl(2)+dxpin
                if((tyu-1.0)<=epspin) then
                  itny1=lny
                  yu(1)=tyu
                  if((txu-1.0)<=epspin) then
                    itnx1=lnx
                    xu(1)=txu
                    tphi=phihav(xl(1),xu(1),yl(1),yu(1),latol(la)%ij(itnx1,itny1),k,m)
                    phih(ipx,ipy,la,k,m)=tphi
                    xl(2)=xu(1)
                  elseif(txu.gt.1.0) then
                    itnx1=lnx
                    itnx2=lnx+1
                    xu(1)=1
                    exx=txu-xu(1) !extra x
                    xl(2)=-one
                    xu(2)=-one+exx
                    tsum=0
                    do iw=1,2
                      tarea(iw)=xu(iw)-xl(iw)
                      tsum=tsum+tarea(iw)
                    enddo
                    phit=0
                    if(tarea(1) .ne. 0) phit(1)=phihav(xl(1),xu(1),yl(1),yu(1),latol(la)%ij(itnx1,itny1),k,m)
                    if(tarea(2) .ne. 0) phit(2)=phihav(xl(2),xu(2),yl(1),yu(1),latol(la)%ij(itnx2,itny1),k,m)                    
                    trsum=1/tsum
                    tphi=0
                    do iw=1,2
                      tphi=tphi+tarea(iw)*phit(iw)
                    enddo
                    tphi=tphi*trsum
                    phih(ipx,ipy,la,k,m)=tphi
                    xl(2)=xu(2)
                    lnx=lnx+1
                  endif ! of txu
                  yl(2)=yu(1)
                elseif(tyu.gt.1) then
                  itny1=lny
                  itny2=lny+1
                  yu(1)=1
                  exy=tyu-yu(1) !extra y
                  yl(2)=-one
                  yu(2)=-one+exy
                  if((txu-1.0)<=epspin) then
                    itnx1=lnx
                    xu(1)=txu
                    tsum=0
                    do iw=1,2
                      tarea(iw)=yu(iw)-yl(iw)
                      tsum=tsum+tarea(iw)
                    enddo
                    phit=0
                    if(tarea(1) .ne. 0) phit(1)=phihav(xl(1),xu(1),yl(1),yu(1),latol(la)%ij(itnx1,itny1),k,m)
                    if(tarea(2) .ne. 0) phit(2)=phihav(xl(1),xu(1),yl(2),yu(2),latol(la)%ij(itnx1,itny2),k,m)
                    
                    trsum=1/tsum
                    tphi=0
                    do iw=1,2
                      tphi=tphi+tarea(iw)*phit(iw)
                    enddo
                    tphi=tphi*trsum
                    phih(ipx,ipy,la,k,m)=tphi
                    xl(2)=xu(1)
                  elseif(txu.gt.1) then
                    itnx1=lnx
                    itnx2=lnx+1
                    xu(1)=1
                    exx=txu-xu(1) !extra x
                    xl(2)=-one
                    xu(2)=-one+exx
                    tarea(1)=(xu(1)-xl(1))*(yu(1)-yl(1))  
                    tarea(2)=(xu(2)-xl(2))*(yu(1)-yl(1))  
                    tarea(3)=(xu(1)-xl(1))*(yu(2)-yl(2))  
                    tarea(4)=(xu(2)-xl(2))*(yu(2)-yl(2))
                    tsum=0
                    do iw=1,4
                      tsum=tsum+tarea(iw)
                    enddo
                    phit=0
                    if(tarea(1) .ne. 0) phit(1)=phihav(xl(1),xu(1),yl(1),yu(1),latol(la)%ij(itnx1,itny1),k,m)
                    if(tarea(2) .ne. 0) phit(2)=phihav(xl(2),xu(2),yl(1),yu(1),latol(la)%ij(itnx2,itny1),k,m)
                    if(tarea(3) .ne. 0) phit(3)=phihav(xl(1),xu(1),yl(2),yu(2),latol(la)%ij(itnx1,itny2),k,m)
                    if(tarea(4) .ne. 0) phit(4)=phihav(xl(2),xu(2),yl(2),yu(2),latol(la)%ij(itnx2,itny2),k,m)

                    trsum=1/tsum
                    tphi=0
                    do iw=1,4
                      tphi=tphi+tarea(iw)*phit(iw)
                    enddo
                    tphi=tphi*trsum
                    phih(ipx,ipy,la,k,m)=tphi
                    xl(2)=xu(2)
                    lnx=lnx+1
                  endif ! of txu
                  if(ipx.eq.npin) lny=lny+1
                  yl(2)=yu(2)
                endif ! of tyu
              enddo ! of x pin
            enddo ! of y pin
          enddo ! of la
        enddo ! of k
      enddo ! of ng
!
!
!!FIXME : undefine HOMOFLUX
!#ifdef MATLAB_DEBUG
!!collapse
!!      do k=1, nz
!!        do la=1, nxya
!!          do ix=1,npin
!!            do iy=1,npin
!!              do m2=1, ng2
!!                phih2(m2,ix,iy,la,k) = 0
!!                do m=mgb(m2), mge(m2)
!!                  phih2(m2,ix,iy,la,k) = phih2(m2,ix,iy,la,k) + phih(ix,iy,la,k,m)
!!                enddo
!!              enddo !m2
!!            enddo !iy
!!          enddo !ix
!!        enddo !la
!!      enddo !k
!
!! pincell wise flux and power
!      do k=kfmb(ka),kfme(ka)
!        do la=1,nxya
!          do ix=1,npin
!            do iy=1,npin
!              powh1=0            
!              do m=1,ng
!                powh1=powh1+phih(ix,iy,la,k,m)*xskpf(m,latol(la)%fm(1),k)
!              enddo !m
!              powh(ix,iy,la,k)=powh1
!            enddo !iy
!          enddo !ix
!        enddo !la
!      enddo !k     
!      
!      iodbg=2001
!      write(filename(1), '("debug/renus_power_",i2.2,"box.out")') (nx/nxa)**2
!      open(iodbg, file=trim(filename(1)),status='unknown')
!      do k=1,nz
!        do ja=1,nya
!          do iy=1,npin
!            do ia=1,nxa
!              la=nodela(ia,ja)
!              do ix=1,npin
!                write(iodbg,600) powh(ix,iy,la,k)
!              enddo
!            enddo !ia
!            write(iodbg,'("")')
!          enddo !iy
!        enddo !ja
!      enddo !k
!      close(iodbg)
!      
!      do m=1,ng
!        write(filename(1), '("debug/renus_flux_G",i1,"_",i2.2,"box.out")') m, (nx/nxa)**2
!        open(iodbg, file=trim(filename(1)),status='unknown')
!        do k=1,nz
!          do ja=1,nya
!            do iy=1,npin
!              do ia=1,nxa
!                la=nodela(ia,ja)
!                do ix=1,npin
!                  write(iodbg,600) phih(ix,iy,la,k,m)
!                enddo
!              enddo !ia
!              write(iodbg,'("")')
!            enddo !iy
!          enddo !ja
!        enddo !k
!        close(iodbg)
!      enddo
! 600  format(1x,1p,e12.5,$)
!#endif
!
      return     
      
    end subroutine
