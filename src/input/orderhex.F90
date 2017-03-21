! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
    subroutine orderhex
!
! generate ordering of hex nodes, surfaces, and points, and neighbor info.

      use param
      use allocs
!
      include 'global.h'
      include 'itrcntl.h'
      include 'files.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'xsec.h'
      include 'defhex.h'
      include 'defsfc.h'
      include 'defpnt.h'
!
      dimension value(100)
      integer mp(2)
      data mp/2,1/
!=====================================================================
!                          hex assembly ordering
!=====================================================================
      if(isolang.lt.360) then
        do iy=icordys,-1
          do ix=icordxs,icordxe
            iaass(ix,iy)=0
          enddo
        enddo
        do ix=icordxs,-1
          iaass(ix,0)=0
        enddo
      endif
      if(isolang.eq.120) then
        do iy=1,icordye
          do ix=icordxs,-iy
            iaass(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.60) then
        do iy=1,icordye
          ix2=iy-1
          if(isymtype.eq.1) ix2=iy
          do ix=icordxs,ix2
            iaass(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.30) then
        do iy=1,icordye
          ix2=iy*3-1
          if(ix2.gt.nxfc) ix2=nxfc
          do ix=icordxs,ix2
            iaass(ix,iy)=0
          enddo
        enddo
      endif
!
! origins in natural number assembly coordinate
      if(isolang.eq.360) then
         i0=2*nring+1
         j0=nring+1
      elseif(isolang.eq.120) then
         i0=nring+1
         j0=1
      else
         i0=1
         j0=1
      endif
!
! consider missing-vertice case
      if(i0.ne.1 .and. icordxe/4.ne.nring) i0=i0-1
!
      call dmalloc0(ltoi,0,nxy)
      call dmalloc0(ltoj,0,nxy)
      call dmalloc(nrnx,ny+1)
      call dmalloc0(iaass2d,-3,nx+4,-1,ny+2)   ! 2014_08_23 . scb
      call dmalloc(iaass1d,nxy)                ! 2014_08_23 . scb
      call dmalloc(iaass1dinv,nxy)             ! 2014_08_23 . scb

      inum=0
      imvalx=mod(icordxe,4)
      do icy=1,3
        do iy=icordys,icordye,2
          imval=mod(iy,4)
          if(imvalx.eq.0) then
            if(imval.eq.0) then
              ittt=icordxs
            else
              ittt=icordxs-6
            endif
          else
            if(imval.eq.0) then
              ittt=icordxs-6
            else
              ittt=icordxs
            endif
          endif
          ittt=ittt+(icy-1)*4
!
! y-location on natural number assembly coordinate
          j=iy/2+j0
!
          do ix=ittt,icordxe,12
            if(iaass(ix,iy).ne.0) then
              inum=inum+1
! x-location on natural number assembly coordinate
              i=ix/2+i0
              l=inum
              nodel(i,j)=l
              ltoi(l)=i
              !ltoi(l)=j
              ltoj(l)=j  ! 2014_04_08 . scb
              
              if(iaass(ix-4,iy).eq.0) nxs(j)=i
              if(iaass(ix+4,iy).eq.0) nxe(j)=i
!
              iassytyp(inum)=iaass(ix,iy)  
              iaass(ix,iy)=inum
              iaass2d(i,j)=inum   ! 2014_08_23 . scb
              
              if(isolang.eq.30) then
                if(iy.eq.0) wtass(inum)=0.5
                if(ix.eq.(3*iy)) wtass(inum)=0.5
              endif
              if(isolang.eq.60.and.isymtype.eq.2) then
                if(iy.eq.0) wtass(inum)=0.5
                if(ix.eq.iy) wtass(inum)=0.5
              endif
            endif
          enddo
        enddo 
      enddo 
      
! 2014_08_23 . scb
      inum=0
      do iy=-1,ny+2
        do ix=-3,nx+4
          if(iaass2d(ix,iy).ne.0) then
            inum=inum+1
            iaass1d(inum)=iaass2d(ix,iy)
            iaass1dinv(iaass2d(ix,iy))=inum
          endif          
        enddo 
      enddo
! added end

      do j=1,ny  !hj 10mar04
        nrnx(j)=nxe(j)-nxs(j)+1
      enddo
!
      if(isolang.eq.60) wtass(iaass(0,0))=1./6.d0
      if(isolang.eq.30) wtass(iaass(0,0))=1./12.d0
      if(isolang.eq.120) wtass(iaass(0,0))=1./3.d0
      nassy=inum
!=====================================================================
!                           surface ordering
!=====================================================================
      if(isolang.ne.360) then
        do iy=icordys-1,-1
          do ix=icordxs-2,icordxe+2
            iasurf(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.120) then
        do iy=0,icordye+1
          do ix=icordxs-2,-iy
            iasurf(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.60) then
        do iy=0,icordye+1
          ix2=iy-1
          if(isymtype.eq.1) ix2=iy
          do ix=icordxs-2,ix2
            iasurf(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.30) then
        do iy=0,icordye+1
          ix2=iy*3-2
          if(ix2.gt.nxfc) ix2=nxfc
          do ix=icordxs-2,ix2
            iasurf(ix,iy)=0
          enddo
        enddo
      endif
!
      inum=0 
      do iy=icordys-1,icordye+1 
        do ix=icordxs-2,icordxe+2 
          item=iasurf(ix,iy)
          if(item.ne.0) then
            inum=inum+1 
            iasurf(ix,iy)=inum 
          endif
        enddo
      enddo
      nxsfc=inum
!=====================================================================
!                           point ordering
!=====================================================================
      if(isolang.ne.360) then
        do iy=icordys-1,-1
          do ix=icordxs-2,icordxe+2
            iapoint(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.120) then
        do iy=0,icordye+1
          do ix=icordxs-2,-iy-1
            iapoint(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.60) then
        do iy=0,icordye+1
          do ix=icordxs-2,iy-1
            iapoint(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.30) then
        do iy=0,icordye+1
          ix2=iy*3-2
          if(ix2.gt.nxfc) ix2=nxfc
          do ix=icordxs-2,ix2
            iapoint(ix,iy)=0
          enddo
        enddo
      endif
! Set Value at Boundary Point and Point numbering 
      do m=1,ng 
         value(m)=4964.d0*rt3*hside*albr(m)/74165.d0
      enddo
      if(hside.lt.5.0.or.ndivhs.ge.2) then
         fac=1.
         fac2=1.
      else
         fac=1.-rt3/2.
         fac2=1.
      endif
      inum=0 
      do iy=icordys-1,icordye+1 
        do ix=icordxs-2,icordxe+2 
          nbdy=iapoint(ix,iy)
          if(nbdy.ne.0) then
            inum=inum+1 
            codpnt(inum)=0
            if(nbdy.eq.1) then
              do m=1,ng
                 pbdv(m,inum)=value(m)*fac2
              enddo
              codpnt(inum)=1
            endif
            if(nbdy.eq.2) then
              do m=1,ng
                 pbdv(m,inum)=value(m)*fac
              enddo
              codpnt(inum)=0.5
            endif
            iapoint(ix,iy)=inum 
          endif
        enddo
      enddo
      nxpnt=inum
!=====================================================================
!                     generate neighbor information
!=====================================================================
      do iy=-nyfc,nyfc
        do ix=-nxfc,nxfc
          layh(ix,iy)=0
          layp(ix,iy)=0
          lays(ix,iy)=0
        enddo
      enddo
      do iy=icordys-1,icordye+1
        do ix=icordxs-2,icordxe+2
          layh(ix,iy)=iaass(ix,iy)
          layp(ix,iy)=iapoint(ix,iy)
          lays(ix,iy)=iasurf(ix,iy)
        enddo
      enddo
#ifdef DBG
      write(111,*) 'ipoint'
      do iy=icordys-1,icordye+1,2
        write(111,'(100i3)') (layp(ix,iy),ix=icordxs-2,icordxe+2,2)
      enddo
      write(111,*)  
#endif
      if(isolang.eq.30) then
        do iy=0,icordye,2
          do ix=3*iy,icordxe,4
            layh(ix/2+3*iy/2,ix/2-iy/2)=layh(ix,iy)
          enddo
        enddo
        do iy=0,icordye,2
          do ix=3*iy+2,icordxe+2,4
            lays(ix/2+3*iy/2,ix/2-iy/2)=lays(ix,iy)
          enddo
        enddo
        do iy=1,icordye+1,2
          do ix=3*iy,icordxe+2,4
           lays((ix-3)/2+3*(iy-1)/2+3,(ix-3)/2-(iy-1)/2+1)=lays(ix,iy)
          enddo
          do ix=3*iy+2,icordxe+2,4
           lays((ix-5)/2+3*(iy-1)/2+4,(ix-5)/2-(iy-1)/2+2)=lays(ix,iy)
          enddo
        enddo
        do iy=1,icordye+1,2
          do ix=3*iy-1,icordxe+2,4
           layp((ix-2)/2+3*(iy-1)/2+2,(ix-2)/2-(iy-1)/2+1)=layp(ix,iy)
          enddo
          do ix=3*iy+1,icordxe+2,4
           layp((ix-4)/2+3*(iy-1)/2+4,(ix-4)/2-(iy-1)/2+1)=layp(ix,iy)
          enddo
        enddo
      endif
#ifdef DBG
      write(111,*) 'ipoint'
      do iy=icordys-1,icordye+1,2
        write(111,'(100i3)') (layp(ix,iy),ix=icordxs-2,icordxe+2,2)
      enddo
      write(111,*)  
#endif
      if(isolang.le.60) then
        do iy=0,icordye,2
          do ix=iy,icordxe,4
            iyto=ix/2+iy/2
            if(isymtype.eq.1) then
              ixto=ix/2-3*iy/2
            else
              ixto=-ix/2+3*iy/2
            endif
            if(indbetwn(ixto,iyto).eq.1) layh(ixto,iyto)=layh(ix,iy)
          enddo
        enddo
        do iy=0,icordye,2
          do ix=iy+2,icordxe+2,4
            iyto=(ix-2)/2+iy/2+1
            if(isymtype.eq.1) then
              ixto=(ix-2)/2-3*iy/2+1
            else
              ixto=-(ix-2)/2+3*iy/2-1
            endif
            if(indbetwn(ixto,iyto).eq.1) lays(ixto,iyto)=lays(ix,iy)
          enddo
        enddo
        do iy=1,icordye+1,2
          do ix=iy+2,icordxe+2,4
            iyto=(ix-3)/2+(iy-1)/2+2
            if(isymtype.eq.1) then
              ixto=(ix-3)/2-3*(iy-1)/2
            else
              ixto=-(ix-3)/2+3*(iy-1)/2
            endif
            if(indbetwn(ixto,iyto).eq.1) lays(ixto,iyto)=lays(ix,iy)
          enddo
          do ix=iy+4,icordxe+2,4
            iyto=(ix-5)/2+(iy-1)/2+3
            if(isymtype.eq.1) then
              ixto=(ix-5)/2-3*(iy-1)/2+1
            else
              ixto=-(ix-5)/2+3*(iy-1)/2-1
            endif
            if(indbetwn(ixto,iyto).eq.1) lays(ixto,iyto)=lays(ix,iy)
          enddo
        enddo
        do iy=1,icordye+1,2
          do ix=iy+1,icordxe+2,4
            iyto=ix/2+(iy-1)/2
            if(isymtype.eq.1) then
              ixto=ix/2-3*(iy-1)/2-1
            else
              ixto=-ix/2+3*(iy-1)/2+1
            endif
            if(indbetwn(ixto,iyto).eq.1) layp(ixto,iyto)=layp(ix,iy)
          enddo
          do ix=iy+3,icordxe+2,4
            iyto=ix/2+(iy-1)/2+1
            if(isymtype.eq.1) then
              ixto=ix/2-3*(iy-1)/2-2
            else
              ixto=-ix/2+3*(iy-1)/2+2
            endif
            if(indbetwn(ixto,iyto).eq.1) layp(ixto,iyto)=layp(ix,iy)
          enddo
        enddo
      endif
#ifdef DBG
      write(111,*) 'ipoint'
      do iy=icordys-1,icordye+1,2
        write(111,'(100i3)') (layp(ix,iy),ix=icordxs-2,icordxe+2,2)
      enddo
      write(111,*)
#endif  
      if(isolang.le.120) then
        do iy=0,icordye,2
          do ix=-iy+4,icordxe,4
            ixto=-ix/2-3*iy/2
            iyto=ix/2-iy/2
            if(indbetwn(ixto,iyto).eq.1) layh(ixto,iyto)=layh(ix,iy)
            ixto=-ix/2+3*iy/2
            iyto=-ix/2-iy/2
            if(indbetwn(ixto,iyto).eq.1) layh(ixto,iyto)=layh(ix,iy)
          enddo
        enddo
        do iy=0,icordye+1,2
          do ix=-iy+2,icordxe+2,4
            ixto=-(ix-2)/2-3*iy/2-1
            iyto=(ix-2)/2-iy/2+1
            if(indbetwn(ixto,iyto).eq.1) lays(ixto,iyto)=lays(ix,iy)
            ixto=-(ix-2)/2+3*iy/2-1
            iyto=-(ix-2)/2-iy/2-1
            if(indbetwn(ixto,iyto).eq.1) lays(ixto,iyto)=lays(ix,iy)
          enddo
        enddo
        do iy=1,icordye+1,2
          do ix=-iy+2,icordxe+2,4
            ixto=-(ix-1)/2-3*(iy-1)/2-2
            iyto=(ix-1)/2-(iy-1)/2
            if(indbetwn(ixto,iyto).eq.1) lays(ixto,iyto)=lays(ix,iy)
            ixto=-(ix-1)/2+3*(iy-1)/2+1
            iyto=-(ix-1)/2-(iy-1)/2-1
            if(indbetwn(ixto,iyto).eq.1) lays(ixto,iyto)=lays(ix,iy)
          enddo
          do ix=-iy+4,icordxe+2,4
            ixto=-(ix-3)/2-3*(iy-1)/2-3
            iyto=(ix-3)/2-(iy-1)/2+1
            if(indbetwn(ixto,iyto).eq.1) lays(ixto,iyto)=lays(ix,iy)
            ixto=-(ix-3)/2+3*(iy-1)/2
            iyto=-(ix-3)/2-(iy-1)/2-2
            if(indbetwn(ixto,iyto).eq.1) lays(ixto,iyto)=lays(ix,iy)
          enddo
        enddo
        do iy=1,icordye+1,2
          do ix=-iy+3,icordxe+2,4
            ixto=-ix/2-3*(iy-1)/2-1
            iyto=ix/2-(iy-1)/2
            if(indbetwn(ixto,iyto).eq.1) layp(ixto,iyto)=layp(ix,iy)
            ixto=-ix/2+3*(iy-1)/2+1
            iyto=-ix/2-(iy-1)/2
            if(indbetwn(ixto,iyto).eq.1) layp(ixto,iyto)=layp(ix,iy)
          enddo
          do ix=-iy+1,icordxe+2,4
            ixto=-ix/2-3*(iy-1)/2-2
            iyto=ix/2-(iy-1)/2-1
            if(indbetwn(ixto,iyto).eq.1) layp(ixto,iyto)=layp(ix,iy)
            ixto=-ix/2+3*(iy-1)/2+2
            iyto=-ix/2-(iy-1)/2-1
            if(indbetwn(ixto,iyto).eq.1) layp(ixto,iyto)=layp(ix,iy)
          enddo
        enddo
      endif
#ifdef DBG
      write(111,*) 'ipoint'
      do iy=icordys-1,icordye+1,2
        write(111,'(100i3)') (layp(ix,iy),ix=icordxs-2,icordxe+2,2)
      enddo
      write(111,*)
#endif
!
! neighbor node linkage  
      do iy=icordys,icordye
        do ix=icordxs,icordxe
          icel=iaass(ix,iy) 
          if(icel.ne.0) then 
            neignd(1,icel)=layh(ix-4,iy)  
            neignd(2,icel)=layh(ix-2,iy-2)  
            neignd(3,icel)=layh(ix+2,iy-2)  
            neignd(4,icel)=layh(ix+4,iy)  
            neignd(5,icel)=layh(ix+2,iy+2)  
            neignd(6,icel)=layh(ix-2,iy+2)  
            neigpt(1,icel)=layp(ix-2,iy+1) 
            neigpt(2,icel)=layp(ix-2,iy-1) 
            neigpt(3,icel)=layp(ix,iy-1) 
            neigpt(4,icel)=layp(ix+2,iy-1)
            neigpt(5,icel)=layp(ix+2,iy+1) 
            neigpt(6,icel)=layp(ix,iy+1)
            neigsfc(1,icel)=lays(ix-2,iy) 
            neigsfc(2,icel)=lays(ix-1,iy-1) 
            neigsfc(3,icel)=lays(ix+1,iy-1) 
            neigsfc(4,icel)=lays(ix+2,iy) 
            neigsfc(5,icel)=lays(ix+1,iy+1) 
            neigsfc(6,icel)=lays(ix-1,iy+1) 
            neigjin(1,icel)=4
            neigjin(2,icel)=5
            neigjin(3,icel)=6
            neigjin(4,icel)=1
            neigjin(5,icel)=2
            neigjin(6,icel)=3
            wtdhat(1,icel)=-1.
            wtdhat(2,icel)=-1.
            wtdhat(3,icel)=-1.
            wtdhat(4,icel)=1.
            wtdhat(5,icel)=1.
            wtdhat(6,icel)=1.
          endif
        enddo
      enddo
!
      do iy=icordys,icordye,2
        do ix=icordxs-2,icordxe+2,2
          isfc=iasurf(ix,iy) 
          if(isfc.ne.0) then 
            neigsnd(1,isfc)=layh(ix-2,iy)  
            neigsnd(2,isfc)=layh(ix+2,iy) 
            neigsnd(4,isfc)=4  
            neigsnd(5,isfc)=1 
          endif
        enddo
      enddo
      do iy=icordys-1,icordye+1,2
        do ix=icordxs-1,icordxe+1,2
          isfc=iasurf(ix,iy) 
          if(isfc.ne.0) then 
            neigsnd(1,isfc)=layh(ix-1,iy-1)  
            neigsnd(2,isfc)=layh(ix+1,iy+1)
            neigsnd(4,isfc)=5  
            neigsnd(5,isfc)=2 
            if(neigsnd(1,isfc).eq.0) then
              neigsnd(1,isfc)=layh(ix+1,iy-1)
              neigsnd(4,isfc)=6  
            endif
            if(neigsnd(2,isfc).eq.0) then
              neigsnd(2,isfc)=layh(ix-1,iy+1)
              neigsnd(5,isfc)=3 
            endif
          endif
        enddo
      enddo
      do isfc=1,nxsfc
        neigsnd(3,isfc)=12
        do id=1,2
          if(neigsnd(id,isfc).eq.0) neigsnd(3,isfc)=mp(id)
        enddo
      enddo
!
      do iz=1,nz
        neigz(1,iz)=iz-1
        neigz(2,iz)=iz+1
      enddo 
      do iz=1,nz
        neigsfcz(1,iz)=iz
        neigsfcz(2,iz)=iz+1
      enddo 
      neigz(1,1)=0 
      neigz(2,nz)=0 
      do iz=1,nz+1
        neigsndz(1,iz)=iz-1
        neigsndz(2,iz)=iz
        neigsndz(3,iz)=12
        neigsndz(4,iz)=2
        neigsndz(5,iz)=1
      enddo 
      neigsndz(1,1)=0
      neigsndz(3,1)=2
      neigsndz(2,nz+1)=0
      neigsndz(3,nz+1)=1
!
      if(isolang.eq.30) then
        do ix=4,icordxe,4
          icel=layh(ix,0) 
          if(icel.ne.0) then 
            neigjin(2,icel)=3
            neigjin(3,icel)=2
            wtdhat(2,icel)=1.
            wtdhat(3,icel)=1.
          endif
        enddo
        do iy=2,icordye,2
          icel=0
          ix=3*iy
          if((nxfc-ix)*(ix+nxfc).ge.0) icel=layh(ix,iy) 
          if(icel.ne.0) then 
            neigjin(1,icel)=5
            neigjin(5,icel)=1
            neigjin(6,icel)=6
            wtdhat(6,icel)=-1.
          endif
        enddo
        do iy=0,icordye,2
          icel=0
          ix=3*iy+4
          if((nxfc-ix)*(ix+nxfc).ge.0) icel=layh(ix,iy) 
          if(icel.ne.0) then 
            neigjin(6,icel)=6
          endif
        enddo
        do it=1,6
          neigjin(it,layh(0,0))=1
          wtdhat(it,layh(0,0))=1.
        enddo
        icel=layh(4,0) 
        if(icel.ne.0) then 
          neigjin(2,icel)=6
          wtdhat(2,icel)=1.
        endif
        do iy=1,icordye+1,2
          isfc=0
          ix=3*iy
          if((nxfc-ix)*(ix+nxfc).ge.0) isfc=iasurf(3*iy,iy) 
          if(isfc.ne.0) then 
            neigsnd(5,isfc)=6  
          endif
        enddo
      endif
!
      if(isolang.eq.60) then
        if(isymtype.eq.1) then
          do ix=0,icordxe,4
            icel=layh(ix,0) 
            if(icel.ne.0) then 
              neigjin(2,icel)=6
              neigjin(3,icel)=1
              wtdhat(3,icel)=1.
            endif
          enddo
          icel=layh(4,0) 
          if(icel.ne.0) then 
            neigjin(6,icel)=2
          endif
          do iy=2,icordye,2
            icel=layh(iy+4,iy) 
            if(icel.ne.0) then 
              neigjin(1,icel)=3
              neigjin(6,icel)=2
            endif
          enddo
          icel=layh(0,0) 
          do it=1,6
            neigjin(it,icel)=1
            wtdhat(it,icel)=1.
          enddo
          do iy=2,icordye,2
            isfc=iasurf(iy+2,iy) 
            if(isfc.ne.0) then 
              neigsnd(4,isfc)=3  
            endif
          enddo
          do iy=1,icordye+1,2
            isfc=iasurf(iy+2,iy) 
            if(isfc.ne.0) then 
              neigsnd(5,isfc)=2  
            endif
          enddo
        endif
        if(isymtype.eq.2) then
          do ix=0,icordxe,4
            icel=layh(ix,0) 
            if(icel.ne.0) then 
              neigjin(2,icel)=3
              neigjin(3,icel)=2
              wtdhat(2,icel)=1.
              wtdhat(3,icel)=1.
            endif
          enddo
          do iy=0,icordye,2
            icel=layh(iy,iy) 
            if(icel.ne.0) then 
              neigjin(1,icel)=6
              neigjin(6,icel)=1
            endif
          enddo
          do it=1,6,2
            neigjin(it,layh(0,0))=2
            neigjin(it+1,layh(0,0))=1
            wtdhat(it,layh(0,0))=1.
          enddo
        endif
      endif
!
      if(isolang.eq.120) then
        do ix=0,icordxe,4
          icel=layh(ix,0) 
          if(icel.ne.0) then 
            neigjin(2,icel)=1
            neigjin(3,icel)=2
            wtdhat(2,icel)=1.
            wtdhat(3,icel)=1.
          endif
        enddo
        icel=layh(2,2) 
        if(icel.ne.0) then 
          neigjin(1,icel)=2
        endif
        do iy=4,icordye,2
          icel=layh(-iy+4,iy) 
          if(icel.ne.0) then 
            neigjin(1,icel)=2
            neigjin(2,icel)=3
          endif
        enddo
        icel=layh(0,0) 
        do it=1,6,2
          neigjin(it,icel)=2
          neigjin(it+1,icel)=1
          wtdhat(it,icel)=1.
        enddo
        do iy=2,icordye,2
          isfc=iasurf(-iy+2,iy) 
          if(isfc.ne.0) then 
            neigsnd(4,isfc)=2  
          endif
        enddo
        do iy=3,icordye+1,2
          isfc=iasurf(-iy+2,iy) 
          if(isfc.ne.0) then 
            neigsnd(4,isfc)=3  
          endif
        enddo
      endif
!=====================================================================
!            set parameters for point flux solver(pntslv.F)
!=====================================================================
      do ih=1,nassy
        do it=1,6
          iptmy=neigpt(it,ih)
          do i2=1,6
            if(i2.ne.it) then
              iptyou=neigpt(i2,ih)
              do i3=1,neignpt(iptmy)
                if(neigppt(i3,iptmy).eq.iptyou) goto 100
              enddo
              neignpt(iptmy)=neignpt(iptmy)+1
              neigppt(neignpt(iptmy),iptmy)=iptyou
  100         continue
            endif
          enddo
        enddo
      enddo
      do ih=1,nassy
        do it=1,6
          iptmy=neigpt(it,ih)
          do i2=1,6
            if(i2.ne.it) then
              idif=i2-it
              if(idif.lt.0) idif=idif+6
              iptyou=neigpt(i2,ih)
              do i3=1,neignpt(iptmy)
                if(neigppt(i3,iptmy).eq.iptyou) then
                  imatid(idif,it,ih)=i3
                  goto 200
                endif
              enddo
              print *, "Error when searching position of point"
  200         continue
            endif
          enddo
        enddo
      enddo
!
! generation of core edit data
! icoreys,icoreye,icorexs(-mxnyfc:mxnyfc),icorexe(-mxnyfc:mxnyfc),icoref
      icoreys=icordye
      icoreye=icordys
      icoref=icordxe
      do iy=icordys,icordye,2
        icorexs(iy)=icordxe
        icorexe(iy)=icordxs
        do ix=icordxs,icordxe,2
          if(iaass(ix,iy).ne.0) then
            if(iy.lt.icoreys) icoreys=iy
            if(iy.gt.icoreye) icoreye=iy
            if(ix.lt.icorexs(iy)) icorexs(iy)=ix
            if(ix.gt.icorexe(iy)) icorexe(iy)=ix
            if(ix.lt.icoref) icoref=ix
          endif
        enddo
      enddo
!
! generation of core edit data
! icoreyps,icoreype,icorexps(-mxnyfc:mxnyfc),icorexpe(-mxnyfc:mxnyfc),icorepf
      icoreyps=icordye+2
      icoreype=icordys-2
      icorepf=icordxe+2
      do iy=icordys-1,icordye+1,2
        icorexps(iy)=icordxe+2
        icorexpe(iy)=icordxs-2
        do ix=icordxs-2,icordxe+2,2
          if(iapoint(ix,iy).ne.0) then
            if(iy.lt.icoreyps) icoreyps=iy
            if(iy.gt.icoreype) icoreype=iy
            if(ix.lt.icorexps(iy)) icorexps(iy)=ix
            if(ix.gt.icorexpe(iy)) icorexpe(iy)=ix
            if(ix.lt.icorepf) icorepf=ix
          endif
        enddo
      enddo
!
! generation of core edit data
! icoreys,icoreye,icorexs(-mxnyfc:mxnyfc),icorexe(-mxnyfc:mxnyfc),icoref
      icoreyss=icordye+1
      icoreyse=icordys-1
      icoresf=icordxe+1
      do iy=icordys-1,icordye+1
        icorexss(iy)=icordxe+1
        icorexse(iy)=icordxs-1
        do ix=icordxs-2,icordxe+2
          if(iasurf(ix,iy).ne.0) then
            if(iy.lt.icoreyss) icoreyss=iy
            if(iy.gt.icoreyse) icoreyse=iy
            if(ix.lt.icorexss(iy)) icorexss(iy)=ix
            if(ix.gt.icorexse(iy)) icorexse(iy)=ix
            if(ix.lt.icoresf) icoresf=ix
          endif
        enddo
      enddo
!
! generation of core edit data for fuel assembly
! icoreyfs,icoreye,icorexs(-mxnyfc:mxnyfc),icorexe(-mxnyfc:mxnyfc),icoref
!
      icoreyfs=icordye
      icoreyfe=icordys
      icoreff=icordxe

      do iy=icordys,icordye,2
        icorexfs(iy)=icordxe
        icorexfe(iy)=icordxs
        do ix=icordxs,icordxe,2
          if(iaass(ix,iy).ne.0) then
            if(iffuela(iassytyp(iaass(ix,iy)))) then
              if(iy.lt.icoreyfs) icoreyfs=iy
              if(iy.gt.icoreyfe) icoreyfe=iy
              if(ix.lt.icorexfs(iy)) icorexfs(iy)=ix
              if(ix.gt.icorexfe(iy)) icorexfe(iy)=ix
              if(ix.lt.icoreff) icoreff=ix
            endif
          endif
        enddo
      enddo
!
! Generation of information for Condensation of Matrix cmat
!
      call premat
!
      return
    end subroutine
!=======================================================================================!
    function indbetwn(ix,iy)

      use param

      include 'global.h'
      include 'geomh.h'
      include 'geomhfc.h'
!
      indbetwn=0
      indx=(nxfc-ix)*(ix+nxfc)
      indy=(nyfc-iy)*(iy+nyfc)
      if(indx.ge.0.and.indy.ge.0) indbetwn=1
!
      return
    end function
!=======================================================================================!
    subroutine premat
!
! Generation of information for Condensation of Matrix cmat
!   ipntr(0,ia) : The number of elements in cmat(ig,nn,ia,iz)
!
      use param

      include 'global.h'
      include 'geom.h'
      include 'geomh.h'
      include 'xsec.h'
      include 'defhex.h'
      include 'defpnt.h'
      include 'lscoefh.h'
      include 'deffg.h'

      dimension irtmp(6)
! radial direction condensation
      do ih=1,nassy
        ipntr(0,ih)=0
        do ir=0,6
          ineigcond(0,ir,ih)=0
        enddo
        do irfr=1,6
          inum=neignd(irfr,ih)
          if(inum.ne.0) then
            if(inum.eq.ih) then
              inext=ineigcond(0,0,ih)+1
              ineigcond(inext,0,ih)=irfr
              ineigcond(0,0,ih)=inext
            else
              do irto=1,ipntr(0,ih)
                if(inum.eq.ipntr(irto,ih)) then
                  inext=ineigcond(0,irto,ih)+1
                  ineigcond(inext,irto,ih)=irfr
                  ineigcond(0,irto,ih)=inext
                  goto 100
                endif
              enddo
              inext=ipntr(0,ih)+1
              ipntr(inext,ih)=inum
              ineigcond(0,inext,ih)=1
              ineigcond(1,inext,ih)=irfr
              ipntr(0,ih)=inext
            endif
  100       continue
          endif
        enddo
        do ir=1,ipntr(0,ih)
          imin=ir
          do ir2=ir+1,ipntr(0,ih)
            if(ipntr(ir2,ih).lt.ipntr(imin,ih)) imin=ir2
          enddo
          if(imin.ne.ir) then
            do ir2=1,ineigcond(0,ir,ih)
              irtmp(ir2)=ineigcond(ir2,ir,ih)
            enddo
            do ir2=1,ineigcond(0,imin,ih)
              ineigcond(ir2,ir,ih)=ineigcond(ir2,imin,ih)
            enddo
            do ir2=1,ineigcond(0,ir,ih)
              ineigcond(ir2,imin,ih)=irtmp(ir2)
            enddo
            ii=ineigcond(0,ir,ih)
            ineigcond(0,ir,ih)=ineigcond(0,imin,ih)
            ineigcond(0,imin,ih)=ii
            ii=ipntr(ir,ih)
            ipntr(ir,ih)=ipntr(imin,ih)
            ipntr(imin,ih)=ii
          endif
        enddo
        ilubnd(ih)=0
        do ir=1,ipntr(0,ih)
          if(ipntr(ir,ih).lt.ih) ilubnd(ih)=ir
          iastopnt(ipntr(ir,ih),ih)=ir
        enddo
        iastopnt(ih,ih)=7
      enddo
!
      return
    end subroutine