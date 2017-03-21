! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
    subroutine readlay(ifile,isyminp)
!
! Read Lay(loading pattern) Card and Expand to 360 degree loading Pattern
!
      use param
      USE MASTERXSL    ! 2014_10_30 . pkb

      include 'global.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'cntl.h'
      include 'itrcntl.h'
      include 'files.h'
      include 'cards.h'
      include 'defhex.h'
!
      iy1=0
      icordys=0
      icordye=2*nxrow-2
      if(isyminp.eq.360) then
        nhalfs=nxrow/2
        iy1=-2*nhalfs
        icordys=iy1
        icordye=-iy1
      endif
      icordxs=0
      icordxe=0
      ir=0
      !idbg=1
      do while(ir.lt.nxrow)
        read(ifile,'(a512)') oneline
        write(io8,form7(oneline,ncol)) oneline
        nxcol=nfields(oneline)
        !print *, idbg,nxcol
        !idbg=idbg+1
        if(sc(1).eq.'!' .or. nxcol.eq.0) go to 200 !skip a line
        ir=ir+1
        radnum = radnum + nxcol   ! 2014_10_17 . pkb        
        read(oneline,*) (idat(ic),ic=1,nxcol) 
        if(isyminp.eq.360) then
          nhalfcol=nxcol/2
          ix1=-4*nhalfcol
          imodv=mod(nxcol,2)
          if(imodv.eq.0) ix1=ix1+2
        endif
        if(isyminp.eq.120) then
          ix1=-(ir-1)*2
          ind=-(2*ix1)/4+1
          if(nxcol.lt.ind) then
            nhalfcol=nxcol/2
            ix1=-4*nhalfcol
            imodv=mod(nxcol,2)
            if(imodv.eq.0) ix1=ix1+2
          endif
        endif
        if(isyminp.eq.60) ix1=(ir-1)*2
        if(isyminp.eq.30) ix1=(ir-1)*6
        do ic=1,nxcol
          if(ix1.lt.icordxs) icordxs=ix1
          if(ix1.gt.icordxe) icordxe=ix1
          iy2=(ix1+iy1)/2
          if(iy2.gt.icordye) icordye=iy2
          iaass(ix1,iy1)=idat(ic)
          ix1=ix1+4 
        enddo
        iy1=iy1+2
  200   continue        
      enddo
      icordxs=-icordxe
      icordys=-icordye
      !pause
!
! generate loading pattern for 360 degree
!
      if(isyminp.eq.30) then
        do iy=0,icordye,2
          do ix=3*iy,icordxe,4
            ixt=ix/2+3*iy/2
            iyt=ix/2-iy/2
            if(indbetwn(ixt,iyt).eq.1) iaass(ixt,iyt)=iaass(ix,iy)
          enddo
        enddo
      endif
!
      if(isyminp.le.60) then
        do iy=0,icordye,2
          do ix=iy,icordxe,4
            iyt=ix/2+iy/2
            if(isymtype.eq.1) then
              ixt=ix/2-3*iy/2
            else
              ixt=-ix/2+3*iy/2
            endif
            if(indbetwn(ixt,iyt).eq.1) iaass(ixt,iyt)=iaass(ix,iy)
          enddo
        enddo
      endif
!
      if(isyminp.le.120) then
        do iy=0,icordye,2
          do ix=-iy,icordxe,4
            ixt=-ix/2-3*iy/2
            iyt=ix/2-iy/2
            if(indbetwn(ixt,iyt).eq.1) iaass(ixt,iyt)=iaass(ix,iy)
            ixt=-ix/2+3*iy/2
            iyt=-ix/2-iy/2
            if(indbetwn(ixt,iyt).eq.1) iaass(ixt,iyt)=iaass(ix,iy)
          enddo
        enddo
      endif
!
! surface and point generation for 360 degree
      do iy=icordys,icordye
        do ix=icordxs,icordxe
          if(iaass(ix,iy).ne.0) then
            iapoint(ix-2,iy-1)=iapoint(ix-2,iy-1)+1 
            iapoint(ix,iy-1)=iapoint(ix,iy-1)+1 
            iapoint(ix+2,iy-1)=iapoint(ix+2,iy-1)+1 
            iapoint(ix-2,iy+1)=iapoint(ix-2,iy+1)+1 
            iapoint(ix,iy+1)=iapoint(ix,iy+1)+1 
            iapoint(ix+2,iy+1)=iapoint(ix+2,iy+1)+1 
            iasurf(ix-1,iy-1)=1 
            iasurf(ix+1,iy-1)=1 
            iasurf(ix-2,iy)=1 
            iasurf(ix+2,iy)=1 
            iasurf(ix-1,iy+1)=1 
            iasurf(ix+1,iy+1)=1 
          endif
        enddo
      enddo
!
      return
    end subroutine
