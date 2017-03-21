! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
    subroutine scanlay(ifile,isyminp)
!
      use param
      USE MASTERXSL  ! 2014.07.28_PKB

!
! scan lay(loading pattern) card and generate maximum value for radial grid
!
      include 'global.h'     
      include 'files.h'
      include 'cards.h'      
      include 'xsec.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'itrcntl.h'  
      include 'ff.h'    
      include 'thcntl.inc'
      include 'perturb.inc'
      include 'thlink.inc' ! FREK/MARS
      logical ifnumeric

      character*65 amesg
      character*4  astr
      logical,save :: first=TRUE
      dimension ionerow(512)

      if(first) then
        allocate(layh(-nxfc:nxfc,-nyfc:nyfc))
        allocate(layp(-nxfc:nxfc,-nyfc:nyfc))
        allocate(lays(-nxfc:nxfc,-nyfc:nyfc))
        layh=0
        layp=0
        lays=0
        ndivhs=1
        first=.FALSE.
      endif
      if(isyminp.eq.30) then
         nxrow=nring/2+1
      else if(isyminp.le.120) then
         nxrow=nring+1
      else
         nxrow=2*nring+1
      endif
! PKB ADDITION
      ALLOCATE(APNUM(NXROW))
! ADDED END

      nat=0 ! initialize number of assembly types to zero
!
! read radial core configuation
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
      nassyinp=0
      do while(ir.lt.nxrow)
        read(ifile,'(a512)') oneline
        nxcol=nfields(oneline)
        if(sc(1).eq.'!' .or. nxcol.eq.0) go to 200 !skip a line
!
! check if next card is read
        do i=1,mxncol
          if(sc(i).ne.' ') then
             iascii=ichar(sc(i))
             if((iascii-48)*(iascii-57) .le. 0) then
               go to 100 !numeric
             else
               backspace(ifile)
               go to 300 !string
             endif
          endif
        enddo
  100   continue
        ir=ir+1
        read(oneline,*) (ionerow(ic),ic=1,nxcol) 
        APNUM(IR) = NXCOL
        nassyinp=nassyinp+nxcol
        do ic=1,nxcol
           nat=max(nat,ionerow(ic))
        enddo
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
          layh(ix1,iy1)=ionerow(ic)
          ix1=ix1+4 
        enddo
        iy1=iy1+2
  200   continue        
      enddo
  300 continue
      NASSEM = NASSYINP
      nxrow=ir        
      icordxs=-icordxe
      icordys=-icordye
!
! generate load patterb for 360 degree
!
      if(isyminp.eq.30) then
        do iy=0,icordye,2
          do ix=3*iy,icordxe,4
            layh(ix/2+3*iy/2,ix/2-iy/2)=layh(ix,iy)
          enddo
        enddo
      endif
!
      if(isyminp.le.60) then
        do iy=0,icordye,2
          do ix=iy,icordxe,4
            if(isymtype.eq.1) then
              layh(ix/2-3*iy/2,ix/2+iy/2)=layh(ix,iy)
            else
              layh(-ix/2+3*iy/2,ix/2+iy/2)=layh(ix,iy)
            endif
          enddo
        enddo
      endif
!
      if(isyminp.le.120) then
        do iy=0,icordye,2
          do ix=-iy,icordxe,4
            layh(-ix/2-3*iy/2,ix/2-iy/2)=layh(ix,iy)
            layh(-ix/2+3*iy/2,-ix/2-iy/2)=layh(ix,iy)
          enddo
        enddo
      endif
!
! surface and point generation for 360 degree
      do iy=icordys,icordye
        do ix=icordxs,icordxe
          if(layh(ix,iy).ne.0) then
            layp(ix-2,iy-1)=1 
            layp(ix,iy-1)=1 
            layp(ix+2,iy-1)=1 
            layp(ix-2,iy+1)=1 
            layp(ix,iy+1)=1 
            layp(ix+2,iy+1)=1 
            lays(ix-1,iy-1)=1 
            lays(ix+1,iy-1)=1 
            lays(ix-2,iy)=1 
            lays(ix+2,iy)=1 
            lays(ix-1,iy+1)=1 
            lays(ix+1,iy+1)=1 
          endif
        enddo
      enddo
!
! hex node generation and numbering 
!
      if(isolang.lt.360) then
        do iy=icordys-1,-1
          do ix=icordxs-2,icordxe+2
            layh(ix,iy)=0
            layp(ix,iy)=0
            lays(ix,iy)=0
          enddo
        enddo
        do ix=icordxs-2,-1
          layh(ix,0)=0
        enddo
      endif
      if(isolang.eq.120) then
        do iy=1,icordye
          do ix=icordxs,-iy
            layh(ix,iy)=0
          enddo
        enddo
        do iy=0,icordye+1
          do ix=icordxs-2,-iy-1
            layp(ix,iy)=0
          enddo
          do ix=icordxs-2,-iy
            lays(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.60) then
        do iy=1,icordye
          ix2=iy-1
          if(isymtype.eq.1) ix2=iy
          do ix=icordxs,ix2
            layh(ix,iy)=0
          enddo
        enddo
        do iy=0,icordye+1
          do ix=icordxs-2,iy-1
            layp(ix,iy)=0
          enddo
          ix2=iy-1
          if(isymtype.eq.1) ix2=iy
          do ix=icordxs-2,ix2
            lays(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.30) then
        do iy=1,icordye
          do ix=icordxs,iy*3-1
            if(ix.le.icordxe) layh(ix,iy)=0
          enddo
        enddo
        do iy=0,icordye+1
          do ix=icordxs-2,iy*3-2
            if(ix.le.icordxe) then
              layp(ix,iy)=0
              lays(ix,iy)=0
            endif
          enddo
        enddo
      endif
!
      inumh=0
      inums=0
      inump=0
!
      nx=2*nring+1
      ny=nring+1
      if(isolang.eq.360) then
         nx=4*nring+1
         ny=nx
      elseif(isolang.eq.120) then
         nx=3*nring+1
      elseif(isolang.eq.30) then
         ny=nring/2+1
      endif
      ibegh=0
      iendh=0
      jbegh=0
      jendh=0
!
      do iy=icordys-1,icordye+1
        do ix=icordxs-2,icordxe+2
          if(layh(ix,iy).ne.0) then
            inumh=inumh+1
            ibegh=min(ibegh,ix)
            iendh=max(iendh,ix)
            jbegh=min(jbegh,iy)
            jendh=max(jendh,iy)
          endif
          if(layp(ix,iy).ne.0) then
            inump=inump+1 
          endif
          if(lays(ix,iy).ne.0) then
            inums=inums+1 
          endif
        enddo 
      enddo 
      nassy=inumh
      ncorn=inump
      nsurf=inums
!
      nasyx=(iendh-ibegh)/2+1
      nasyy=(jendh-jbegh)/2+1
      nassy=inumh
      nx=nasyx
      ny=nasyy
      nxy=nassy
      nxya=nassy
      nchan=nassy
!      
      return
    end subroutine




