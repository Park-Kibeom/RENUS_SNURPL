! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
!23456789112345678921234567893123456789412345678951234567896123456789712
    subroutine readbankh(iroute,ifile,isyminp)
!
! Read Lay(loading pattern) Card and Expand to 360 degree loading Pattern
!
      use param

      include 'global.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'cntl.h'
      include 'itrcntl.h'
      include 'files.h'
      include 'cards.h'
      include 'defhex.h'
      include 'trancntl.inc'

      integer,allocatable,save :: laybank(:,:)
      logical,save :: first=TRUE
!
      if(first) then
         allocate(laybank(-nxfc:nxfc,-nyfc:nyfc))
         laybank=0
         first=.FALSE.
      endif
 
      if(iroute.eq.2) go to 1000
      
      ival=0
      layh=ival
      iy1=0
      if(isyminp.eq.360) then
        nhalfs=nxrow/2
        iy1=-2*nhalfs
      endif
      
      ir=0
      do while(ir.lt.nxrow)
        read(ifile,'(a512)') oneline
        write(io8,form7(oneline,ncol)) oneline
        nxcol=nfields(oneline)
        if(sc(1).eq.'!' .or. nxcol.eq.0) go to 200 !skip a line
        ir=ir+1
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
          iy2=(ix1+iy1)/2
          laybank(ix1,iy1)=idat(ic)
          ix1=ix1+4 
        enddo
        iy1=iy1+2
  200   continue        
      enddo
!
! generate loading pattern for 360 degree
!
      if(isyminp.eq.30) then
        do iy=0,icordye,2
          do ix=3*iy,icordxe,4
            ixt=ix/2+3*iy/2
            iyt=ix/2-iy/2
            if(indbetwn(ixt,iyt).eq.1) laybank(ixt,iyt)=laybank(ix,iy)
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
            if(indbetwn(ixt,iyt).eq.1) laybank(ixt,iyt)=laybank(ix,iy)
          enddo
        enddo
      endif
!
      if(isyminp.le.120) then
        do iy=0,icordye,2
          do ix=-iy,icordxe,4
            ixt=-ix/2-3*iy/2
            iyt=ix/2-iy/2
            if(indbetwn(ixt,iyt).eq.1) laybank(ixt,iyt)=laybank(ix,iy)
            ixt=-ix/2+3*iy/2
            iyt=-ix/2-iy/2
            if(indbetwn(ixt,iyt).eq.1) laybank(ixt,iyt)=laybank(ix,iy)
          enddo
        enddo
      endif
!
      if(isolang.lt.360) then
        do iy=icordys,-1
          do ix=icordxs,icordxe
            laybank(ix,iy)=0
          enddo
        enddo
        do ix=icordxs,-1
          laybank(ix,0)=0
        enddo
      endif
      if(isolang.eq.120) then
        do iy=1,icordye
          do ix=icordxs,-iy
            laybank(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.60) then
        do iy=1,icordye
          ix2=iy-1
          if(isymtype.eq.1) ix2=iy
          do ix=icordxs,ix2
            laybank(ix,iy)=0
          enddo
        enddo
      endif
      if(isolang.eq.30) then
        do iy=1,icordye
          ix2=iy*3-1
          if(ix2.gt.nxfc) ix2=nxfc
          do ix=icordxs,ix2
            laybank(ix,iy)=0
          enddo
        enddo
      endif
!
      return
!
 1000 continue
      ipend=0
      lcrbptr(0)=0
      do icr=1,nrodtyp
        ncrbasy=0
        do iy=icordys,icordye,2
          do ix=icordxs,icordxe,2
            icrb=laybank(ix,iy)
            if(icr.eq.icrb) then
              ncrbasy=ncrbasy+1
	        iasy=iaass(ix,iy)
              irod=ipend+ncrbasy
!              lcrbtola(ipend+ncrbasy)=iaass(ix,iy)
			rodtola(irod)=iasy
	        irodtyp(iasy)=icrb
            endif
          enddo
        enddo
        ipend=ipend+ncrbasy
        lcrbptr(icr)=ipend
      enddo
!
      return
    end subroutine
