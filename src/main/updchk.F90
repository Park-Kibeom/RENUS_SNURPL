! added in ARTOS ver. 0.2 ( for acceleration ). 2012_07_04 by SCB
    subroutine updchk
! check the conditional nodal update and multi-group CMFD acceleration

      use param
      use timer
      use sfam_cntl, only : cmfdupd, nodalupd

      include 'global.h'
      include 'times.h'
      include 'xsec.h'
      include 'geom.h'
      include 'files.h'
      include 'trancntl.inc'

      character*3 :: onoff
      logical,save :: first=TRUE
      real,save :: dxst=0.

      dxsn=0.
      dxsc=0.
      do k=1,nz
        do l=1,nxy
          if(k.ge.kfbeg .and. k.le.kfend) then
            do m=1,ng 									
              dxs1=abs(xstf(m,l,k)/xstfn(m,l,k)-1)
              dxsn=max(dxsn,dxs1)

              dxs2=abs(xstf(m,l,k)/xstfc(m,l,k)-1)
              dxsc=max(dxsc,dxs2)
            enddo
          endif
        enddo
      enddo

      if(dxsn.gt.epsxsec) then
        nodalupd=TRUE
        onoff='ON'
      else
        nodalupd=FALSE
        onoff='OFF'		
      endif

      if(dxsc.gt.0.5*epsxsec) then
        cmfdupd=TRUE
      else
        cmfdupd=FALSE	
      endif

      write(mesg, '(a,2f10.6,5x,a)'), ' Delta XSEC & Nodal Update', dxsn, dxsc, onoff
      call message(true,true,mesg)

      return
    end subroutine
