    subroutine validateinput
      
      use param
      
      include 'global.h'
      include 'geom.h'
      include 'thgeom.inc'      
      include 'thcntl.inc'
      
      integer :: nchanx,nchany
      real :: fnchan1a
      
! t/h mesh validation
      if(fdbk) then
        nchanx=0;nchany=0
        do ia=1,nxa
          nrest=mod(nmeshx(ia),nthx(ia))
          if(nrest .ne. 0) then
            call terminate("T/H meshes doesn't corresponde to neutronic meshs")
          endif
        enddo
        do ja=1,nya
          nrest=mod(nmeshy(ja),nthy(ja))
          if(nrest .ne. 0) then
            call terminate("T/H meshes doesn't corresponde to neutronic meshs")
          endif
        enddo
      endif

    end subroutine