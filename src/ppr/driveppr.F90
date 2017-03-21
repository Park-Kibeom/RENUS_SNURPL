    subroutine driveppr(iftran)
    
      use param
      use allocs
      use ppr,  only : mallocppr, initppr, calhomo

      include 'global.h'
      include 'xsec.h'
      include 'geom.h'
      include 'ffdm.h'
      include 'nodal.h'
      include 'itrcntl.h'
      include 'files.h'
      include 'pinpr.h'
      include 'ff.h'

      logical   :: iftran
      logical,save :: first=TRUE

      call message(TRUE,TRUE,'Calculating Pin Power...')

      if(first) then
        first = FALSE
        call mallocppr(npin)
        call dmalloc(phih,npin,npin,nxya,nz,ng)
        call dmalloc(powvalr,npin,npin,nxya)
        call dmalloc(powval1,npin,npin)
        call dmalloc(powval,npin,npin,ng)
        call dmalloc0(pppeak,1,nxya,0,nz)
        call dmalloc(powtot,ng)        
      endif
      call initppr(iftran)
      
      nmx=nx/nxa
      nmy=ny/nya
      if(nmx .ne. npin .or. nmy.ne.npin) then
        call calhomo
        call homopin
      else
        do k=1,nz
          do ja=1,nya
            do iy=1,npin
              j=(ja-1)*npin+iy
              do ia=nxsa(ja),nxea(ja)
                la=nodela(ia,ja)
                do ix=1,npin
                  i=(ia-1)*npin+ix
                  l=nodel(i,j)
                  do m=1,ng
                    phih(ix,iy,la,k,m)=phif(m,l,k)
                  enddo !m
                enddo !ix
              enddo !ia
            enddo !iy
          enddo !ja
        enddo !k
      endif
!
!FIXME undefine MATLAB_DEBUG       
!#define MATLAB_DEBUG     
!#ifdef MATLAB_DEBUG
!        iodbg=7001
!        write(filename(1), '("debug/",a,"_FLUX.out")') trim(caseid)
!        open(iodbg, file=trim(filename(1)),status='unknown')
!        do m=1,ng
!          write(iodbg, *) 'GROUP : ', m
!          do ja=nya/2+1,nya
!            do iy=1,npin
!              do ia=(nxea(ja)+nxsa(ja))/2+1,nxea(ja)
!                la=nodela(ia,ja)
!                do ix=1,npin
!                  write(iodbg,'(f7.4,$)') phih(ix,iy,la,1,m)
!                enddo
!              enddo !ia
!              write(iodbg,*)
!            enddo !iy
!          enddo !ja
!        enddo
!        close(iodbg)
!#endif
!
      call multiplyff
!          
      return
!
    end subroutine
