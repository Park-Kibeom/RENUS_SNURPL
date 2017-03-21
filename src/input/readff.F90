    subroutine readff

      use param
      use allocs
      
      include 'global.h'
      include 'geom.h'
      include 'xsec.h'
      include 'nodal.h'
      include 'pinpr.h'
      include 'ff.h'
      include 'first.h'
      include 'cards.h'
      include 'files.h'
!
      real,allocatable,save :: dumb(:,:) !dmm
      character*4 ntype
      character*4 astr
      character*1024 tmpstr, dum
      save ipff,nsize
      logical flagnpin,flagcomp
      logical ifnumeric
!
! set defaults
      if(first) then
!
        ntype='NODE' ! pff for a quarter assy
        do i=1,ncomp !dmm
          iffset(i)=1
        enddo
!
        do i=1,mxnfcard
          iffcard(i)=FALSE
        enddo
! to check the existence of cardname NPIN_SIDE, PFF_COMP
        flagnpin=FALSE
        flagcomp=FALSE
!
! allocate temporary memory
        allocate(dumb(npin,npin))
!
        first=FALSE
      endif
!
      iffile=FALSE
      indev=io5    
!
  100 continue
      do while(probe.ne.DOT)
        read(indev,'(a512)',end=1000) oneline
        write(io8,'(a)') trim(oneline)  
        if(probe.eq.BANG .or. oneline.eq.BLANK .or.ifnumeric(oneline)) cycle
        if(probe.eq.DOT .or. probe.eq.SLASH) exit
        if(probe.ne.BLANK) then
          backspace(indev)
          backspace(io8)
          return
        endif
!#ifdef MOXBENCH
!        npin=17
!        fntpin=npin*npin
!        nffset=ncomp
!        do icomp=1,ncomp
!            iffset(icomp)=icomp
!        enddo
!        call dmalloc(ff,npin,npin,ng,nffset)  ! ff for unrodded set
!        backspace(indev)
!        imox =1000
!! FIXME 6 is hard coded.
!        do i=1,6
!          read(indev,'(a512)',end=1000) oneline      
!          read(oneline,*) cardname, icbase
!          call getfn(oneline,3,filename)
!          call openfile(imox,TRUE,FALSE,filename)
!          call readbenchff(imox, icbase)
!          close(imox)
!        enddo
!        go to 1000
!#else
!        read(oneline,*) cardname
!        call toupper(cardname)
!        if(cardname.eq.'FILE') then
!          indev=io5+100
!          call openlf(indev,oneline)
!          iffile=TRUE
!          go to 100
!        endif
!#endif        
        read(oneline,*) cardname
        call toupper(cardname)
        if(cardname.eq.'FILE') then
          indev=io5+100
          call openlf(indev,oneline)
          iffile=TRUE
          go to 100
        endif
        ndataf=nfields(oneline)-1
        select case(cardname) 
          case('NPIN_SIDE')
              flagnpin=TRUE
              read(oneline,*) cardname,npin
              fntpin=npin*npin
              nsize=(npin+1)/2
! 2013_12_10 . commented              
          !case('PFF_GEOM')
          !    if(ndataf.gt.1) then
          !      read(oneline,*) cardname,ipff,(idum(i),i=1,ndataf-1)
          !    else
          !      mesg='Composition Numbers Not Specified for PFF'
          !      call terminate(mesg)
          !    endif
          !    if(ipff.lt.0) then
          !      mesg='Negative PFF Set Number'
          !      call terminate(mesg)
          !    elseif(ipff.gt.nffset) then
          !      write(astr,'(i4)') nffset
          !      mesg='Composition Number > '//astr
          !      call terminate(mesg)
          !    endif
          !    call intexpand(ndataf-1,ipff,idum,iffset,ncomp)
          !    nffset=max(nffset,ipff)
! end              
          case('PFF_COMP')
              flagcomp=TRUE
              if(ndataf.gt.1) then
                read(oneline,*) cardname,ipff,(idum(i),i=1,ndataf-1)
              else
                mesg='Composition Numbers Not Specified for PFF'
                call terminate(mesg)
              endif
              if(ipff.lt.0) then
                mesg='Negative PFF Set Number'
                call terminate(mesg)
              elseif(ipff.gt.nffset) then
                write(astr,'(i4)') nffset
                mesg='Composition Number > '//astr
                call terminate(mesg)
              endif
              call intexpand(ndataf-1,ipff,idum,iffset,ncomp)
              nffset=max(nffset,ipff)
          case('PFF_UNRODD','PFF_RODDED')
              ntype='NODE' ! pff for a quarter assy
              if(flagnpin.and.flagcomp) then
                if(ndataf .gt. 1) then
                  read(oneline,*) cardname,m,ntype
                else
                  read(oneline,*) cardname,m
                endif
                call toupper(cardname)
                call toupper(ntype)
              else
                mesg='HFF Set Number or Npin Not Specified'
                call terminate(mesg)
              endif
              if(ntype.eq.'FULL') then
                if(cardname.eq.'PFF_UNRODD') then ! unrodded
                  do j=1,npin
                    call read1more(indev,oneline,ndatafnew)
                    read(oneline,*) (ff(i,j,m,ipff),i=1,npin)
                  enddo
                else                ! rodded
                  do j=1,npin
                    call read1more(indev,oneline,ndatafnew)
                    read(oneline,*) (ffr(i,j,m,ipff),i=1,npin)
                  enddo
                endif
              else
                do j=1,nsize
                  call read1more(indev,oneline,ndatafnew)
                  read(oneline,*) (dumb(i,j),i=1,nsize)
                enddo
                if(cardname.eq.'PFF_UNRODD') then ! unrodded
!!define DBG
!#ifdef DBG
!                  write(111,*)
!                  do j=1,nsize  
!                    write(111,'(100(f6.4,1x))') (dumb(i,j),i=1,nsize)
!                  enddo
!#endif
                  call ffexpand(m,nsize,dumb,ipff)
!#ifdef DBG      
!                  write(111,*) "cardname=", cardname
!                  write(111,*) "ipff=",ipff,"group=",m
!                  do ipy=1,npin
!                    write(111,'(100(f6.4,1x))') (ff(ipx,ipy,m,ipff),ipx=1,npin)
!                  enddo
!                  continue
!#endif
                else                               !rodded
                  call ffrexpand(m,nsize,dumb,ipff)
!#ifdef DBG      
!                  write(111,*) "cardname=", cardname
!                  write(111,*) "ipff=",ipff,"group=",m
!                  do ipy=1,npin
!                    write(111,'(100(f6.4,1x))') (ff(ipx,ipy,m,ipff),ipx=1,npin)
!                  enddo
!                  continue
!#endif
                endif
                continue
              endif
          case default
              call terminate(trim(cardname)//' Card Not Allowed')
        end select
      enddo
 1000 continue
!
! reset the input device after done with local file
      if(iffile) then
         close(indev)
         indev=io5
         iffile=FALSE
!
! return to the next card in the input file
         goto 100
      endif
!
      return
      
    end subroutine

!----------------------------------------------------------
    subroutine ffexpand(m,nsize,dumb,ipff)

      use param
      
      include 'global.h'
      include 'geom.h'
      include 'xsec.h'
      include 'nodal.h'
      include 'pinpr.h'
      include 'ff.h'
      
      integer,intent(in) :: m,nsize,ipff
      real,intent(in) :: dumb(npin,npin)
! if form function is given in upper left 
      if((nsize*2-1).eq.npin) then
        do i=1,nsize-1
            do j=1,nsize-1
               ff(i,j,m,ipff)=dumb(i,j)
               ff(2*nsize-i,j,m,ipff)=dumb(i,j)
               ff(i,2*nsize-j,m,ipff)=dumb(i,j)
               ff(2*nsize-i,2*nsize-j,m,ipff)=dumb(i,j)
            enddo
         enddo
         do i=1,nsize-1
            ff(nsize,2*nsize-i,m,ipff)=dumb(nsize,i)
            ff(2*nsize-i,nsize,m,ipff)=dumb(i,nsize)
            ff(nsize,i,m,ipff)=dumb(nsize,i)
            ff(i,nsize,m,ipff)=dumb(i,nsize)
         enddo
         ff(nsize,nsize,m,ipff)=dumb(nsize,nsize)
      elseif((nsize*2).eq.npin) then
         do i=1,nsize
            do j=1,nsize
               ff(i,j,m,ipff)=dumb(i,j)
               ff(2*nsize+1-i,j,m,ipff)=dumb(i,j)
               ff(i,2*nsize+1-j,m,ipff)=dumb(i,j)
               ff(2*nsize+1-i,2*nsize+1-j,m,ipff)=dumb(i,j)
            enddo
         enddo
      endif

!!define lowerright
!#ifdef lowerright
!      if((nsize*2-1).eq.npin) then
!        do j=0,nsize-1
!          do i=0,nsize-1
!            ff(i+nsize,j+nsize,m,ipff)=dumb(i+1,j+1)
!            ff(i+nsize-2*i,j+nsize,m,ipff)=dumb(i+1,j+1)
!            ff(i+nsize,j+nsize-2*j,m,ipff)=dumb(i+1,j+1)
!            ff(i+nsize-2*i,j+nsize-2*j,m,ipff)=dumb(i+1,j+1)
!          enddo
!        enddo
!      elseif((nsize*2).eq.npin) then
!        do j=1,nsize
!          do i=1,nsize
!            ff(i+nsize,j+nsize,m,ipff)=dumb(i,j)
!            ff(nsize-(i-1),j+nsize,m,ipff)=dumb(i,j)
!            ff(i+nsize,nsize-(j-1),m,ipff)=dumb(i,j)
!            ff(nsize-(i-1),nsize-(j-1),m,ipff)=dumb(i,j)
!          enddo
!        enddo        
!      endif
!#endif
 !
      return
      
    end subroutine

!----------------------------------------------------------
    subroutine ffrexpand(m,nsize,dumb,ipff)

      use param
      
      include 'global.h'
      include 'geom.h'
      include 'xsec.h'
      include 'nodal.h'
      include 'pinpr.h'
      include 'ff.h'
      
      integer,intent(in) :: m,nsize,ipff
      real,intent(in) :: dumb(npin,npin)


      if((nsize*2-1).eq.npin) then
        do i=1,nsize-1
            do j=1,nsize-1
               ffr(i,j,m,ipff)=dumb(i,j)
               ffr(2*nsize-i,j,m,ipff)=dumb(i,j)
               ffr(i,2*nsize-j,m,ipff)=dumb(i,j)
               ffr(2*nsize-i,2*nsize-j,m,ipff)=dumb(i,j)
            enddo
         enddo
         do i=1,nsize-1
            ffr(nsize,2*nsize-i,m,ipff)=dumb(nsize,i)
            ffr(2*nsize-i,nsize,m,ipff)=dumb(i,nsize)
            ffr(nsize,i,m,ipff)=dumb(nsize,i)
            ffr(i,nsize,m,ipff)=dumb(i,nsize)
         enddo
         ffr(nsize,nsize,m,ipff)=dumb(nsize,nsize)
      elseif((nsize*2).eq.npin) then
         do i=1,nsize
            do j=1,nsize
               ffr(i,j,m,ipff)=dumb(i,j)
               ffr(2*nsize+1-i,j,m,ipff)=dumb(i,j)
               ffr(i,2*nsize+1-j,m,ipff)=dumb(i,j)
               ffr(2*nsize+1-i,2*nsize+1-j,m,ipff)=dumb(i,j)
            enddo
         enddo
      endif

!!define lowerright2
!#ifdef lowerright2
!      if((nsize*2-1).eq.npin) then
!        do j=0,nsize-1
!          do i=0,nsize-1
!            ffr(i+nsize,j+nsize,m,ipff)=dumb(i+1,j+1)
!            ffr(i+nsize-2*i,j+nsize,m,ipff)=dumb(i+1,j+1)
!            ffr(i+nsize,j+nsize-2*j,m,ipff)=dumb(i+1,j+1)
!            ffr(i+nsize-2*i,j+nsize-2*j,m,ipff)=dumb(i+1,j+1)
!          enddo
!        enddo
!      elseif((nsize*2).eq.npin) then
!        do j=1,nsize
!          do i=1,nsize
!            ffr(i+nsize,j+nsize,m,ipff)=dumb(i,j)
!            ffr(nsize-(i-1),j+nsize,m,ipff)=dumb(i,j)
!            ffr(i+nsize,nsize-(j-1),m,ipff)=dumb(i,j)
!            ffr(nsize-(i-1),nsize-(j-1),m,ipff)=dumb(i,j)
!          enddo
!        enddo  
!      endif
!#endif
!
      return
      
    end subroutine
