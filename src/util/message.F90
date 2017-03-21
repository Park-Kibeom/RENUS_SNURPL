    subroutine message(ifcpu,ifdisp,amesg)
!
      use param
      use timer
!
      include 'global.h'
      include 'cntl.h'
      include 'files.h'
      include 'times.h'
!
      character  amesg*(*)
      character*1 sc(128)
      character*10 form
      character*12 atime
      integer :: timevals(8)
      real :: firsttm
      logical ifcpu,ifdisp,iffdisp
      logical :: first=TRUE
      
      data tprinted/0./
      save tprinted,first,firsttm
      
      real :: tlap
!
      iffdisp = ifdisp
      if(prtscrn .eq. FALSE) then
        iffdisp=FALSE
      endif
      
      if(first) then
        call date_and_time(values=timevals)
        firsttm=chglong(timevals)
        first=FALSE
      endif
!
      if(ifcpu) then
        call date_and_time(values=timevals)
        tm=chglong(timevals)-firsttm
        tm=tm-mod(tm,0.001)
        itm=tm
        mss=mod(tm,1.0)*1000
        ihh=itm/3600
        itm=itm-ihh*3600
        imm=itm/60
        itm=itm-imm*60
        iss=itm
        
        write(atime,'(2(i2.2,":"),i2.2,".",i3.3)') ihh,imm,iss,mss
        if(iffdisp) print 600, atime,trim(amesg)
        write(io8,600) atime,trim(amesg)
      else
        if(iffdisp) print 601,trim(amesg)
        write(io8,601) trim(amesg)
      endif
      
  600 format(a,1x,a)
  601 format(a)
  
    end subroutine
