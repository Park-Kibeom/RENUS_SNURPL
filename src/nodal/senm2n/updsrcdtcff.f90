! 2014_09_26 . scb
!#define printratio  
  subroutine updsrcdtcff(idir,l,k,m,srccff)    
    use allocs   ! 2014_11_06 . scb
    use const
    use senm2n
    use bdf,        only : bdfcoef,bdfordersave,dtcff
    use tranxsec,   only : rvelof
    use xsec,       only : xstrf
    use geom,       only : hmesh
    use geom,       only : nxy,nz   ! 2014_11_06 . scb
    
    integer                 :: idir,l,k,m
    real                    :: srccff(0:2)
    
    real                    :: srcdtcff(0:2), coef
    real                    :: dtcff0(0:2)
    integer :: nbdf
    
! 2014_09_27 . scb
    real,save :: time=0.d0
    logical,save :: first=.true.
    
    real :: ratio(0:2), dum(6)
    real,save :: maxratio(0:2)=0.d0
    
    integer :: iodbg
    character*256 :: oneline
    
    real :: timesrcdtdbg  ! 2014_09_27 . scb
    common / srcdtdbg / timesrcdtdbg
    
    integer,save :: istepdbg=0
! added end      

! 2014_11_06 . scb
    real,pointer :: rat(:,:,:,:), ratmax(:,:), ratavg(:,:)
    real :: dum2
    common / srcdtdbg2 / rat, ratmax, ratavg
    
    integer :: nt
    
#ifdef printratio    
    if(first) then
      call dmalloc(rat,ng,nxy,nz,3)
      call dmalloc(ratmax,ng,3)
      call dmalloc(ratavg,ng,3)
      
      first=.false.
      open(unit=1106,file='srcdt_ratio',status='unknown')
      
      write(1106,'(a)')  "   Time(sec)  ratmax(1:ng,1:ndir)    ratavg(1:ng,1:ndir)"
      
      time=timesrcdtdbg
      
    elseif(time.ne.timesrcdtdbg) then      
            
      nt=nxy*nz
      do ir=1,3
        do ig=1,ng
          ratavg(ig,ir)=(ratavg(ig,ir)/nt)**0.5d0
        enddo
      enddo
      
      write(1106,'(f12.4,12es13.4)') time, &
        ratmax(1,1), ratmax(1,2), ratmax(1,3), ratmax(2,1), ratmax(2,2), ratmax(2,3), &
        ratavg(1,1), ratavg(1,2), ratavg(1,3), ratavg(2,1), ratavg(2,2), ratavg(2,3) 
            
      time=timesrcdtdbg              
    endif
    
    if(l.eq.1 .and. k.eq.1) then
      ratmax(m,idir)=0.d0
      ratavg(m,idir)=0.d0
    endif
#endif    
    
! added end
    
    nbdf=bdfordersave
    rh=1.d0/hmesh(idir,l,k)
    dtcff0(0:2) = dtcff(:,m,l,k,idir,0)    
    
    dtcff0=0.d0
    do icff=0,2
      do iorder=1,nbdf
        dtcff0(icff) = dtcff0(icff) + bdfcoef(iorder)*dtcff(icff,m,l,k,idir,iorder)
      enddo
      dtcff0(icff) = rvelof(m,l,k)*dtcff0(icff)
    enddo
    
    dum2=4.d0/(3.d0*rh*rh*xstrf(m,l,k))*(3.d0*phicff(2,m,l,k,idir)+10.d0*phicff(4,m,l,k,idir))  ! 2014_11_06 . scb
        
    dtcff0(0) = dtcff0(0) + 2.d0/3.d0*rh*(3.d0*phicff(2,m,l,k,idir)+10.d0*phicff(4,m,l,k,idir))
    dtcff0(1) = dtcff0(1) + 10.d0*rh*phicff(3,m,l,k,idir)
    dtcff0(2) = dtcff0(2) + 70.d0/3.d0*rh*phicff(4,m,l,k,idir)
    
    coef=-1.d0/(xstrf(m,l,k) + bdfcoef(0)*rvelof(m,l,k))
    do icff=0,2
      dtcff0(icff) = coef*dtcff0(icff)
    enddo
    dtcff(:,m,l,k,idir,0) = dtcff0(0:2)  
    
    coef=2.d0*rvelof(m,l,k)*rh/xstrf(m,l,k)
    srcdtcff=0.d0
    do icff=0,2
      do iorder=0,nbdf
        srcdtcff(icff) = srcdtcff(icff) +  bdfcoef(iorder)*dtcff(icff,m,l,k,idir,iorder)
      enddo
      srcdtcff(icff) = coef*srcdtcff(icff)
    enddo
    
! 2014_09_27 . scb        
    !do icff=0,2
    !  ratio(icff)=abs(srcdtcff(icff))/abs(srccff(icff))
    !enddo    
    !
    !if(first) then
    !  first=.false.
    !  open(unit=9270,file='srcdt0_dbg',status='unknown')
    !  open(unit=9271,file='srcdt1_dbg',status='unknown')
    !  open(unit=9272,file='srcdt2_dbg',status='unknown')
    !  
    !  write(9270,'(a)')  "   Time(sec)  srcdtcff(0)    srccff(0)  maxratio(0)"
    !  write(9271,'(a)')  "   Time(sec)  srcdtcff(1)    srccff(1)  maxratio(1)"
    !  write(9272,'(a)')  "   Time(sec)  srcdtcff(2)    srccff(2)  maxratio(2)"
    !endif    
    !
    !if(time.ne.timesrcdtdbg) then
    !  time=timesrcdtdbg
    !  do icff=0,2
    !    iodbg=9270+icff
    !    write(iodbg,'(f12.4,3es13.4)') timesrcdtdbg, srcdtcff(icff), srccff(icff), ratio(icff)
    !    maxratio(icff)=ratio(icff)
    !  enddo      
    !else
    !  do icff=0,2
    !    iodbg=9270+icff
    !    if(maxratio(icff) .lt. ratio(icff)) then
    !      backspace(iodbg)
    !      write(iodbg,'(f12.4,3es13.4)') time, srcdtcff(icff), srccff(icff), ratio(icff)
    !    endif
    !  enddo                
    !endif    
! added end
    
    srccff(:)=srccff(:)+srcdtcff(:)
    
#ifdef printratio  
! 2014_11_06 . scb
    rat(m,l,k,idir)=1.d0/(abs(dum2)/abs(srcdtcff(0)) + 1.d0)
    if(ratmax(m,idir) .lt. rat(m,l,k,idir)) ratmax(m,idir) = rat(m,l,k,idir)
    ratavg(m,idir) = ratavg(m,idir) + rat(m,l,k,idir)*rat(m,l,k,idir)
! added end
#endif    
    
  end subroutine
  
