!    program ARTOS
!    
!      use param
!      use timer
!      use allocs
!      use geom,     only : setgeom
!      use xsec,     only : setxsec
!      use sfam,     only : mallocsfam, initsfam, eigv     
!      use itrinfo,  only : tnodal,tcmfd,tother,nnodal,ncmfd2g,ncmfdmg
!            
!      include 'global.h'
!      include 'itrcntl.h'
!      include 'ffdm.h'
!      include 'xsec.h'
!      include 'geom.h'
!      include 'nodal.h'
!      include 'pinpr.h'
!      include 'times.h'
!      include 'cntl.h'
!      include 'ff.h'
!      include 'srchppm.h'
!      include 'thcntl.inc'
!      integer :: ncmfdmgtot=0,ncmfd2gtot=0,ninmgtot=0,nin2gtot=0
!
!      call timeron()
!
!! initialize
!      call timeron()
!      call preproc
!      call readinput
!      call init   
!      call timeroff(tinit)
!      
!      call setgeom(ng,nx,ny,nz,nxy,ndir,symopt,isymang, &
!                   isymloc,nxs,nxe,nys,nye,nxsf,nxef, &
!                   kfbeg,kfend,nodel,hmesh,volnode,volcore,volfuel, &
!                   (/albxl(1),albxr(1),albyl(1),albyr(1),albzb(1),albzt(1)/) )
!
!      call setxsec(ng,mgb(2),nxy,nz,xstf,xsaf,xsdf,xsnff,xskpf,xschif,  &
!                   xbeta,xsadf,xscdf,xssf,xssfs,xssfe)
!     
!      call dmalloc(curil,ndirmax,nxy,nz,ng)
!      call dmalloc(curir,ndirmax,nxy,nz,ng)
!      call dmalloc(curol,ndirmax,nxy,nz,ng)
!      call dmalloc(curor,ndirmax,nxy,nz,ng)
!      
!      call mallocsfam(ng,nxy,nz,phif,phisfc,jnet,curol,curor,curil,curir)
!      call initsfam(TRUE)
!      
!! standard code
!      call runss(ncmfdmgtot,ninmgtot,ncmfd2gtot,nin2gtot)
!      call normalize
!
!      if(printff) call writeff
!            
!! pin power calculation
!      if(pinpower) call driveppr(FALSE)
!
!! transient calculation
!      if(transient) call runtr
!
!      call timeroff(ttotal)      
!            
!! write results to file	
!      call editout
!
!      write(mesg, '("    ")')
!      call message(FALSE,true,mesg) 
!      write(mesg, '("=======================================    Result    =======================================")')
!      call message(FALSE,true,mesg) 
!      write(mesg,'(a18,f10.5)')       'K-EFF           : ', eigv
!      call message(FALSE,true,mesg)
!      if(srchppm .or. fdbk) then
!          write(mesg,'(a18,f10.3)')       'BORON           : ', ppm
!          call message(FALSE,true,mesg)
!      endif
!
!      write(mesg, '("    ")')
!      write(mesg, '("=================================    Performance Summary    ================================")')
!      call message(FALSE,true,mesg) 
!      write(mesg,'(a18,f10.3)')       'Initial    Time : ',tinit
!      call message(FALSE,true,mesg)
!      write(mesg,'(a18,f10.3)')       'XSEC       Time : ',txsec
!      call message(FALSE,true,mesg)      
!      write(mesg,'(a18,f10.3,a9,i4)') 'Nodal-Calc Time : ',tnodal, '##: ',nnodal
!      call message(FALSE,true,mesg)
!      write(mesg,'(a18,f10.3,a9,i4,a9,i4)')   'CMFD-Calc  Time : ',tcmfd, '2G: ',ncmfd2g,'MG: ',ncmfdmg
!      call message(FALSE,true,mesg)
!      write(mesg,'(a18,f10.3)')       'T/H        Time : ',tth
!      call message(FALSE,true,mesg)
!      write(mesg,'(a18,f10.3)')       'PPR        Time : ',tppr
!      call message(FALSE,true,mesg)      
!      write(mesg,'(a18,f10.3)')       'OTHER      Time : ',tother
!      call message(FALSE,true,mesg)      
!      write(mesg,'(a18,f10.3)')       'Total-CPU  Time : ',ttotal
!      call message(FALSE,true,mesg)
!
!      negpnt=0
!      do m=1, ng
!        do k=1,nz
!          do l=1,nxy
!            if(phif(m,l,k) .lt. 0) negpnt=negpnt+1
!          enddo
!        enddo
!      enddo
!      
!      if(negpnt .gt. 0) print *, negpnt," NEGATIVE POINTS / ",nxy*nz*ng," POINTS"
!      
!      continue
!   
!    end program
