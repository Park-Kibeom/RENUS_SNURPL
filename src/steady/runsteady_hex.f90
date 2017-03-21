! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
   subroutine runsteady_hex(flag2g,                &
                           noutbegmg,  ninbegmg,  &
                           noutbeg2g,  ninbeg2g,  &
                           epsl2, erreig,     errl2)

#define mg_correction    ! 2014_03_12 . scb
!#define cmfdmg_test   ! 2014_02_24 . scb
!#define sp3_dbg3 ! 2014_02_24 . scb

      use const
      use sfam
      use timer
      use itrinfo
      use sfam_cntl,  only : nmaxcmfd2g                             
      use geom,       only : ng
      use xsec,       only : colxs, colxs_sp3
      use tpen,       only : drivetpen
      use tpen_sp3,   only : drivetpen_sp3
      use geomhex,    only : ifhexfdm, ifhexsp3
      use cmfd2g,     only : drivecmfd2g
      use cmfdmg,     only : drivecmfdmg   ! 2014_02_24 . scb
      use cmfdhex2g_sp3,  only : upddtilhex2g_sp3, setlshex2g_sp3, ilufachex_sp3
      use fluxyoyo,   only : expphi      ! 2014_02_24 . scb
      use sfam_cntl,  only : ifcmfd    ! 2014_03_20 . scb
      use cmfd2g,     only : ifcmfd2g  ! 2014_12_05 . scb      
       
      implicit none

      logical                 :: flag2g
      integer                 :: noutbegmg,  ninbegmg
      integer                 :: noutbeg2g,  ninbeg2g
      real                    :: epsl2,erreig,errl2
      	      
      real*8 :: tinitscb=0.d0, tendscb=0.d0
      
      integer,save :: first=.true. ! 2013_04_25 . scb for debugging
      
      if(noutbeg2g.ne.0) then
        if(ifhexsp3) then
          call drivetpen_sp3(FALSE) 
        else
          call drivetpen(FALSE)   ! 2013_10_11 . scb 
        endif
        if(.not.ifcmfd) return  ! 2014_03_20 . scb
      endif
      
      call timeron()

! 2014_02_24 . scb      
#ifdef cmfdmg_test
      call drivecmfdmg(FALSE,TRUE,nmaxcmfdmg,noutbegmg,ninbegmg,eigv,reigv,phif,epsl2,erreig,errl2)  
! added end      
#else
  ! 2014_02_24 . scb  for dbg
#ifdef sp3_dbg3
      if(noutbeg2g.ge.nmaxcmfd2g) then
         continue
         goto 224
      endif
#endif      
            
      if(ng.ne.2) then
         call colxs(fphi)
         if(ifhexsp3) call colxs_sp3(fphi) 
         !call upddhathex2g_sp3   ! 2014_03_10 . scb
      endif
      if(ifhexsp3) then 
         call upddtilhex2g_sp3
         call setlshex2g_sp3(FALSE,1)
      else
         call upddtilhex       
         call setlshex2g(FALSE,1)
      endif
      call faciluhex

      ncmfd2g = ncmfd2g - noutbeg2g  
! 2012_12_28 . scb      
      !if(ifhexsp3) then 
      !  !call drivecmfd2g_sp3(FALSE,TRUE,nmaxcmfd2g,noutbeg2g,ninbeg2g,eigv,reigv,phi,epsl2,erreig,errl2) 
      !  call drivecmfd2g_sp3(FALSE,TRUE,nmaxcmfd2g,noutbeg2g,ninbeg2g,eigv,reigv,epsl2,erreig,errl2) 
      !else  
      !nmaxcmfd2g = 10  ! 2014_02_13 . scb
      
      if(ifcmfd2g .or. ng.eq.2) call drivecmfd2g(FALSE,TRUE,nmaxcmfd2g,noutbeg2g,ninbeg2g,eigv,reigv,phi,epsl2,erreig,errl2)  
        
        !if(first) call drivecmfd2g(FALSE,TRUE,nmaxcmfd2g,noutbeg2g,ninbeg2g,eigv,reigv,phi,epsl2,erreig,errl2)  
        !first=.false.
        !noutbeg2g=noutbeg2g+1
      !endif  
! added end    
      ncmfd2g = ncmfd2g + noutbeg2g

      call timeroff(tcmfd2g) 
      
! 2014_02_24 . scb
#ifdef mg_correction
      if(ng.ne.2) call expphi(phi,fphi,phif) 
#endif      
! added end

224   continue        
#endif 

!#ifdef DEBUG
!      if(noutbeg2g.ne.0) then
!         if(ifhexsp3) then
!            call tpenbc_sp3(false)    ! 2012_08_14 . scb
!         else
!            call tpenbc(false)    ! 2012_08_14 . scb
!         endif
!      endif
!#endif

      			 
   end subroutine
