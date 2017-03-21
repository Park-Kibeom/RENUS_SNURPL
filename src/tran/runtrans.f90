subroutine runtrans( flag2g,flagnodal,      &
                     noutbegmg,  ninbegmg,  &
                     noutbeg2g,  ninbeg2g,  &
                     epsl2tr,errl2)
    use const
    use sfam
    use itrinfo
    use timer
    use sfam_cntl,  only : nmaxcmfdmg, nmaxcmfd2g, &
                            nodalupd, cmfdupd, ifrect    ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
    use tran,       only : chibetapf, rveloftm, &
                            deltm    ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
    use geom,       only : ng, nxy,nz, volnode
    use geomhex,    only : ifhexsp3   ! 2012_12_06 . scb
    use xsec,       only : xschif,colxs

    use cmfdmg,     only : setlsmg,     drivecmfdmg, diagf

    use cmfd2g,     only : upddtil2g,   setls2g,        &
                           upddhat2g,   drivecmfd2g
                           
    use fluxyoyo,   only : colphi,      expphi
                               
    use bicgmg,     only : facilumg
    use bicg2g,     only : facilu2g
    
    use nodal,      only : drivenodal
    
    use tpen,       only : drivetpen    ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
    use bdf   ! 2012_09_26 . scb
    use geom,       only : ifrecsp3  ! 2013_08_13 . scb
    use cmfd2g,     only : ifcmfd2g  ! 2014_12_05 . scb    
        
    implicit none
    
    logical                 :: flag2g, flagnodal
    integer                 :: noutbegmg,  ninbegmg
    integer                 :: noutbeg2g,  ninbeg2g
    real                    :: epsl2tr,errl2
    
    real                    :: erreig,errl2temp
    integer                 :: m,l,k
    
! added in ARTOS ver. 0.2 . 2012_07_24 by SCB
    real,save :: time=0.
    
    time=time+deltm    
        
! 2012_09_26 . scb
!    if(flagnodal) then
!        call drivenodal(TRUE,reigv,phif,psi,jnet,phisfc)
!    endif

!    if(flagnodal .and. nodalupd) then
    if(flagnodal .and. nodalupd .and. flagmain) then
! 2013_08_13 . scb      
      if(ifrect) then
        if(ifrecsp3) then
          call drivesenm_sp3(TRUE)
        else      
! added end          
          call drivenodal(TRUE,reigv,phif,psi,jnet,phisfc)
        endif        
      else
! 2012_12_06 . scb
        if(ifhexsp3) then
          call drivetpen_sp3(TRUE) 
        else
! added end        
          call drivetpen(true)   ! 2013_10_11 . scb 
        endif
      endif
    endif
    
    !if(.not.ifrect) goto 150  ! 2013_10_11 . scb
! added end

!    call timeron()     ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
    if(.not.flag2g) then
      if(ifrect) then   ! 2013_10_11 . scb
        call timeron()     ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
        call setlsmg(TRUE)
        do k=1,nz
        do l=1,nxy
        do m=1,ng
            diagf(l,k,m) = diagf(l,k,m)+rveloftm(m,l,k)*volnode(l,k)
        enddo
        enddo
        enddo
        
        call facilumg
      
! 2012_09_26 . scb      
!        ncmfdmg = ncmfdmg - iterinnermg
        if(flagmain)  ncmfdmg = ncmfdmg - iterinnermg
!        call drivecmfdmg(TRUE,TRUE,nmaxcmfdmg,noutbegmg,ninbegmg,eigv,reigv,chibetapf,phif,epsl2tr,erreig,errl2)
        call drivecmfdmg(TRUE,TRUE,nmaxcmfdmg,iterinnermg,ninbegmg,eigv,reigv,chibetapf,phif,epsl2tr,erreig,errl2)
!        ncmfdmg = ncmfdmg + iterinnermg
        if(flagmain)  ncmfdmg = ncmfdmg + iterinnermg 
! added end        
        
        if(errl2.le.epsl2tr) goto 100
        
        call colphi(phif,phi,fphi)
        call colxs(fphi)
        call upddtil2g             
        call upddhat2g(phif,phi) 
        
        call timeroff(tcmfd)    ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
      else
! 2013_10_11 . scb        
        call colxs(fphi)
        if(ifhexsp3) call colxs_sp3(fphi)
! added end        
      endif
        
    endif

!150 continue     ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
    
    call timeron()     ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   

! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
!    call setls2g(TRUE)
    if(ifrect) then
      call setls2g(TRUE)
    else
! 2012_12_06 . scb      
      if(ifhexsp3) then 
        call upddtilhex2g_sp3
        call setlshex2g_sp3(TRUE,1)
! added end         
      else
        call upddtilhex
        call setlshex2g(TRUE,1)   ! 2013_10_11 . scb must be true
      endif
    endif
! added end
    
! 2012_09_26 . scb    
!    ncmfd2g = ncmfd2g - iterinner2g
    if(flagmain)  ncmfd2g = ncmfd2g - iterinner2g
    if(ifcmfd2g .or. flag2g)  call drivecmfd2g(TRUE,TRUE,nmaxcmfd2g,iterinner2g,ninbeg2g,eigv,reigv,phi,epsl2tr,erreig,errl2)
!    call drivecmfd2g(TRUE,TRUE,nmaxcmfd2g,noutbeg2g,ninbeg2g,eigv,reigv,phi,epsl2tr,erreig,errl2)
    !if(ifhexsp3) then
    !  !call drivecmfd2g_sp3(TRUE,TRUE,nmaxcmfd2g,iterinner2g,ninbeg2g,eigv,reigv,phi,epsl2tr,erreig,errl2)   ! 2012_12_06 . scb
    !  !call drivecmfd2g_sp3(TRUE,TRUE,nmaxcmfd2g,iterinner2g,ninbeg2g,eigv,reigv,epsl2tr,erreig,errl2) 
    !  call drivecmfd2g(TRUE,TRUE,nmaxcmfd2g,iterinner2g,ninbeg2g,eigv,reigv,phi,epsl2tr,erreig,errl2)  ! 2012_12_06 . scb
    !else
    !  call drivecmfd2g(TRUE,TRUE,nmaxcmfd2g,iterinner2g,ninbeg2g,eigv,reigv,phi,epsl2tr,erreig,errl2)
    !endif
!    ncmfd2g = ncmfd2g + iterinner2g
    if(flagmain)  ncmfd2g = ncmfd2g + iterinner2g
! added end    
    
    call timeroff(tcmfd2g)    ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
    
    if(.not.ifrect) return    ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
    
    if(flag2g .and. errl2.le.epsl2tr) goto 100

    if(.not.flag2g) then
!to achieve stability
        call timeron()     ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
        
        call expphi(phi,fphi,phif)

! 2012_09_26 . scb      
!        ncmfdmg = ncmfdmg - iterinnermg
        if(flagmain)  ncmfdmg = ncmfdmg - iterinnermg
!        call drivecmfdmg(TRUE,TRUE,1,noutbegmg,ninbegmg,eigv,reigv,chibetapf,phif,epsl2tr,erreig,errl2)
        call drivecmfdmg(TRUE,TRUE,1,iterinnermg,ninbegmg,eigv,reigv,chibetapf,phif,epsl2tr,erreig,errl2)
!        ncmfdmg = ncmfdmg + iterinnermg
        if(flagmain)  ncmfdmg = ncmfdmg + iterinnermg 
! added end    
        
        call timeroff(tcmfd)    ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
    endif
    
100 continue

    if(.not.flag2g) then
        call timeron()     ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
        do k=1,nz
        do l=1,nxy
        do m=1,ng
            diagf(l,k,m) = diagf(l,k,m)-rveloftm(m,l,k)*volnode(l,k)
        enddo
        enddo
        enddo    
        call timeroff(tcmfd)    ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
    endif
!    call timeroff(tcmfd)        ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
    
    return
end subroutine