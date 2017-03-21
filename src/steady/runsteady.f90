subroutine runsteady(flag2g,                &
                     noutbegmg,  ninbegmg,  &
                     noutbeg2g,  ninbeg2g,  &
                     epsl2, erreig,     errl2)
    use const
    use sfam
    use timer
    use itrinfo
    use sfam_cntl,  only : nmaxcmfdmg, nmaxcmfd2g, &
                            ifrect       ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
    use geom,       only : ng, isymloc
    use geom,       only : ifrecsp3   ! 2013_07_05 . scb
    use xsec,       only : xschif,colxs

    use cmfdmg,     only : setlsmg,     drivecmfdmg

    use cmfd2g,     only : upddtil2g,   setls2g,        &
                           upddhat2g,   drivecmfd2g
                           
    use fluxyoyo,   only : colphi,      expphi
                           
    use bicgmg,     only : facilumg
    use bicg2g,     only : facilu2g
    
    use nodal,      only : drivenodal
    
    use cmfd2g,     only : ifcmfd2g  ! 2014_12_05 . scb
        
    implicit none
    
    logical                 :: flag2g
    integer                 :: noutbegmg,  ninbegmg
    integer                 :: noutbeg2g,  ninbeg2g
    real                    :: epsl2,erreig,errl2
    
    integer                 :: i,j,l
     
    
    !if(noutbeg2g.ne.0) then
    if(noutbeg2g+noutbegmg.ne.0) then
! 2013_07_05 . scb
      if(ifrecsp3) then
        call drivesenm_sp3(FALSE)
      else
! added end        
        call drivenodal(FALSE,reigv,phif,psi,jnet,phisfc)
      endif      
    endif
    
!    call timeron()       ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
    !if(.not.flag2g) then   ! 2013_07_22 . scb
    if(.not.flag2g .or. ifrecsp3) then
        call timeron()       ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
        if(updls) then
            call setlsmg(FALSE)
            call facilumg
        endif
        ncmfdmg = ncmfdmg - noutbegmg
        call drivecmfdmg(FALSE,TRUE,nmaxcmfdmg,noutbegmg,ninbegmg,eigv,reigv,xschif,phif,epsl2,erreig,errl2)
        ncmfdmg = ncmfdmg + noutbegmg
        
        if(ifrecsp3) goto 300
        
        goto 300   ! 2015_01_23 . scb added

        call colphi(phif,phi,fphi)
        call colxs(fphi)
        call upddtil2g              
        call upddhat2g(phif,phi) 
        updls = TRUE !to update 2G linear system 
        call timeroff(tcmfd)      ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
    endif
    
    call timeron()       ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
    if(updls) then  
        call setls2g(FALSE)
        call facilu2g ! initialize ilu.      
    endif

    ncmfd2g = ncmfd2g - noutbeg2g
    if(ifcmfd2g .or. flag2g) call drivecmfd2g(FALSE,TRUE,nmaxcmfd2g,noutbeg2g,ninbeg2g,eigv,reigv,phi,epsl2,erreig,errl2)
    ncmfd2g = ncmfd2g + noutbeg2g
    call timeroff(tcmfd2g)      ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB

    if(.not.flag2g) then
        call timeron()       ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
!to achieve stability
        call expphi(phi,fphi,phif)
        ncmfdmg = ncmfdmg - noutbegmg
        call drivecmfdmg(FALSE,TRUE,1,noutbegmg,ninbegmg,eigv,reigv,xschif,phif,epsl2,erreig,errl2) 
        ncmfdmg = ncmfdmg + noutbegmg
        call timeroff(tcmfd)      ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB
    endif
!    call timeroff(tcmfd)        ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB 

300 continue
    
end subroutine