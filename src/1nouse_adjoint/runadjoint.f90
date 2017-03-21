!subroutine runadjoint(phiaf0, epsl2)
!    use const
!    use adjoint
!    use geom,       only : ng, nxy, nz
!    use sfam_cntl,  only : nmaxcmfd2g
!    use sfam,       only : resetsfam, phif, phi
!    use cmfdmg,     only : dhatrf,dhatzf
!    use cmfd2g,     only : dhatr,dhatz
!
!    implicit none
!
!    real, pointer           :: phiaf0(:,:,:)
!    real                    :: epsl2
!    
!    integer                 :: noutbegmg = 0 , &
!                               ninbegmg  = 0 , &
!                               noutbeg2g = 0 , &
!                               ninbeg2g  = 0
!
!    integer                 :: k, l, iout
!    character               :: mesg*120    
!    real, pointer           :: phif0(:,:,:)    
!    
!    dhatrf = 0.0
!    dhatzf = 0.0
!    dhatr  = 0.0
!    dhatz  = 0.0
!    
!    do k=1,nz
!    do l=1,nxy
!        phiaf(:,l,k) = phiaf0(:,l,k)
!    enddo
!    enddo
!    
!    phif0 => phif
!    phif  => phiaf
!    phi   => phiaf
!
!    write(mesg,'(a)'), '* ADJOINT ITERATION : BEGIN'
!    call message(true,true,mesg)            
!
!    if(ng.eq.ng2) then
!        call transadj
!        call resetsfam()
!        do iout = 1, noutmax
!            call runsteady( ng.eq.ng2 , iout.gt.nmaxcmfd2g   ,  &
!                            noutbegmg , ninbegmg ,  &
!                            noutbeg2g , ninbeg2g ,  &
!                            epsl2 , erreig , errl2)
!
!            if(errl2 .lt. epsl2) exit
!        enddo
!        call transadj        
!    else
!        call setlsmg(FALSE)
!        call transmg
!        call facilumg
!        call caladjmg(noutmax)
!    endif
!    write(mesg,'(a)'), '* ADJOINT ITERATION : END'
!    call message(true,true,mesg)            
!    
!    phif  => phif0
!    phi   => phif0
!
!    dhatrf = 0.0
!    dhatzf = 0.0
!    dhatr  = 0.0
!    dhatz  = 0.0
!end subroutine