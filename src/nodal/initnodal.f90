!subroutine initnodal
subroutine initnodal(ifrect)  ! 2012_10_12 . scb
    use const
    use nodal
    use sanm2n, only : initsanm2n
    use senm2n, only : initsenm2n
    use anm2n,  only : initanm2n
    use sfam_cntl, only : nodalkrnl
    implicit none

    logical, save           :: first=TRUE
    integer                 :: ndefaultswp=6
    character*80            :: mesg
    logical                 :: ifrect  ! 2012_10_12 . scb

    if(first.eq.FALSE) return
    is1n=FALSE
    
    nmaxswp = min(ng+2,ndefaultswp)
    nlupd=0;nlswp=0
    
    call mallocnodal

    select case(trim(nodalkrnl))
    case('SANM2N')
        call initsanm2n
    case('SENM2N')
        call initsenm2n
    case('ANM2N')
        call initanm2n
    end select

    write(mesg, '("Selected Nodal Kernel :",a)') nodalkrnl
    
    if(.not.ifrect)  write(mesg, '("Selected Nodal Kernel : TPEN-NEM")')  ! 2012_10_12 . scb
    call message(true, true, mesg)
    
    first=FALSE
    
    return   
!    
end subroutine   

    
