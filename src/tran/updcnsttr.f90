subroutine updcnsttr(phi,phif,fphi,plevel)
!   determine initial parameters
    use const
    use geom,       only : ng,nxy,nz
    use xsec,       only : mgb,mge
    use tranxsec,   only : lmbdk
    use tran
    implicit none
    
    real,pointer            :: phi(:,:,:),          &
                               phif(:,:,:),         &
                               fphi(:,:,:)
    real                    :: plevel
    
    character*80            :: mesg
    integer                 :: m,m2,l,k,mp
    real                    :: rtdblr,rtdblrt

!   adaptive theta determination for multigroup
    rtdblr=0
    if(pleveld .ne. 0.) then
        rtdblr=log(plevel/pleveld)/log(2.) !previous doubling time
        
        if(ng .eq. 2 .or. thetak0.gt. 0.75 .or. rtdblr.gt.0) then
            thetak=thetak0
        else      
            if(rtdblr .gt. rtdblrd) then !power changing rate deccreases?
                if(rtdblr .gt. -0.1) then 
                    thetak=ONE
                else
                    thetak=HALF+RTHREE
                endif
            else
                if(rtdblr.gt.-0.2) then !5 times
                    thetak=HALF+RTHREE
                else
                    thetak=thetak0
                endif      
            endif
        endif
    else
        thetak=ONE
    endif

    thetatm=thetak*deltm
!    if(.not.exptheta) &   ! 2012_08_22 . scb for fixing of a bug in exponential transformation
        thetabar=(ONE-thetak)/thetak

    rtdblrd=rtdblr
    pleveld=plevel

    write(mesg,'(a,2f10.3,1p,e12.5,i5)') ' Doubling Time & Theta', thetak,pleveld
    call message(true,true,mesg)
    
    ! calculate constants for precursor calculation
    do k=1,nz
    do l=1,nxy
        cappa(:,l,k)=exp(-lmbdk(:,l,k)*deltm)
    enddo
    enddo                                                          

    ! calculate 1/(theta*velocity*deltm) and exponential term of theta method       
    if(ng.ne.ng2) then
        do k=1,nz
        do l=1,nxy
            do m=1,ng
                rveloftm(m,l,k)=rvelof(m,l,k)/thetatm
                if(exptheta) then
                    pinvf(m,l,k)=log(abs(phif(m,l,k)/phifd(m,l,k)))/deltmd
                    expof(m,l,k)=exp(pinvf(m,l,k)*deltm) ! exponential term of theta method with exponential transform
                    pinvf(m,l,k)=pinvf(m,l,k)*rvelof(m,l,k)
                endif  
            enddo

            do m2=1,ng2
                rvelo(m2,l,k)=0
                do m=mgb(m2),mge(m2)
                    rvelo(m2,l,k)=rvelo(m2,l,k)+rvelof(m,l,k)*fphi(m,l,k)
                enddo !m
            enddo !m2          
        enddo
        enddo   
    endif ! not ng2     

    do k=1,nz
    do l=1,nxy
    do m=1,ng2
        rvelotm(m,l,k)=rvelo(m,l,k)/thetatm
        if(exptheta) then
            pinv(m,l,k)=0
            if(phi(m,l,k).gt.0 .and. phid(m,l,k).gt.0 .and. phi(m,l,k)/phid(m,l,k).lt.10 ) then
                pinv(m,l,k)=log(abs(phi(m,l,k)/phid(m,l,k)))/deltmd
            endif
            expo(m,l,k)=exp(pinv(m,l,k)*deltm) ! exponential term of theta method with exponential transform
            pinv(m,l,k)=pinv(m,l,k)*rvelo(m,l,k)
        endif          
    enddo
    enddo
    enddo        

    return
end subroutine