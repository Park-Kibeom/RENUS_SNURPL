module tranxsec
    use const
    use allocs
    implicit none
    real,pointer,dimension(:,:,:):: rvelof,     &   !(ng,nxy,nz)     multigroup velocity
                                    lmbdk,      &   !(nprec,nxy,nz)  decay constant, lamda_k
                                    betak           !(nprec,nxy,nz)  k-th delayed neutron fraction
    real,pointer,dimension(:,:,:):: xschid,     &
                                    xschifd          ! for transient
! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
    real,pointer,dimension(:)    :: decalpha,   &   ! DECAY HEAT
                                    deczeta         ! DECAY HEAT
! added end
    contains
    
    subroutine malloctranxsec(ng,       nxy,    nz,         &
                              lmbdk0,   betak0, xschifd0,   &
                              rvelof0)
        integer                 :: ng,nxy,nz
        real,   pointer         :: lmbdk0(:,:,:),        &
                                   betak0(:,:,:),        &
                                   xschifd0(:,:,:),      &
                                   rvelof0(:,:,:)
        lmbdk   =>  lmbdk0
        betak   =>  betak0
        rvelof  =>  rvelof0
        xschifd =>  xschifd0
        
        if(ng.ne.ng2) then
            call dmalloc(xschid,ng2,nxy,nz)
        else
            xschid   => xschifd
        endif
    end subroutine
end module