module chebyacc
    use const
    implicit none
    
    logical                 ::  chebyon             &
                              , chebynew

    integer                 ::  ncheby  = 10        &
                              , iorder              &
                              , iorderm1            &
                              , icheby0 = 3         &
                              , ioutcb0
    
    real                    ::  domrmin = 0.5       &
                              , domrmax = 1.0       &
                              , sigma   = 0.0       &
                              , sigbar  = 0.0       &
                              , erfcb
                                
    real,pointer            ::  phiafd(:,:,:)       &
                              , phiafdd(:,:,:)      &
                              , psiafd(:,:,:)       &
                              , psiafdd(:,:,:)

    interface runcheby
        module procedure runcheby1
        module procedure runcheby2
    end interface
    
    contains
    
    subroutine malloccheby
        use allocs
        use geom, only      : ng, nxy, nz
        call dmalloc(phiafd , ng, nxy, nz)
        call dmalloc(phiafdd, ng, nxy, nz)
        call dmalloc(psiafd , ng, nxy, nz)
        call dmalloc(psiafdd, ng, nxy, nz)

        return
    end subroutine

    subroutine runcheby1(iout, domr, phif, psif)
        use geom, only      :  ng, nxy, ny, nz      &
                             , kfbeg, kfend         &
                             , nxsf, nxef           &
                             , nodel
        implicit none      
        
        integer             :: iout
        real                :: domr
        real,pointer        :: phif(:,:,:)         &
                             , psif(:,:)
      
        real                :: alphacb              &
                             , betacb               &
                             , theta
                             
        integer             :: m, l, k, i, j

        if(chebyon) then
            if(chebynew) then
                sigma=sigbar
!                if(iout.le.6) sigma=min(sigma,0.900)
!                if(iout.le.9) sigma=min(sigma,0.950)
!                if(iout.le.12) sigma=min(sigma,0.985)
                alphacb   = 2/(2-sigma)
                betacb    = 0.0
                chebynew  = FALSE
                iorder    = 1
                erfcb     = 1
            else
                iorderm1  = iorder
                iorder    = iorder + 1
                theta     = acosh(2/sigma - 1)
                alphacb   = 4*cosh(iorderm1*theta)/(sigma*cosh(iorder*theta))
                betacb    = (1-sigma*0.5)*alphacb - 1
            endif

!           flux extrapolation
            do k=1,nz
            do l=1,nxy
            do m=1,ng
                phif(m,l,k) = phiafd(m,l,k)                            &
                            + alphacb*(phif(m,l,k)  - phiafd(m,l,k))   &
                            + betacb *(phiafd(m,l,k) - phiafdd(m,l,k))
            enddo
            enddo
            enddo

!           fission source extrapolation
            do k  = kfbeg, kfend
            do j  = 1, ny
            do i  = nxsf(j),nxef(j)
                l=nodel(i,j)
                
                psif(l,k) = psiafd(1,l,k)                        &
                         +  alphacb*(psif(l,k)-psiafd(1,l,k))    &
                         +  betacb*(psiafd(1,l,k)-psiafdd(1,l,k))
            enddo
            enddo
            enddo
        endif
        do  k = 1, nz
        do  l = 1, nxy
            phiafdd(:,l,k) = phiafd(:,l,k)
            phiafd(:,l,k)  = phif(:,l,k)
            psiafdd(1,l,k) = psiafd(1,l,k)
            psiafd(1,l,k)  = psif(l,k)
        enddo 
        enddo
        
        call turnchebyon(iout, domr)
    end subroutine
    
    subroutine runcheby2(iout, domr, phiaf, psiaf)
        use geom, only      :  ng, nxy, ny, nz      &
                             , kfbeg, kfend         &
                             , nxsf, nxef           &
                             , nodel
        implicit none      
        
        integer             :: iout
        real                :: domr
        real,pointer        :: phiaf(:,:,:)         &
                             , psiaf(:,:,:)
      
        real                :: alphacb              &
                             , betacb               &
                             , theta
                             
        integer             :: m, l, k, i, j

        if(chebyon) then
            if(chebynew) then
                sigma=sigbar
!                if(iout.le.6) sigma=min(sigma,0.900)
!                if(iout.le.9) sigma=min(sigma,0.950)
!                if(iout.le.12) sigma=min(sigma,0.985)
                alphacb   = 2/(2-sigma)
                betacb    = 0.0
                chebynew  = FALSE
                iorder    = 1
                erfcb     = 1
            else
                iorderm1  = iorder
                iorder    = iorder + 1
                theta     = acosh(2/sigma - 1)
                alphacb   = 4*cosh(iorderm1*theta)/(sigma*cosh(iorder*theta))
                betacb    = (1-sigma*0.5)*alphacb - 1
            endif

!           flux extrapolation
            do k=1,nz
            do l=1,nxy
            do m=1,ng
                phiaf(m,l,k) = phiafd(m,l,k)                            &
                            + alphacb*(phiaf(m,l,k)  - phiafd(m,l,k))   &
                            + betacb *(phiafd(m,l,k) - phiafdd(m,l,k))
            enddo
            enddo
            enddo

!           fission source extrapolation
            do k  = kfbeg, kfend
            do j  = 1, ny
            do i  = nxsf(j),nxef(j)
                l=nodel(i,j)
                
                psiaf(:,l,k) = psiafd(:,l,k)                        &
                         +  alphacb*(psiaf(:,l,k)-psiafd(:,l,k))    &
                         +  betacb*(psiafd(:,l,k)-psiafdd(:,l,k))
            enddo
            enddo
            enddo
        endif
        do  k = 1, nz
        do  l = 1, nxy
            phiafdd(:,l,k) = phiafd(:,l,k)
            phiafd(:,l,k)  = phiaf(:,l,k)
            psiafdd(:,l,k) = psiafd(:,l,k)
            psiafd(:,l,k)  = psiaf(:,l,k)
        enddo 
        enddo
        
        call turnchebyon(iout, domr)
    end subroutine
    
    subroutine turnchebyon(iout, domr)
        implicit none
        integer             ::  iout
        
        real                ::  domr          &
                              , theta         &
                              , erfcb0        &
                              , erfcbr

        if(       .not.chebyon                &
            .and. domr.gt.domrmin             &
            .and. (iout-ioutcb0).ge.icheby0   &
          ) then
          
            chebyon   = TRUE
            chebynew  = TRUE

        endif
        
        sigbar  = domr
        
        if(chebyon .and. iorder.gt.1) then
            erfcb   =   erfcb*domr
            theta   =   acosh(2/sigma-1)

!           theoretical error reduction factor
            erfcb0  =   1/cosh(iorderm1*theta)

!           relative             
            erfcbr  =   erfcb/erfcb0
            
            if(erfcbr.lt.1.0) then
                sigbar=sigma*(1+cos(acos(erfcbr)/iorderm1))*0.5
            else
                sigbar=sigma*(1+cosh(acosh(erfcbr)/iorderm1))*0.5
            endif
        endif

!       reset polynomial if no further improvement is expected
        if(iorder.ge.3) then
            if(iorder.ge.ncheby .or. erfcbr.gt.1.0) chebynew = .TRUE.
            
            if(sigbar.ge.domrmax) then
                chebynew  = TRUE
                chebyon   = FALSE
                ioutcb0   = iout
            endif
        endif
    end subroutine
end module