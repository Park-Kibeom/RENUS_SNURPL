subroutine solbicg2g(r20, r2,phi)
    use const
    use bicg2g
    use cmfd2g,     only : axb2g
    use scalprods,  only : scalprod
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    use cmfdhex2g,  only : axbhex
    use sfam_cntl, only : ifrect
! added end
    implicit none
    
    real                    :: r20
    real                    :: r2
    real,pointer           :: phi(:,:,:)

! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
!    real                    :: crhod, l, k, r0v, pts, ptt
    real                    :: crhod, r0v, pts, ptt
    real                    :: calpha2,cbeta2,crho2,comega2    ! 2012_11_26 . scb for debugging
    integer                 :: l, k
! added end
    
    calpha2=calpha
    cbeta2=cbeta
    crho2=crho
    comega2=comega
! solves the linear system by preconditioned BiCGSTAB Algorithm
    crhod=crho
    crho=scalprod(ng2,vr0,vr)
    cbeta=crho*calpha/(crhod*comega)

    do k=1,nz
    do l=1,nxy
        vp(1,l,k)=vr(1,l,k)+cbeta*(vp(1,l,k)-comega*vv(1,l,k))
        vp(2,l,k)=vr(2,l,k)+cbeta*(vp(2,l,k)-comega*vv(2,l,k))
    enddo
    enddo



! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
!    call minv2g(vp,vy)
!    call axb2g(vy,vv)
    if(ifrect) then
        call minv2g(vp,vy)
        call axb2g(vy,vv)
    else
        call minvhex(vp,vy)
        call axbhex(vy,vv)
    endif
! added end

    r0v=scalprod(ng2,vr0,vv)
    !calpha=crho/r0v   ! 2014_05_20 . scb
    calpha=crho/(r0v+1.e-15) 

    do k=1,nz
    do l=1,nxy
        vs(1,l,k)=vr(1,l,k)-calpha*vv(1,l,k)
        vs(2,l,k)=vr(2,l,k)-calpha*vv(2,l,k)
    enddo
    enddo

! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
!    call minv2g(vs,vz)
!    call axb2g(vz,vt)
    if(ifrect) then
        call minv2g(vs,vz)
        call axb2g(vz,vt)
    else
        call minvhex(vs,vz)
        call axbhex(vz,vt)
    endif
! added end

    pts=scalprod(ng2,vs,vt)
    ptt=scalprod(ng2,vt,vt)

    comega=0
    if(ptt .ne. 0) comega=pts/ptt

    r2=0
    do k=1,nz
    do l=1,nxy
        phi(1,l,k)=phi(1,l,k)+calpha*vy(1,l,k)+comega*vz(1,l,k)
        phi(2,l,k)=phi(2,l,k)+calpha*vy(2,l,k)+comega*vz(2,l,k)
        vr(1,l,k)=vs(1,l,k)-comega*vt(1,l,k)
        vr(2,l,k)=vs(2,l,k)-comega*vt(2,l,k)
        r2=r2+vr(1,l,k)*vr(1,l,k)+vr(2,l,k)*vr(2,l,k)
    enddo
    enddo
    if (r20/=0.d0) then                        ! 2016.9.20. jjh
      r2=sqrt(r2)/r20
    else if (r20==0.d0 .and. r2==0.d0) then    ! 2016.9.20. jjh
      r2=0.d0                                  ! 2016.9.20. jjh
    end if                                     ! 2016.9.20. jjh

    return
end subroutine