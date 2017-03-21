!subroutine updsrccff(idir,l,k,m,psisrc,reigv,isrcsign,srccff)
subroutine updsrccff(idir,l,k,m,psisrc,reigv,isrcsign,srccff,iftran)   ! 2014_09_26 . scb
    use const
    use senm2n
    use nodal,  only : trlcff0,trlcff1,trlcff2
    use xsec,   only : xssf,xssfs,xssfe,xschif
    use bdf,    only : flagsrcdt2    ! 2014_09_26 . scb
    
    integer                 :: idir,l,k,m,isrcsign
    real                    :: psisrc(0:4)
    real                    :: reigv
    real                    :: srccff(0:4)
    logical                 :: iftran     ! 2014_09_26 . scb
    
!   local var
    integer                 :: ms
    real                    :: sc1g1n(0:4),fissrc(0:4) ! scattering, fission
    real                    :: srcdtcff(0:2)    ! 2014_09_26 . scb

    ! construct scattering src
    sc1g1n=0.0      
    do ms=xssfs(m,l,k),xssfe(m,l,k)
    sc1g1n(0:4)=sc1g1n(0:4)+xssf(ms,m,l,k)*phicff(0:4,ms,l,k,idir)
    enddo

    ! total src
    srccff(0:4)=reigv*xschif(m,l,k)*psisrc(0:4)+sc1g1n(0:4)
    srccff(0)=srccff(0)-trlcff0(m,l,k,idir)
    srccff(1)=srccff(1)-trlcff1(m,l,k,idir)
    srccff(2)=srccff(2)-trlcff2(m,l,k,idir)

! 2014_09_26 . scb    
    if(iftran .and. flagsrcdt2) call updsrcdtcff(idir,l,k,m,srccff(0:2))
! added end      
    
    srccff(1)=isrcsign*srccff(1)
    srccff(3)=isrcsign*srccff(3)
    

    return 

end subroutine 