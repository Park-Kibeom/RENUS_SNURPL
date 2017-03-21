subroutine calby1n(idir,l,k,m,srccff,phifs,psol,hsolb)
    use const
    use senm2n
    use xsec,   only : xstf
    implicit none
    
    integer,intent(in)      :: idir,l,k,m
    real,   intent(in)      :: srccff(0:4)
    real,   intent(in)      :: phifs
    real,   intent(out)     :: psol(0:4),hsolb
    
    real                    :: kp2s,rkp2s,rxstfs,kps,cschkps ! local var
    
    kp2s=kp2(m,l,k,idir)
    rkp2s=rkp2(m,l,k,idir)
    rxstfs=1/xstf(m,l,k)
    cschkps=cschkp(m,l,k,idir)      
    kps=kp(m,l,k,idir)

!   obtain particular solution
    if(nsrcexp .eq. 2) then 
        psol(0)=(kp2s*srccff(0)+3*srccff(2))*rkp2s*rxstfs
        psol(1)=srccff(1)*rxstfs
        psol(2)=srccff(2)*rxstfs
        psol(3)=0
        psol(4)=0
    else
        psol(0)=(   kp2s*kp2s*srccff(0)+105*srccff(4)+kp2s*(        &
                        3*srccff(2)+10*srccff(4)                    &
                    )                                               &
                )*rkp2s*rkp2s*rxstfs

        psol(1)=(kp2s*srccff(1)+15*srccff(3))*rkp2s*rxstfs
        psol(2)=(kp2s*srccff(2)+35*srccff(4))*rkp2s*rxstfs
        psol(3)=srccff(3)*rxstfs
        psol(4)=srccff(4)*rxstfs
    endif

!   obtain homogeneous solution B.
!   Since sinh is odd function, B can be obtained without RIGHT node.
    hsolb=cschkps*kps*(phifs-psol(0))

end subroutine
