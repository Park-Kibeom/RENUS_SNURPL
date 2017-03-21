subroutine updphicff(idir,l,k,m,psol,hsola,hsolb)
    use const
    use senm2n
    
    integer                 :: idir,l,k,m
    real, intent(in)        :: psol(0:4),hsola,hsolb


    phicff(0,m,l,k,idir)=phicnst(0,m,l,k,idir)*hsolb+psol(0)
    phicff(2,m,l,k,idir)=phicnst(2,m,l,k,idir)*hsolb+psol(2)
    phicff(4,m,l,k,idir)=phicnst(4,m,l,k,idir)*hsolb+psol(4)
    phicff(1,m,l,k,idir)=phicnst(1,m,l,k,idir)*hsola+psol(1)
    phicff(3,m,l,k,idir)=phicnst(3,m,l,k,idir)*hsola+psol(3)

    return
end subroutine
