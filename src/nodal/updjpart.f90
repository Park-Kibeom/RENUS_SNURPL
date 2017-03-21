subroutine updjpart(jnet,phisfc)
    use const
    use nodal
    use sfam,  only : curil,curir,curol,curor
    use geom,   only : ndir,ng,nxy,nz

    real,pointer                :: jnet(:,:,:,:,:)
    real,pointer                :: phisfc(:,:,:,:,:)    

    integer                 :: m,l,k,idir
    
    do m=1, ng
    do k=1,nz
    do l=1,nxy
    do idir=1,ndir
        curil(idir,l,k,m)=RFOUR*phisfc(LEFT,m,l,k,idir)+HALF*jnet(LEFT,m,l,k,idir)
        curol(idir,l,k,m)=RFOUR*phisfc(LEFT,m,l,k,idir)-HALF*jnet(LEFT,m,l,k,idir)

        curir(idir,l,k,m)=RFOUR*phisfc(RIGHT,m,l,k,idir)-HALF*jnet(RIGHT,m,l,k,idir)
        curor(idir,l,k,m)=RFOUR*phisfc(RIGHT,m,l,k,idir)+HALF*jnet(RIGHT,m,l,k,idir)
    enddo
    enddo
    enddo
    enddo

return
end subroutine
