module decusping3n
    use const
    use allocs
    implicit none
    
    integer,parameter         :: nmesh = 10, ntfine = 3*10
    integer,parameter         :: IXSD=1,IXSTR=2,IXSA=3,IXST=4,IXSNF=5,IXSKF=6,IXSF=7,IXSS=8,NXS=9
    real                      :: xsunrodded(ng2,NXS), xsrodded(ng2,NXS)
    
    integer                   :: nfine(4)
    real                      :: hfine(4)
    real,pointer              :: trlfine(:,:)
    real,pointer              :: diag1d(:,:), bot1d(:,:), top1d(:,:), src1d(:,:)
    real,pointer              :: phi1d(:,:)
    real,pointer              :: fwght(:,:)
    
    contains
    subroutine mallocdecusping()
        use geom,     only : nxy
        call dmalloc(trlfine,ng2,ntfine)
        call dmalloc(diag1d,4,ntfine)
        call dmalloc(bot1d,ng2,ntfine)
        call dmalloc(top1d,ng2,ntfine)
        call dmalloc(src1d,ng2,ntfine)
        call dmalloc(phi1d,ng2,ntfine)
        call dmalloc(fwght,ng2,nxy)
        fwght = 1
    end subroutine
end module