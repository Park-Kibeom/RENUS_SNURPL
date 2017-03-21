module geom
! dimension parameters
    use const
    use sfam_cntl, only : ifrect      ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    implicit none

    integer               :: ng                   ! number of group
    integer               :: nx, ny, nz, nxy      ! number of nodes in x, y and z direction and total nodes
    integer               :: nxya, nxa, nya
    integer               :: nza                  ! 2013_07_15 . scb
    integer               :: nzp1                 ! 1 + number of assembly-wise or node-wise planes
    integer               :: nsurf                ! number of node-wise surfaces
    integer               :: ndir                 ! the current dimension.
    character*4           :: symopt
    integer               :: isymang, isymloc
    integer,parameter     :: SURFDIR(2)     = (/RIGHT,LEFT/)
    integer,parameter     :: SURFSGN(2)     = (/PLUS,PLUS/)
    integer               :: rotsurfdir(2)  = (/RIGHT,LEFT/)  ! surface direction (left or rght) with rotational geometry
    integer               :: rotsurfsgn(2)  = (/PLUS,PLUS/)
    integer               :: rotdir(2)      = (/XDIR, YDIR/)

    integer,pointer,dimension(:) :: nxs,nxe,nys,nye,nxsf,nxef

    integer :: kfbeg,kfend
    integer :: jfbeg,jfend     ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB  
    integer,pointer,dimension(:,:)  :: nodel, neibz, neibr, lsfc

    ! albedo
    real :: albedo(2,ndirmax)

    ! related to node size 
    real                            :: volfuel, volcore
    real,pointer,dimension(:,:)     :: volnode
    real,pointer,dimension(:,:,:)   :: hmesh
    real,pointer,dimension(:,:,:)   :: rhmesh   ! 2013_07_15 . scb
    
! added in ARTOS ver. 0.2 . 2012_07_20 by SCB  
    integer,pointer,dimension(:)    :: ltola,ktoka,iassytyp ! DECAY HEAT
    logical,pointer,dimension(:)    :: iffuela
! added end

    integer,pointer,dimension(:) :: ltoi, ltoj  ! 2013_07_16 . scb

    logical  :: ifrecsp3   ! 2013_07_05 . scb

    contains

    subroutine setgeom( ngl,nxl,nyl,nzl,nxyl,ndirl,symoptl,isymangl,  &
                        isymlocl,nxsl,nxel,nysl,nyel,nxfsl,nxfel,       &
                        kfbegl,kfendl,                                  &
                        jfbegl,jfendl,                                  &   ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB  
                        nodell,                                         &
                        hmeshl,volnodel,volcorel,volfuell,              &
                        albedol,                                        &
                        ltolal,ktokal,iassytypl,iffuelal,               &   ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB  ! DECAY HEAT
                        ifrecsp3l,                                      &   ! 2013_07_05 . scb
                        nxal,nyal,nzal,rhmeshl,ltoil,ltojl)     ! 2013_07_15 . scb
        use allocs
        
        integer                             :: ngl,nxl,nyl,nzl,nxyl,ndirl
        integer                             :: isymlocl,isymangl
        character*4                         :: symoptl
        integer ,pointer,dimension(:)       :: nxsl,nxel,nysl,nyel,nxfsl,nxfel
        integer                             :: kfbegl,kfendl
        integer                             :: jfbegl,jfendl     ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB  
        integer ,pointer,dimension(:,:)     :: nodell
        real    ,pointer,dimension(:,:,:)   :: hmeshl
        real    ,pointer,dimension(:,:,:)   :: rhmeshl    ! 2013_07_15 . scb
        real    ,pointer,dimension(:,:)     :: volnodel
        real                                :: volcorel,volfuell
        real                                :: albedol(2,ndirmax)
! added in ARTOS ver. 0.2 . 2012_07_20 by SCB  
        integer ,pointer,dimension(:)       :: ltolal,ktokal,iassytypl ! DECAY HEAT
        logical,pointer,dimension(:)        :: iffuelal
! added end
        logical                             :: ifrecsp3l   ! 2013_07_05 . scb
        integer,pointer,dimension(:)        :: ltoil, ltojl  ! 2013_07_16 . scb        
        
        integer                             :: l,i,j,lrot,k,kb,is
        
        integer :: nxal,nyal,nzal    ! 2013_07_15 . scb
        
        nxa=nxal    ! 2013_07_15 . scb
        nya=nyal    ! 2013_07_15 . scb
        nza=nzal    ! 2013_07_15 . scb

        ng=ngl;nx=nxl;ny=nyl;nz=nzl;nzp1=nz+1;nxy=nxyl;ndir=ndirl
        symopt=symoptl
        isymang=isymangl
        isymloc=isymlocl
        albedo=albedol

        nxs     =>  nxsl
        nxe     =>  nxel
        nys     =>  nysl
        nye     =>  nyel
        nxsf    =>  nxfsl
        nxef    =>  nxfel
                
        kfbeg   =   kfbegl
        kfend   =   kfendl
! added in ARTOS ver. 0.2 . 2012_07_20 by SCB  
        jfbeg   =   jfbegl
        jfend   =   jfendl
! added end
        nodel   =>  nodell

        volnode =>  volnodel
        volcore =   volcorel
        volfuel =   volfuell
        hmesh   =>  hmeshl 
        rhmesh  =>  rhmeshl   ! 2013_07_15 . scb
		  
! added in ARTOS ver. 0.2 . 2012_07_20 by SCB 
        ltola   =>  ltolal    ! DECAY HEAT
        ktoka   =>  ktokal    ! DECAY HEAT
        iassytyp=>  iassytypl ! DECAY HEAT
        iffuela =>  iffuelal
! added end
        ifrecsp3 = ifrecsp3l  ! 2013_07_05 . scb
        
        ltoi => ltoil  ! 2013_07_16 . scb
        ltoj => ltojl  ! 2013_07_16 . scb
		
		    if(.not.ifrect) goto 100   ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB  
! assigning neighbor node
        call dmalloc(neibr,nrdir2,nxy)
!
        do j=1,ny
        do i=nxs(j),nxe(j)
            l=nodel(i,j)
            if(i.ne.nxs(j)) neibr(1,l)=nodel(i-1,j)
            if(i.ne.nxe(j)) neibr(2,l)=nodel(i+1,j)
            if(j.ne.1)      neibr(3,l)=nodel(i,j-1)
            if(j.ne.ny)     neibr(4,l)=nodel(i,j+1)
        enddo
        enddo

        rotsurfdir    = SURFDIR
        rotsurfsgn(:) = SURFSGN
        if(symopt .eq. 'ROT') then
            call setrotgeom
        endif

    100 continue       ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB              
                
        call dmalloc(neibz,2,nz)        
        k=1
        neibz(1,k)=0
        kb=k
        do k=2,nz
            neibz(2,kb)=k
            neibz(1,k)=kb
            kb=k
        enddo 
        k=nz
        neibz(2,k)=0 

        if(.not.ifrect) return       ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB

! assigning surface index.
        call dmalloc(lsfc,nrdir2,nxy)        
        is=0
        do j=1,ny
            is=is+1
            do i=nxs(j),nxe(j)
                l=nodel(i,j)
                lsfc(1,l)=is
                is=is+1
                lsfc(2,l)=is
            enddo
        enddo
        do i=1,nx
            is=is+1
            do j=nys(i),nye(i)
                l=nodel(i,j)
                lsfc(3,l)=is
                is=is+1
                lsfc(4,l)=is
            enddo
        enddo   

! count the number of surfaces         
        nsurf=0
        do j=1,ny
            nsurf=nsurf+nxe(j)-nxs(j)+2
        enddo
        do i=1,nx
            nsurf=nsurf+nye(i)-nys(i)+2
        enddo      
    end subroutine
    
    subroutine setrotgeom
        integer                             :: l,i,j,lrot
        
        if(isymang.eq.90) then
            rotdir(XDIR)=YDIR
            rotdir(YDIR)=XDIR
            select case (isymloc) 
            case(1)
                i=1
                do j=nys(i), nye(i)
                    l=nodel(i,j)
                    lrot=nodel(nx-j+1,ny)
                    neibr(1,l)=lrot
                    neibr(4,lrot)=l
                enddo
            case(2)
                j=ny
                do i=nxs(j), nxe(j)
                    l=nodel(i,j)
                    lrot=nodel(ny,i)
                    neibr(4,l)=lrot
                    neibr(2,lrot)=l
                enddo
                rotsurfdir(LEFT)    = RIGHT
                rotsurfdir(RIGHT)   = RIGHT                
                rotsurfsgn(RIGHT)   = MINUS
            case(3)
                i=nx
                do j=nys(i), nye(i)
                    l=nodel(i,j)
                    lrot=nodel(ny-j+1,1)
                    neibr(2,l)=lrot
                    neibr(3,lrot)=l
                enddo
            case(4)
                j=1
                do i=nxs(j), nxe(j)
                    l=nodel(i,j)
                    lrot=nodel(1,i)
                    neibr(3,l)=lrot
                    neibr(1,lrot)=l
                enddo
                rotsurfdir(LEFT)   = LEFT
                rotsurfdir(RIGHT)  = LEFT
                rotsurfsgn(LEFT)   = MINUS
            end select   
        elseif(isymang.eq.180) then
            select case (isymloc) 
            case(1)
                j=1
                do i=nxs(j), nxe(j)
                    l=nodel(i,j)
                    lrot=nodel(nx-i+1,j)
                    neibr(NORTH,l)=lrot
                    neibr(NORTH,lrot)=l
                enddo
                rotsurfdir(LEFT)    = LEFT
                rotsurfdir(RIGHT)   = LEFT                
                rotsurfsgn(LEFT)   = MINUS
            case(2)
                i=1
                do j=nys(i), nye(i)
                    l=nodel(i,j)
                    lrot=nodel(i,ny-j+1)
                    neibr(WEST,l)=lrot
                    neibr(WEST,lrot)=l
                enddo
                rotsurfdir(LEFT)    = LEFT
                rotsurfdir(RIGHT)   = LEFT                
                rotsurfsgn(LEFT)   = MINUS
            case(3)
                j=ny
                do i=nxs(j), nxe(j)
                    l=nodel(i,j)
                    lrot=nodel(nx-i+1,j)
                    neibr(SOUTH,l)=lrot
                    neibr(SOUTH,lrot)=l
                enddo
                rotsurfdir(LEFT)    = RIGHT
                rotsurfdir(RIGHT)   = RIGHT                
                rotsurfsgn(RIGHT)   = MINUS
            case(4)
                i=nx
                do j=nys(i), nye(i)
                    l=nodel(i,j)
                    lrot=nodel(i,ny-j+1)
                    neibr(EAST,l)=lrot
                    neibr(EAST,lrot)=l
                enddo
                rotsurfdir(LEFT)    = RIGHT
                rotsurfdir(RIGHT)   = RIGHT                
                rotsurfsgn(RIGHT)   = MINUS
            end select           
        endif 
    end subroutine
end module