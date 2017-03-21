    module allocs

    integer(8) :: nbytesf = 0

    interface dmalloc
    module procedure mallocf1
    module procedure mallocd1
    module procedure malloci1
    module procedure mallocl1
    module procedure mallocf2
    module procedure mallocd2
    module procedure malloci2
    module procedure mallocl2
    module procedure mallocf3
    module procedure mallocd3
    module procedure malloci3
    module procedure mallocf4
    module procedure mallocd4
    module procedure malloci4
    module procedure mallocf5
    module procedure mallocd5
    module procedure malloci5
    end interface
    interface dmalloc0
    module procedure mallocf01
    module procedure mallocd01
    module procedure malloci01
    module procedure mallocf02
    module procedure mallocd02
    module procedure malloci02
    module procedure mallocf03
    module procedure mallocd03
    module procedure malloci03
    module procedure mallocf04
    module procedure mallocd04
    module procedure malloci04
    module procedure mallocf05
    module procedure mallocd05
    module procedure malloci05
    module procedure mallocf06
    module procedure mallocd06
    module procedure malloci06
! 2014_10_06 . scb    
    module procedure mallocf07
    module procedure mallocd07
    module procedure malloci07
! added end    
    end interface
    !
    contains

    subroutine mallocf1(a,n1)

    real(4),pointer :: a(:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd1(a,n1)

    real(8),pointer :: a(:)

    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci1(a,n1)

    integer,pointer :: a(:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocl1(a,n1)

    logical,pointer :: a(:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1))
    a=.FALSE.
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocf2(a,n1,n2)

    real(4),pointer :: a(:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd2(a,n1,n2)

    real(8),pointer :: a(:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci2(a,n1,n2)

    integer,pointer :: a(:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocl2(a,n1,n2)

    logical,pointer :: a(:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2))
    a=.FALSE.
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocf3(a,n1,n2,n3)

    real(4),pointer :: a(:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd3(a,n1,n2,n3)

    real(8),pointer :: a(:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci3(a,n1,n2,n3)

    integer,pointer :: a(:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocf4(a,n1,n2,n3,n4)

    real(4),pointer :: a(:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd4(a,n1,n2,n3,n4)

    real(8),pointer :: a(:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci4(a,n1,n2,n3,n4)

    integer,pointer :: a(:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine

    subroutine mallocf5(a,n1,n2,n3,n4,n5)

    real(4),pointer :: a(:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4,n5))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd5(a,n1,n2,n3,n4,n5)

    real(8),pointer :: a(:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4,n5))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci5(a,n1,n2,n3,n4,n5)

    integer,pointer :: a(:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4,n5))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    
! 2014_10_06 . scb
    subroutine mallocf6(a,n1,n2,n3,n4,n5,n6)

    real(4),pointer :: a(:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4,n5,n6))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd6(a,n1,n2,n3,n4,n5,n6)

    real(8),pointer :: a(:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4,n5,n6))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci6(a,n1,n2,n3,n4,n5,n6)

    integer,pointer :: a(:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4,n5,n6))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    
    subroutine mallocf7(a,n1,n2,n3,n4,n5,n6,n7)

    real(4),pointer :: a(:,:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4,n5,n6,n7))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd7(a,n1,n2,n3,n4,n5,n6,n7)

    real(8),pointer :: a(:,:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4,n5,n6,n7))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci7(a,n1,n2,n3,n4,n5,n6,n7)

    integer,pointer :: a(:,:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(n1,n2,n3,n4,n5,n6,n7))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
! added end

    !
    subroutine mallocf01(a,nb1,ne1)

    real(4),pointer :: a(:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd01(a,nb1,ne1)

    real(8),pointer :: a(:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci01(a,nb1,ne1)

    integer,pointer :: a(:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocf02(a,nb1,ne1,nb2,ne2)

    real(4),pointer :: a(:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd02(a,nb1,ne1,nb2,ne2)

    real(8),pointer :: a(:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci02(a,nb1,ne1,nb2,ne2)

    integer,pointer :: a(:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocf03(a,nb1,ne1,nb2,ne2,nb3,ne3)

    real(4),pointer :: a(:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd03(a,nb1,ne1,nb2,ne2,nb3,ne3)

    real(8),pointer :: a(:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci03(a,nb1,ne1,nb2,ne2,nb3,ne3)

    integer,pointer :: a(:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocf04(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4)

    real(4),pointer :: a(:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd04(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4)

    real(8),pointer :: a(:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci04(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4)

    integer,pointer :: a(:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocf05(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5)

    real(4),pointer :: a(:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd05(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5)

    real(8),pointer :: a(:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5))
    a=0
    nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci05(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5)

    integer,pointer :: a(:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
    allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5))
    a=0
    nbytesf=nbytesf+size(a)
    end subroutine
  
    subroutine mallocf06(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5,nb6,ne6)
      real(4),pointer :: a(:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
      allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5,nb6:ne6))
      a=0
      nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd06(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5,nb6,ne6)
      real(8),pointer :: a(:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
      allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5,nb6:ne6))
      a=0
      nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci06(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5,nb6,ne6)
      integer,pointer :: a(:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
      allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5,nb6:ne6))
      a=0
      nbytesf=nbytesf+size(a)
    end subroutine

! 2014_10_06 . scb
    subroutine mallocf07(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5,nb6,ne6,nb7,ne7)
      real(4),pointer :: a(:,:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
      allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5,nb6:ne6,nb7:ne7))
      a=0
      nbytesf=nbytesf+size(a)
    end subroutine
    !
    subroutine mallocd07(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5,nb6,ne6,nb7,ne7)
      real(8),pointer :: a(:,:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
      allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5,nb6:ne6,nb7:ne7))
      a=0
      nbytesf=nbytesf+2*size(a)
    end subroutine
    !
    subroutine malloci07(a,nb1,ne1,nb2,ne2,nb3,ne3,nb4,ne4,nb5,ne5,nb6,ne6,nb7,ne7)
      integer,pointer :: a(:,:,:,:,:,:,:)
    
    if(associated(a)) return   ! 2014_12_18 . scb
    
      allocate(a(nb1:ne1,nb2:ne2,nb3:ne3,nb4:ne4,nb5:ne5,nb6:ne6,nb7:ne7))
      a=0
      nbytesf=nbytesf+size(a)
    end subroutine
! added end    
    end module