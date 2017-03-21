! 2012_08_22 . scb . module for using BDF method
    module bdf

      type bdftype
        sequence
        real,pointer :: bdfcoef(:),phibdf(:,:,:,:),phifbdf(:,:,:,:)
        real,pointer :: phibdf2(:,:,:,:),phifbdf2(:,:,:,:)
      end type

      type(bdftype) :: bdfarray

      integer :: maxiter, imaxiter, ibound, bdforder_mod=5

      real,pointer,dimension(:) :: deltmarray
      real,pointer,dimension(:,:,:) :: phisave,phifsave,difphi
      real,pointer,dimension(:,:) :: psisave,difpsi
      real(8) :: invorder
      real(8) :: unittm

      ! power update? => is used for step size control method
      logical :: automain=.false.,flagmain=.true.
      logical :: flagbdf=.true.   ! default : true . 2013_01_28 . scb
      logical :: deltmupd=.FALSE.

! 2012_08_24  . scb      
      integer :: bdfordersave
      real :: deltmsave,error
      logical :: updneed=.true.
      
      real, pointer, dimension(:,:,:,:) :: phibdf,phifbdf
      real, pointer, dimension(:) :: bdfcoef
      real :: cont
      real, pointer, dimension(:,:) :: psi0
      
! 2012_12_06 . scb
      real, pointer, dimension(:,:,:,:) :: phibdf2,phifbdf2   ! 2nd momentum flux storage
      
! 2014_09_18 . scb
      real, pointer, dimension(:,:,:,:) :: srcdt, spsrcdt
      real, pointer, dimension(:,:,:,:,:) :: srcdtbdf
      logical  :: flagsrcdt=.false.
! added end      

! 2014_09_24 . scb
      real, pointer :: dtcff(:,:,:,:,:,:)
      real, pointer :: dtcff2(:,:,:,:,:,:,:)  ! 2014_10_06 . scb
      logical  :: flagsrcdt2=.false.    ! 2014_09_24 . scb
! added end     
      
    end module      