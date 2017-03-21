! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
module hexfdm

! Set FDM Linear System for Hexagonal Geometry

   use const
   use geomhex
   use xsec
   use allocs
   implicit none

   real, pointer, dimension(:,:,:)     :: dtilrfdm, psifdm
   real, pointer, dimension(:,:,:)     :: dtilzfdm   
   real, pointer, dimension(:,:,:,:)   :: diagfdm, phifdm, srcfdm
   real, pointer, dimension(:,:,:,:,:) :: cczfdm
   real, pointer, dimension(:,:,:,:,:) :: ccrfdm

   interface

   subroutine drive_hexfdm(epseig,epsl2)
      real :: epseig,epsl2
   end

   subroutine setls_hexfdm

   end

   subroutine upddtil_hexfdm

   end

   subroutine solvels_hexfdm(epseig,epsl2)
      real :: epseig,epsl2
   end

   subroutine sol_power(reigv,phi,errl2,erreig,domr,psi)
      real :: errl2,erreig,domr
      real,pointer :: phi(:,:,:,:)
      integer :: iout
      real :: reigv
      real,pointer :: psi(:,:,:)   
   end

   subroutine sol_gs(reigv,psi,phi)
      real,pointer :: psi(:,:,:)
      real :: reigv
      real,pointer :: phi(:,:,:,:)
   end
                                   
   end interface


   contains

   subroutine malloc_hexfdm
      ! dtil
		call dmalloc(dtilrfdm,ng,nxsfc,nz)
		call dmalloc(dtilzfdm,ng,nassy,nz+1)
      ! setls
		call dmalloc(diagfdm,ng,nsub,nassy,nz)
      call dmalloc(ccrfdm,3,ng,nsub,nassy,nz)
      call dmalloc(cczfdm,2,ng,nsub,nassy,nz)
      !
      call dmalloc0(psifdm,0,nsub,0,nassy,0,nz)
      call dmalloc0(phifdm,1,ng,0,nsub,0,nassy,0,nz)
      call dmalloc0(srcfdm,1,ng,0,nsub,0,nassy,0,nz)
      return
   end subroutine


end module
