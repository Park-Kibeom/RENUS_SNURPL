! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
module hexfdm3

! Set FDM Linear System for Hexagonal Geometry

   use const
   use geomhex
   use xsec
   use allocs
   implicit none

   real, pointer, dimension(:,:,:)     :: dtilrfdm, dtilrfdm2, psifdm
   real, pointer, dimension(:,:,:)     :: dtilzfdm, dtilzfdm2
   real, pointer, dimension(:,:,:,:)   :: phifdm, phifdm2, srcfdm, srcfdm2
   real, pointer, dimension(:,:,:,:,:) :: ccrfdm, ccrfdm2, diagfdm
   real, pointer, dimension(:,:,:,:,:) :: cczfdm, cczfdm2

   interface

   subroutine drive_hexfdm3(epseig,epsl2)
      real :: epseig,epsl2
   end

   subroutine setls_hexfdm3

   end

   subroutine upddtil_hexfdm3

   end

   subroutine solvels_hexfdm3(epseig,epsl2)
      real :: epseig,epsl2
   end

   subroutine sol_power3(reigv,phi,phi2,errl2,erreig,psi)
      real :: errl2,erreig
      real,pointer :: phi(:,:,:,:), phi2(:,:,:,:)
      integer :: iout
      real :: reigv
      real,pointer :: psi(:,:,:)   
   end

   subroutine sol_gs3(reigv,psi,phi,phi2)
      real,pointer :: psi(:,:,:)
      real :: reigv
      real,pointer :: phi(:,:,:,:), phi2(:,:,:,:)
   end
                               
   end interface


   contains

   subroutine malloc_hexfdm3
      ! dtil
		call dmalloc(dtilrfdm,ng,nxsfc,nz)
		call dmalloc(dtilrfdm2,ng,nxsfc,nz)

		call dmalloc(dtilzfdm,ng,nassy,nz+1)
      call dmalloc(dtilzfdm2,ng,nassy,nz+1)
      ! setls
		call dmalloc(diagfdm,4,ng,nsub,nassy,nz)
      call dmalloc(ccrfdm,3,ng,nsub,nassy,nz)
      call dmalloc(ccrfdm2,3,ng,nsub,nassy,nz)

      call dmalloc(cczfdm,2,ng,nsub,nassy,nz)
      call dmalloc(cczfdm2,2,ng,nsub,nassy,nz)
      !
      call dmalloc0(psifdm,0,nsub,0,nassy,0,nz)
      call dmalloc0(phifdm,1,ng,0,nsub,0,nassy,0,nz)
      call dmalloc0(phifdm2,1,ng,0,nsub,0,nassy,0,nz)

      call dmalloc0(srcfdm,1,ng,0,nsub,0,nassy,0,nz)
      call dmalloc0(srcfdm2,1,ng,0,nsub,0,nassy,0,nz)
      return
   end subroutine


end module
