! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
module fdm_solver

   use const
   use geomhex
   use xsec
   use allocs
   implicit none

   real :: domr
! cheby
   real, pointer, dimension(:,:,:,:)   :: phifdmd, phifdmdd
   real, pointer, dimension(:,:,:)     :: psifdmd, psifdmdd
   logical :: chebyon=FALSE
   logical :: newcheby=FALSE
   logical :: firstdom
   real :: sigma,sigbar,erfratio,erfcb,erfcb0
   real :: alphacb,betacb
   integer :: mcp,mincheby

   real, parameter :: domrlo=0.5, domrhi=0.99999
   integer, parameter :: ncheby=10, icheby0=3

! lsor


   interface

   subroutine cheby(iffront,iout)
      logical :: iffront
      integer :: iout
   end subroutine


   subroutine sol_cheby(phi,psi,iout)
      real, pointer :: phi(:,:,:,:)
      real, pointer :: psi(:,:,:)
      integer :: iout
   end subroutine
                                    
   end interface


   contains

   subroutine malloc_fdmsolver
      ! cheby
      call dmalloc0(phifdmd,1,ng,0,nsub,0,nassy,0,nz)
      call dmalloc0(phifdmdd,1,ng,0,nsub,0,nassy,0,nz)
      call dmalloc0(psifdmd,0,nsub,0,nassy,0,nz)
      call dmalloc0(psifdmdd,0,nsub,0,nassy,0,nz)
      return
   end subroutine

end module
