! 2012_09_07 . scb
    subroutine initautodt(deltm0,bdforder,ng,nxy,nz)
      
      use allocs
      use bdf
      
      real :: deltm0
      integer :: bdforder, ng, nxy, nz
      
      integer,parameter :: ng2=2
      
      deltmsave=deltm0
      maxiter=bdforder
      imaxiter=0
      ibound=1
      
      call dmalloc0(phisave,1,ng2,0,nxy,0,nz)
      call dmalloc(psisave,nxy,nz)
      call dmalloc(psi0,nxy,nz)
      if(plevel.ge.0.9999999*1.e-2)  psi0=psif
      call dmalloc(difpsi,nxy,nz)
      invorder=1./bdforder
      if (ng.ne.ng2) then
        call dmalloc0(phifsave,1,ng,0,nxy,0,nz)
      ! difphi : difference between n-1th and nth bdf solutions
        call dmalloc0(difphi,1,ng,0,nxy,0,nz)
      else
        call dmalloc0(difphi,1,ng2,0,nxy,0,nz)
      endif             
      
    end subroutine