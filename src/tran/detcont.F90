! 2012_09_07 . scb
    subroutine detcont(ng,bdforder,deltm,errtol,safefac)
    
      use sfam,     only : phi, phif, psi
      use tran,     only : psit
      use bdf
    
      integer :: ng, bdforder
      real*8 :: deltm,errtol,safefac
      logical :: flag2g
      
      integer,parameter :: ng2=2
      
      flag2g=ng2.eq.ng
      
      flagmain=.false.
                
      phisave=phi
      psisave=psi
      
      psi=psit
      
      if(ng.ne.ng2)  phifsave=phif
      
      bdforder=bdforder-1

      call calbdfcoef(bdforder)
      
      call runtrtm(flag2g,deltm)
      
      call calcont(errtol,safefac)
      
      phi=phisave
      psi=psisave
      
      if(ng.ne.ng2) phif=phifsave
      
      bdforder=bdforder+1
      
      flagmain=.true.
      
!      if(cont>1.0.and.imaxiter.ge.maxiter.and.deltm.ne.maxdeltm) then
        
      
    end subroutine