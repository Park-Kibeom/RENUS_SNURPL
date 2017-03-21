! 2012_09_07 . scb
    subroutine calcont(errtol,safefac)
    
      use sfam,     only : phi, psi
      use geom,     only : nz, nxy
      use bdf
      
      real(8) :: psisum,errtol,safefac
      
      error=0
      psisum=0
      
      difpsi=psi-psisave
      
      do k=1,nz
        do j=1,nxy
            error=error+(difpsi(j,k)*difpsi(j,k))
            psisum=psisum+(psi(j,k)-psi0(j,k))*(psi(j,k)-psi0(j,k))
        enddo
      enddo

      error=(error/psisum)**0.5
      
      cont=(errtol*safefac/error)**invorder
      
    end subroutine
