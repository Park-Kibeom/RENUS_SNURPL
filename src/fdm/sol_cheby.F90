! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      subroutine sol_cheby(phi,psi,iout)

      use const
      use geomhex
      use hexfdm
      use fdm_solver
      implicit none
!
! extrapolate fission source using the chebyshev acceleration parameters
!
      integer :: iout
      real, pointer    :: phi(:,:,:,:), psi(:,:,:)
      integer :: k, l, m, n
      real    :: psil1, err, err2

! determine extrapolation paramters
      call cheby(TRUE,iout)
!
! extrapolate spectrum eigenvector
      if(chebyon) then
        psil1=0
        do k=1,nz
          do l=1,nxy
            do n=1,nsub
               do m=1,ng
                  phifdmdd(m,n,l,k)=phifdmd(m,n,l,k)
                  phifdmd(m,n,l,k)=phi(m,n,l,k)
                   
                  phi(m,n,l,k)=phifdmd(m,n,l,k) &
                  +alphacb*(phi(m,n,l,k)-phifdmd(m,n,l,k)) &
                  +betacb*(phifdmd(m,n,l,k)-phifdmdd(m,n,l,k))
               enddo
               psifdmdd(n,l,k)=psifdmd(n,l,k)
               psifdmd(n,l,k)=psi(n,l,k)

               psi(n,l,k)=psifdmd(n,l,k) &
                  +alphacb*(psi(n,l,k)-psifdmd(n,l,k)) &
                  +betacb*(psifdmd(n,l,k)-psifdmdd(n,l,k))

               psil1=psil1+psi(n,l,k)
               err=psifdmd(n,l,k)-psi(n,l,k)
               err2=err2+err*err
            enddo
          enddo
        enddo
      endif
! update and monitor dominance ratio
      call cheby(FALSE,iout)

      return
      end subroutine

