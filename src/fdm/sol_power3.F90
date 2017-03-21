! added in ARTOS ver. 0.2 . 2012_07_24 by SCB     
   subroutine sol_power3(reigv,phi,phi2,errl2,erreig,psi)
                             
! Solve Linear System for Hexagonal Geometry
	   use const
	   use geomhex
      use xsec,       only : xsnff
      use fdm_solver, only : domr
	   implicit none
! in
      real :: errl2,erreig
      real,pointer :: phi(:,:,:,:), phi2(:,:,:,:)
      integer :: iout
! out
      real :: reigv
      real,pointer :: psi(:,:,:)   

      integer :: m, l, k, n
      real    :: err2d, err2, psipsid, psipsi, psid, err
      real    :: eigvd, phid, eigv
      real    :: relresid, resid, resid0

      eigv=1./reigv
      err2d=err2;err2=0;
      psipsid=0;psipsi=0

      do k=1,nz
      do l=1,nassy
      do n=1,nsub
         psid=psi(n,l,k)
         psi(n,l,k)=0
         do m=1,ng
            psi(n,l,k)=psi(n,l,k)+xsnff(m,l,k)*(phi(m,n,l,k)-2.*phi2(m,n,l,k))
         enddo
         psi(n,l,k)=psi(n,l,k)*volnodet(n,l,k)
         psipsi=psipsi+psi(n,l,k)*psi(n,l,k)               
         psipsid=psipsid+psi(n,l,k)*psid

         err=psid-psi(n,l,k)
         err2=err2+err*err
      enddo ! n
      enddo ! l
      enddo ! k

      ! estimate error reduction factor (dominance ratio if no extrapolation)
      if(err2d .ne. 0) domr=sqrt(err2/err2d)

      eigvd=eigv
      ! update eigenvalue
      eigv=eigv*psipsi/psipsid
      reigv=1/eigv        
    
      errl2=sqrt(err2/psipsid)
      erreig=abs(eigv-eigvd)/eigv

      return
   end subroutine
