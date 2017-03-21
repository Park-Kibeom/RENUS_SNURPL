! added in ARTOS ver. 0.2 . 2012_07_24 by SCB     
   subroutine sol_gs(reigv,psi,phi)

! Solve Linear System for Hexagonal Geometry
	   use const
      use geom,   only : neibz
	   use geomhex
	   use hexfdm
      use sfam_cntl,  only : nupmax,maxupscatgr
      use xsec,       only : xschif,xssf,xssfs,xssfe,xsnff
	   implicit none

! in
      real,pointer :: psi(:,:,:)
      real :: reigv
! out
      real,pointer :: phi(:,:,:,:)
! local
      integer :: m, l, k, n, iin
      integer :: mbeg, ms, idir, idirz
      integer :: iout, iup, negative

      real    :: errphi, errphi2, phi2, phid, ss
      real    :: rerrl2, repsl2=1.e-5
      
! inner iteration
      do iin=1,20
         phi2=0.
         errphi2=0.
         do k=1,nz
            do l=1,nassy
               do n=1,nsub
                  do m=1,ng
                     srcfdm(m,n,l,k)=reigv*xschif(m,l,k)*psi(n,l,k)
                     ss=0
                     do ms=xssfs(m,l,k),xssfe(m,l,k)
                        ss=ss+phi(ms,n,l,k)*xssf(ms,m,l,k)
                     enddo
                     srcfdm(m,n,l,k)=srcfdm(m,n,l,k)+ss*volnodet(n,l,k)
                  enddo ! m

                  do m=1,ng
                     do idir=1,3
                        srcfdm(m,n,l,k)=srcfdm(m,n,l,k)  &
                        -ccrfdm(idir,m,n,l,k)*phi(m,neigtri(idir,n,l),neigtria(idir,n,l),k)  
                     enddo
                     do idirz=1,2
                        srcfdm(m,n,l,k)=srcfdm(m,n,l,k)  &
                        -cczfdm(idirz,m,n,l,k)*phi(m,n,l,neibz(idirz,k))
                     enddo

                     phid=phi(m,n,l,k)
                     phi(m,n,l,k)=srcfdm(m,n,l,k)/diagfdm(m,n,l,k)

                     errphi=phid-phi(m,n,l,k)
                     errphi2=errphi2+errphi*errphi
                     phi2=phi2+phi(m,n,l,k)*phi(m,n,l,k)
                  enddo
               enddo ! n
            enddo ! l
         enddo ! k 
         rerrl2=sqrt(errphi2/phi2)
         if(rerrl2.lt.repsl2 .and. iin.ge.10) then
            exit  
         endif     
      enddo ! iin

      return
   end subroutine
