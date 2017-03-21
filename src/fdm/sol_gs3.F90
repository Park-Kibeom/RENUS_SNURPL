! added in ARTOS ver. 0.2 . 2012_07_24 by SCB     
   subroutine sol_gs3(reigv,psi,phi,phi2)

! Solve Linear System for Hexagonal Geometry
	   use const
      use geom,   only : neibz
	   use geomhex
	   use hexfdm3
      use sfam_cntl,  only : nupmax,maxupscatgr
      use xsec,       only : xschif,xssf,xssfs,xssfe,xsnff
	   implicit none

! in
      real,pointer :: psi(:,:,:)
      real :: reigv
! out
      real,pointer :: phi(:,:,:,:), phi2(:,:,:,:)
! local
      integer :: m, l, k, n, iin
      integer :: ms, idir, idirz
      integer :: iout, negative

      real    :: phid, ss
      real,save :: r2o3=2./3.
      real    :: rdet, dm(4), src(ng), src2(ng)

      real    :: errl2, phil2, err
      real    :: rerrl2, repsl2=1.e-5
      
! inner iteration
      do iin=1,20
         errl2=0.
         phil2=0.
         do k=1,nz
            do l=1,nassy
               do n=1,nsub
                  do m=1,ng
                     src(m)=reigv*xschif(m,l,k)*psi(n,l,k)
                     ss=0
                     do ms=xssfs(m,l,k),xssfe(m,l,k)
                        ss=ss+(phi(ms,n,l,k)-2.*phi2(ms,n,l,k))*xssf(ms,m,l,k)
                     enddo
                     src(m)=src(m)+ss*volnodet(n,l,k)
                  enddo ! m

                  ! 2-nd source term                  
                  src2(:)=-r2o3*src(:)
                  do m=1,ng  
                     do idir=1,3
                        ! 0-th
                        src(m)=src(m)  &
                        -ccrfdm(idir,m,n,l,k)*phi(m,neigtri(idir,n,l),neigtria(idir,n,l),k) 
                        ! 2-nd 
                        src2(m)=src2(m)  &
                        -ccrfdm2(idir,m,n,l,k)*phi2(m,neigtri(idir,n,l),neigtria(idir,n,l),k)  
                     enddo

                     do idirz=1,2
                        ! 0-th
                        src(m)=src(m)  &
                        -cczfdm(idirz,m,n,l,k)*phi(m,n,l,neibz(idirz,k))
                        ! 2-nd
                        src2(m)=src2(m)  &
                        -cczfdm2(idirz,m,n,l,k)*phi2(m,n,l,neibz(idirz,k))
                     enddo
                     phid=phi(m,n,l,k)

                     dm(:)=diagfdm(:,m,n,l,k)          
                     rdet=1./(dm(1)*dm(4)-dm(2)*dm(3))

                     phi(m,n,l,k)=rdet*(dm(4)*src(m)-dm(2)*src2(m))
                     phi2(m,n,l,k)=rdet*(-dm(3)*src(m)+dm(1)*src2(m))

                     err=phid-phi(m,n,l,k)
                     errl2=errl2+err*err
                     phil2=phil2+phi(m,n,l,k)*phi(m,n,l,k)
                  enddo
               enddo ! n
            enddo ! l
         enddo ! k       
         rerrl2=sqrt(errl2/phil2)
         if(rerrl2.lt.repsl2 .and. iin.ge.10) then
            exit  
         endif  
      enddo ! iin

      return
   end subroutine
