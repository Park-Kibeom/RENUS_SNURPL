! added in ARTOS ver. 0.2 . 2012_07_24 by SCB         
      subroutine minvhex(b,x)
 
! solve Mx=b for x given b
 
	   use const
	   use geomhex
	   use cmfdhex2g
	   use bicg2g, only : sol2dhex
	   implicit none

	   real,pointer :: x(:,:,:), b(:,:,:)
	   integer :: l, k, kp1

! forward solve

! for bottom plane
      do l=1,nxy
         dumrv(1,l)=b(1,l,1)
         dumrv(2,l)=b(2,l,1)
      enddo
      call sol2dhex(delinv(:,:,1),xlufac(:,:,:,1),dumrv,x(:,1:nxy,1))

      do k=2,nz
         do l=1,nxy
            dumrv(1,l)=b(1,l,k)+cmat(1,7,l,k)*x(1,l,k-1)
            dumrv(2,l)=b(2,l,k)+cmat(2,7,l,k)*x(2,l,k-1)
         enddo

         call sol2dhex(delinv(:,:,k),xlufac(:,:,:,k),dumrv,x(:,1:nxy,k))

      enddo

! backward solve
      do k=nz-1,1,-1
         kp1=k+1
         do l=1,nxy
            dumrv(1,l)=-cmat(1,8,l,k)*x(1,l,kp1)
            dumrv(2,l)=-cmat(2,8,l,k)*x(2,l,kp1)
         enddo

         call sol2dhex(delinv(:,:,k),xlufac(:,:,:,k),dumrv,dumrs)

         do l=1,nxy
            x(1,l,k)=x(1,l,k)-dumrs(1,l)
            x(2,l,k)=x(2,l,k)-dumrs(2,l)
         enddo
      enddo

      return
      end 

      subroutine sol2dhex(df,xlu,b,x)

! solve a 2d problem using precalculated LU factors
! df : diagonal term D of preconditioner
! xlu : off-diagonal term 

	   use const
	   use geomhex
	   use cmfdhex2g
	   implicit none

      real,pointer :: df(:,:),x(:,:),b(:,:)
	   real,pointer :: xlu(:,:,:)

	   integer :: j, l, ifr
      real :: sb(ng2)

!  forward solve
      do l=1,nxy
        sb(1)=0
        sb(2)=0
        do j=1,ilubnd(l)
          ifr=ipntr(j,l)
          sb(1)=sb(1)+xlu(1,j,l)*b(1,ifr)
          sb(2)=sb(2)+xlu(2,j,l)*b(2,ifr)
        enddo
        b(1,l)=b(1,l)-sb(1)
        b(2,l)=b(2,l)-sb(2)
      enddo
!  backward solve
      do l=nxy,1,-1
        sb(1)=b(1,l)
        sb(2)=b(2,l)
        do j=ilubnd(l)+1,ipntr(0,l)
          ifr=ipntr(j,l)
          sb(1)=sb(1)-xlu(1,j,l)*x(1,ifr)
          sb(2)=sb(2)-xlu(2,j,l)*x(2,ifr)
        enddo
        x(1,l)=df(1,l)*sb(1)+df(2,l)*sb(2)
        x(2,l)=df(3,l)*sb(1)+df(4,l)*sb(2)
      enddo

      return
      end
