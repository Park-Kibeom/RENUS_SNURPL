    module xedriver

      use precision, only : dp

      real(dp) :: lxe, lio, lpm
      real(dp) :: ypm
      real(dp) :: rtol=1.D-3, atol=1.D-12
      real(dp), dimension(:), pointer :: arate, frate      
      real(dp), dimension(:), pointer :: yxe, yio

      contains

      subroutine driver_setup(yio_,yxe_,arate_,frate_,lxe_,lio_,dim)
        integer :: dim
        real(dp) :: lxe_, lio_
        real(dp), target :: yxe_(dim), yio_(dim)
        real(dp), target :: arate_(dim), frate_(dim)

        lxe = lxe_
        lio = lio_
        yxe => yxe_
        yio => yio_
        arate => arate_
        frate => frate_

      end subroutine

!      function absorbtion(t)
        ! Calculate neutron absorbtion by xenon
!        real(dp) :: t
!        real(dp), dimension(dim) :: absorbtion

!        return
!      end function

!      function fission_rate(t)
        ! Calculate rate of fission
!        real(dp) :: t
!        real(dp), dimension(dim) :: fission_rate

!        return
!      end function


      function FUNC(x, t, dim)

        implicit none

        integer :: dim, n1
        real(dp) :: t
        real(dp), dimension(dim) :: x, FUNC

        n1 = dim/2

        FUNC(1:n1)     = -lxe*x(1:n1) - arate(:)*x(1:n1) + &
                          yxe(:) * frate(:) + lio*x(n1+1:dim)
        FUNC(n1+1:dim) = -lio*x(n1+1:dim) + yio(:) * frate(:)

      return
      end function FUNC


      function SOLVER(X0, t, alpha, h, F, dim)
      !
      ! [1 - alpha*h*A(t)] * X = F  =>  X = ...
        implicit none

        integer  :: dim, n1, i
        real(dp) :: alpha, h, t, ah
        real(dp) :: det
        real(dp) :: rhs_xe, rhs_io
        real(dp), dimension(dim) :: x0, f, SOLVER

        ah = alpha * h
        n1 = dim/2


        do i=1,dim/2
         
          rhs_xe = (f(i)+ah*yxe(i)*frate(i))
          rhs_io = (f(i+n1)+ah*yio(i)*frate(i))

          det = 1./(1.+ah*lxe+ah*arate(i))/(1.+ah*lio)
          SOLVER(i) = det * ((1.+ah*lio) * rhs_xe + ah*lio * rhs_io)
          SOLVER(i+n1) = det * (1.+ah*(lxe + arate(i))) * rhs_io

        enddo


!        write(*,*) 't, h = ', t, h


        return
      end function SOLVER


      FUNCTION LOCERR(dX, X0, dim)
      !
          IMPLICIT NONE

          EXTERNAL NORM

          INTEGER :: dim, i
          REAL(dp) :: LOCERR, NORM
          REAL(dp), DIMENSION(dim) :: X0, dX, ERR

          DO i = 1, dim
              ERR(i) = dX(i) / (RTOL * ABS(X0(i)) + ATOL)
          END DO

          LOCERR = NORM(ERR, dim, 2) / dim ** 0.5

      RETURN
      END FUNCTION LOCERR



      SUBROUTINE CALLBACK
      !
          IMPLICIT NONE

      END SUBROUTINE CALLBACK


      FUNCTION NORM(X, dim, p)
      !
      ! Âû÷èñëåíèå íîðìû âåêòîðà 
      ! X = ||X||p = Sum(|X(i)|**p, i) ** 1./p
      !
          USE PRECISION, ONLY : dp
          IMPLICIT NONE

          INTEGER :: i, dim, p
          REAL(dp) :: NORM, X(dim)

          NORM = 0.
          DO i = 1, dim
              NORM = NORM + ABS(X(i))**p
          END DO

          NORM = NORM ** (1./p)

      END FUNCTION NORM



    end module xedriver
