    module xedriver

      use precision, only : dp

      real(dp) :: lxe, lio, lpm
      real(dp) :: yxe, yio, ypm
      real(dp) :: absxe
      real(dp) :: rtol, atol

      contains

      subroutine set_pars(yio_,yxe_,absxe_,lxe_,lio_)
        real(dp) :: lxe_, lio_
        real(dp) :: yxe_, yio_
        real(dp) :: absxe_

        yio = yio_
        yxe = yxe_
        absxe = absxe_
        lxe = lxe_
        lio = lio_
        rtol = 1.D-3
        atol = 1.D-12
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

        integer :: dim
        real(dp) :: t
        real(dp), dimension(dim) :: x, FUNC

        FUNC(1) = -lxe*x(1) - absxe*x(1) + yxe + lio*x(2)
        FUNC(2) = -lio*x(2)              + yio

      return
      end function FUNC


      function SOLVER(X0, t, alpha, h, F, dim)
      !
      ! [1 - alpha*h*A(t)] * X = F  =>  X = ...
        implicit none

        integer  :: dim
        real(dp) :: alpha, h, t
        real(dp) :: det
        real(dp), dimension(dim) :: x0, f, SOLVER


        det = 1./(1.+alpha*h*lxe+alpha*h*absxe)/(1.+alpha*h*lio)
        SOLVER(1) = (1.+alpha*h*lio) * (f(1)+alpha*h*yxe) + &
                                  alpha*h*lio * (f(2)+alpha*h*yio)
        SOLVER(2) = (1.+alpha*h*(lxe + absxe)) * (f(2)+alpha*h*yio)

        SOLVER(:) = det * SOLVER(:)

!        write(*,*) 't, h = ', t, h
!        write(*,*) 'X0 = ', X0
!        write(*,*) 'RES= ', SOLVER
!        pause

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
