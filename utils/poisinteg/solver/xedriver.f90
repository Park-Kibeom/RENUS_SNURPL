    module xedriver

      use precision, only : dp
      
      integer :: dim
      real(dp) :: lxe, lio, lpm
      real(dp) :: ypm
      real(dp) :: rtol=1.D-9, atol=1.D-12
      real(dp), dimension(:), pointer :: varate, vfrate
      real(dp), dimension(:), pointer :: dvarate, dvfrate
      real(dp), dimension(:), pointer :: io0
      real(dp), dimension(:), pointer :: yxe, yio

      contains

      subroutine driver_setup(lxe_,lio_,dim_,yxe_,&
        yio_,varate_,vfrate_,dvarate_,dvfrate_,io0_)

        implicit none

        integer :: dim_
        real(dp) :: lxe_, lio_
        real(dp), target :: yxe_(dim_), yio_(dim_)
        real(dp), target :: varate_(dim_), vfrate_(dim_)
        real(dp), target :: dvarate_(dim_), dvfrate_(dim_)
        real(dp), target :: io0_(dim_)

        dim = dim_
        lxe = lxe_
        lio = lio_
        yxe => yxe_
        yio => yio_
        io0 => io0_
        varate => varate_
        vfrate => vfrate_
        dvarate => dvarate_
        dvfrate => dvfrate_

      end subroutine

      function arate(t)
        ! Calculate neutron absorbtion by xenon
        real(dp) :: t
        real(dp), dimension(dim) :: arate

        arate(:) = varate(:) + t * dvarate(:)

        return
      end function

      function frate(t)
        ! Calculate rate of fission
        real(dp) :: t
        real(dp), dimension(dim) :: frate

        frate(:) = vfrate(:) + t * dvfrate(:)

        return
      end function


      function rho_io(t)
        ! Calculate rate of fission
        real(dp) :: t
        real(dp) :: expt
        real(dp), dimension(dim) :: rho_io

        expt = exp(-lio*t)
        rho_io(:) = io0(:)*expt + vfrate(:)*yio(:)/lio * (1-expt)+ &
        dvfrate * (expt + (lio*t-1))/lio**2

!        write(*,*) t, expt
!        write(*,*) io0(1:5)
!        write(*,*) rho_io(1:5)
!        pause

        return
      end function


      function FUNC(x, t, dim)

        implicit none

        integer :: dim, n1
        real(dp) :: t
        real(dp), dimension(dim) :: x, FUNC

        FUNC(:) =-(lxe+arate(t))*x(:) + yxe(:)*frate(t) &
                                             + lio*rho_io(t)

      return
      end function FUNC


      function SOLVER(X0, t, alpha, h, F, dim)
      !
      ! [1 - alpha*h*A(t)] * X = F  =>  X = ...
        implicit none

        integer  :: dim, n1, i
        real(dp) :: alpha, h, t, ah
        real(dp) :: det
        real(dp), dimension(dim) :: x0, f, SOLVER

        ah = alpha * h

        SOLVER(:)=(F(:) + ah*(yxe(:)*frate(t) + lio*rho_io(t))) &
                                      / (1.D0+ah*(lxe + arate(t)))

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
