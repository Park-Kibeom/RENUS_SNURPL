!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  NEXT_STEP_SDBDF
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE NEXT_STEP_SDBDF

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          REAL(dp) :: Ri0, coef

          INTEGER :: i, j, k

!          write(*,*) 'SDBDF'

          ZP = 0

          DO i = 1, nq + 1
              DO j = 1, nq + 1

                  ZP(i,:) = ZP(i,:) + PASC(i,j) * ZC1(j,:)

              END DO
          END DO

          Ri0 = Ri(nq,1)

          VEC0 = ZP(1,:) - Ri0/Ri(nq,3) * ZP(3,:)                          !  f = zp_[0] - r[0]*zp_[1]

          coef = (Ri0/Ri(nq,3)/2.)**0.5

          ZC(1,:) = SOLVER(ZC1(1,:), time+h,-coef, h, VEC0(:), dim) 

          VEC0(:) = SOLVER(ZC1(1,:), time+h, coef, h, ZC(1,:), dim) 
          ZC(1,:) = VEC0(:)

          VEC0    = FUNC(ZC(1,:), time+h, dim)                !  y1 = h * self.fun(y0, t+h)
          ZC(2,:) = h * VEC0
          ZC(3,:) = 0.5D0* h**2 * FUNC(VEC0, time+h, dim)                !  y1 = h * self.fun(y0, t+h)

          VEC0 = ZC(3,:) - ZP(3,:)

!          ZC(2,:) = ZP(2,:) + Ri(nq,2)/Ri(nq,3) * VEC0

          DO k = 4, nq+1
              ZC(k,:) = ZP(k,:) + Ri(nq,k)/Ri(nq,3) * VEC0          ! zc_ = np.array([y0, y1] + [zp_[i] + r[i] * (y1 - zp_[1]) for i in range(2,len(r))])
          END DO

      END SUBROUTINE NEXT_STEP_SDBDF
