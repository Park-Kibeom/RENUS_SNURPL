!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  NEXT_STEP_BDF
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE NEXT_STEP_BDF
      !DEC$ ATTRIBUTES DLLEXPORT :: NEXT_STEP_BDF
      
          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          REAL(dp) :: Ri0, Y(dim)

          INTEGER :: i, j, k

          ZP = 0

          DO i = 1, nq + 1
              DO j = 1, nq + 1

                  ZP(i,:) = ZP(i,:) + PASC(i,j) * ZC1(j,:)

              END DO
          END DO

          Ri0 = Ri(nq,1)

          VEC0 = ZP(1,:) - Ri0 * ZP(2,:)                          !  f = zp_[0] - r[0]*zp_[1]

          Y(:) = ZP(1,:)
          ZC(1,:) = SOLVER(Y(:), time+h, Ri0, h, VEC0, dim)     !  y0 = self.solver(z[0], t+h, h*r[0], f)

          ZC(2,:) = (ZC(1,:) - VEC0(:))/Ri0
!          ZC(2,:) = h*FUNC(ZC(1,:), time+h, dim)
!          k=2532
!          write(*,*) 't, h: ', time, h
!          write(*,*) 'next: ', dim, ZC(1,k), ZC1(1,k), VEC0(k)

          VEC0 = ZC(2,:) - ZP(2,:)

          DO k = 3, nq+1
              ZC(k,:) = ZP(k,:) + Ri(nq,k) * VEC0          ! zc_ = np.array([y0, y1] + [zp_[i] + r[i] * (y1 - zp_[1]) for i in range(2,len(r))])
          END DO

      END SUBROUTINE NEXT_STEP_BDF

