!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE REJECT
      !DEC$ ATTRIBUTES DLLEXPORT :: REJECT

          USE PRECISION, ONLY : dp
          USE PARS_BDF
          IMPLICIT NONE

          EXTERNAL FACTORIAL

          INTEGER :: FACTORIAL, i

          REAL(dp) :: enorm

!          YP = h*FUNC(ZC1(1,:), time+h, dim)
!          EL = ( ZC(nq+1,:) - YP(:) ) / Ri(nq,nq+1)

          EL = ( ZC(nq+1,:) - ZP(nq+1,:) ) / Ri(nq,nq+1)

          enorm = LOCERR(EL(:), ZC1(1,:), dim)

          REL = 1./(nq+1) * Ri(nq,nq+1) * FACTORIAL(nq) * enorm

          fREJ = .false.

          IF ( (REL > 1.).and.(ELocCtrl) ) fREJ = .true. 

!          write(*,*) 'BDF_LU '
!          write(*,*) 'time: ', time
!          write(*,*) 'h: ', h
!          write(*,*) 'ZC1: ', ZC1(1,1:10)
!          write(*,*) 'EL: ', EL(1:10)
!          write(*,*) 'ZC: ', ZC(1,1:5)          
!          write(*,*) 'ZP: ', ZP(1,1:5)          
!          write(*,*) 'REL: ', rel
!          pause

      END SUBROUTINE REJECT

