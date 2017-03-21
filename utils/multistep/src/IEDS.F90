!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  Implicit Euler with double-step local error control
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE INTEGRATOR_TEST
      !DEC$ ATTRIBUTES DLLEXPORT :: INTEGRATOR_TEST

          USE PRECISION, ONLY : dp
          USE PARS_BDF
          IMPLICIT NONE

          INTEGER :: k
          REAL(dp) :: Y1(dim), Y2(dim), Y0(dim), Ri0, a
          REAL(dp) :: Y3(dim)

          k=2532
          Y0(:) = ZC1(1,:)
          Ri0 = 1.0


          fREJ = .false.
          CALL CALLBACK

!          Y1(:) = FUNC(Y0(:), time+h, dim)
          Y1(:) = SOLVER(Y0(:), time+h, Ri0, h, Y0(:), dim)          
          write(*,*) 'Y1: ', Y1(k:k+5)

!          Y3(:) = SOLVER(Y0(:), time+h, Ri0, h, Y0(:), dim)          
!          write(*,*) 'Y3: ', Y3(k:k+5)

          fREJ = .true.
          CALL CALLBACK

!          Y2(:) = FUNC(Y0(:), time+h, dim)
          Y2(:) = SOLVER(Y0(:), time+h, Ri0, h, Y0(:), dim)
          write(*,*) 'Y2: ', Y2(k:k+5)

          a = LOCERR(Y2-Y1, Y0(:), dim)

          write(*,*) 'time: ', time
          write(*,*) 'h: ', h
          write(*,*) 'REL: ', a
          pause          

      END SUBROUTINE


      SUBROUTINE NEXT_STEP_IEDS
      !DEC$ ATTRIBUTES DLLEXPORT :: NEXT_STEP_IEDS
      
          USE PRECISION, ONLY : dp
          USE PARS_BDF
          IMPLICIT NONE

          REAL(dp) :: Ri0, Y1(dim), Y2(dim), Y0(dim)
          INTEGER :: i, j, k

          Ri0 = 1.

          Y0(:) = ZC1(1,:)

          Y1(:) = SOLVER(Y0(:), time+h, Ri0, h, Y0(:), dim)    

          fREJ = .true.
          CALL CALLBACK

          Y2(:)   = SOLVER(Y0(:),time+0.5*h,Ri0, 0.5*h,Y0(:),dim)

          ZC(1,:) = SOLVER(Y2(:), time+h, Ri0, 0.5*h, Y2(:), dim)

          ZC(2,:) = h*FUNC(ZC(1,:), time+h, dim)

          EL(:) = Y1(:) - ZC(1,:)

          REL = LOCERR(EL(:), ZC1(1,:), dim)

          IF ( (REL > 1.).and.(ELocCtrl) ) then
            fREJ = .true.
          else
            fREJ = .false.
          endif

!          write(*,*) 'IEDS '
!          write(*,*) 'time: ', time
!          write(*,*) 'h: ', h
!          write(*,*) 'ZC1: ', ZC1(1,k:k+5)
!          write(*,*) 'EL: ', EL(1:10)
!          write(*,*) 'Y1: ', Y1(1:5)
!          write(*,*) 'Y2: ', Y2(1:5)
!          write(*,*) 'Y3: ', Y3(k:k+5)
!          write(*,*) 'ZC: ', ZC(1,k:k+5)          
!          write(*,*) 'ZC2: ', ZC(2,k:k+5)          
!          write(*,*) 'REL: ', rel
!          pause


      END SUBROUTINE NEXT_STEP_IEDS


      SUBROUTINE NEXT_STEP_IEDS_
      !DEC$ ATTRIBUTES DLLEXPORT :: NEXT_STEP_IEDS_
      
          USE PRECISION, ONLY : dp
          USE PARS_BDF
          IMPLICIT NONE

          REAL(dp) :: Ri0, Y1(dim), Y2(dim), Y0(dim)
          INTEGER :: i, j, k

          Ri0 = 1.

          Y0(:) = ZC1(1,:)

          nq = 1
          nq1 = 1

          Y1(:) = SOLVER(Y0(:), time+2*h, Ri0, 2*h, Y0(:), dim)    

          fREJ = .true.
          CALL CALLBACK

          Y2(:)   = SOLVER(Y0(:),time+h,Ri0, h,Y0(:),dim)

          ZC(1,:) = SOLVER(Y2(:), time+2*h, Ri0, h, Y2(:), dim)

          ZC(2,:) = h*FUNC(ZC(1,:), time+2*h, dim)

          EL(:) = Y1(:) - ZC(1,:)

          REL = LOCERR(EL(:), ZC(1,:), dim)

          IF ( (REL > 1.).and.(ELocCtrl) ) then
            fREJ = .true.
          else
            fREJ = .false.
          endif

          write(*,*) 'IEDS '
          write(*,*) 'time: ', time
          write(*,*) 'h: ', h
!          write(*,*) 'ZC1: ', ZC1(1,k:k+5)
!          write(*,*) 'EL: ', EL(1:10)
!          write(*,*) 'Y1: ', Y1(1:5)
!          write(*,*) 'Y2: ', Y2(1:5)
!          write(*,*) 'Y3: ', Y3(k:k+5)
!          write(*,*) 'ZC: ', ZC(1,k:k+5)          
!          write(*,*) 'ZC2: ', ZC(2,k:k+5)          
          write(*,*) 'REL: ', rel
          pause

      END SUBROUTINE NEXT_STEP_IEDS_      

