      SUBROUTINE SET_DRIVER(dim_, FUNC_, SOLVER_, 
     &   LOCERR_, CALLBACK_)
      !DEC$ ATTRIBUTES DLLEXPORT :: SET_DRIVER

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          INTERFACE 
    
              FUNCTION FUNC_(X, t, dim) !RESULT(Y)
                  USE PRECISION, ONLY : dp    
                  INTEGER :: dim
                  REAL(dp), DIMENSION(dim) :: X, FUNC_
                  REAL(dp), INTENT(IN) :: t
              END FUNCTION FUNC_
    
              FUNCTION SOLVER_(Y, t, alpha, h, F, dim)
                  USE PRECISION, ONLY : dp    
                  INTEGER :: dim
                  REAL(dp), DIMENSION(dim) :: Y, F, SOLVER_
                  REAL(dp), INTENT(IN) :: alpha, t, h
              END FUNCTION SOLVER_
    
    
              FUNCTION LOCERR_(dX, X0, dim)
                  USE PRECISION, ONLY : dp    
                  INTEGER :: dim
                  REAL(dp), INTENT(IN) :: dX(dim), X0(dim)
                  REAL(dp):: LOCERR_
              END FUNCTION LOCERR_
    
   
              SUBROUTINE CALLBACK_
              END SUBROUTINE CALLBACK_
    
          END INTERFACE

          INTEGER :: dim_, NOM

          PT_FUNC      = LOC(FUNC_)
          PT_SOLVER    = LOC(SOLVER_)
          PT_LOCERR    = LOC(LOCERR_)
          PT_CALLBACK  = LOC(CALLBACK_)
          PT_FUNCT = 0

          dim = dim_

          NOM = 10  ! max order in BDF
          ALLOCATE ( VEC0(  dim), VEC1(  dim) ) ! âñïîìîãàòåëüíûå âåêòîðà
          ALLOCATE ( EL  (  dim), EL1 (  dim) ) ! âåêòîð ëîêàëüíûõ îøèáîê
          ALLOCATE ( ZC  (NOM+1,dim), ZC1 (NOM+1,dim) ) ! âåêòîð Íîðäñèêà (êîððåêòîð)
          ALLOCATE ( ZP  (NOM+1,dim), ZP1 (NOM+1,dim) ) ! âåêòîð Íîðäñèêà (ïðåäèêòîð)
          ALLOCATE ( Ri  (NOM,   NOM+1) ) ! 
          ALLOCATE ( PASC(NOM+1, NOM+1) ) ! 
          ALLOCATE ( XINIT(dim) ) ! vector of initial state

          CALL PASCAL(PASC, NOM+1)      ! âû÷èñëÿåì ìàòðèöó Ïàñêàëÿ
          CALL RiMATRIX(Ri, NOM)        ! âû÷èñëÿåì ìàòðèöó 

      END SUBROUTINE SET_DRIVER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SET_FUNCTIONAL(dim_, FUNCT_)
      !DEC$ ATTRIBUTES DLLEXPORT :: SET_FUNCTIONAL

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          INTERFACE 
    
              FUNCTION FUNCT_(X, dim) !RESULT(Y)
                  USE PRECISION, ONLY : dp    
                  INTEGER :: dim
                  REAL(dp), DIMENSION(dim) :: X
                  REAL(dp) :: FUNCT_
              END FUNCTION FUNCT_    
    
          END INTERFACE 

          INTEGER :: dim_, NOM

          PT_FUNCT     = LOC(FUNCT_)

      END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE INTEGRATE
      !DEC$ ATTRIBUTES DLLEXPORT :: INTEGRATE
      ! T1_ : âðåìÿ îêîí÷àíèÿ èíòåãðèðîâàíèÿ

          USE PRECISION, ONLY : dp

          USE PARS_BDF
    
          IMPLICIT NONE

          REAL(dp) :: hsave
          INTEGER :: k, i, NMETHOD_

!          write(*,*) 'INTEGRATE begin '
          NMETHOD_ = NMETHOD
          fstop  = .False.

          DO k = 1, NHMAX
            DO i = 1, 100

                SELECT CASE (NMETHOD)

                CASE(-1)
!                  CALL LSODA_INTEGRATE
                  CALL PRINT_STATE(1)

                CASE(0)
                  CALL CALLBACK
                  CALL NEXT_STEP_IEDS

                CASE(1)
                  CALL CALLBACK
!                  IF (i < 3) THEN
                   CALL NEXT_STEP_BDF
                   IF (NR < NRMAX) CALL REJECT                   
!                    write(*,*) 'bdf step ', i
!                  ELSE
!                    CALL NEXT_STEP_IEDS
!                  ENDIF

                CASE(2)
                  CALL CALLBACK                 
                  CALL NEXT_STEP_SDBDF
                  IF (NR < NRMAX) CALL REJECT

                END SELECT

                IF (fREJ) THEN
                    IF (i == 1) THEN
                        NR = NR + 1                     
                        hsave = h                     
                    END IF 

                    IF ((i == 2) .or.((hsave/h > 1.2)
     &                       .and.(nmethod .ne. 0))) THEN
!                    IF ((hsave/h > 1.5).and.(nmethod .ne. 0)) then

!                    IF (i > 10) then
!                        write(*,*) 'Number of rejects = ', i, h
!                        stop
!                    endif

!                    IF (i == 2) then
                        NMETHOD = 0 
                        nq1 = 1
                        nq = 1
                        CALL RENORM(hsave/h)
                        h = hsave
                        h1 = hsave
!                        h  = 0.1*h
!                        h1 = h
!                        CALL RENORM(0.1)
                        CALL PRINT_STATE(3)                    
                    ELSE 
                        CALL PRINT_STATE(3)                    
                        CALL H_RESIZE                
                    ENDIF

                    IF (fstop) EXIT

                ELSE
                    NMETHOD = NMETHOD_
                    EXIT
                END IF 

            END DO

            CALL NEW_HSIZE
            CALL ACCEPT_TIME_STEP
            CALL PRINT_STATE(1)
            IF (fstop) EXIT
          END DO

      END SUBROUTINE INTEGRATE


      SUBROUTINE SET_INIT_STATE(Y0_, t0_)
      !DEC$ ATTRIBUTES DLLEXPORT :: SET_INIT_STATE

          USE PRECISION, ONLY : dp

          USE PARS_BDF
    
          IMPLICIT NONE

          REAL(dp) :: t0_, Y0_(dim), C0, R1(dim)

!          write(*,*) 'INIT_STATE start..'

          XINIT(:) = Y0_(:)
          ZC1(1,:) = Y0_(:)
!          ZC1(2,:) = h1*FUNC(Y0_, t0_+h1, dim)

          C0 = 1. 
          R1(:) = SOLVER(Y0_(:), t0_, C0, h1, Y0_(:), dim)   
          ZC1(2,:) = R1(:) - ZC1(1,:)

          IF (NMETHOD == 2) THEN
              VEC0 = h1 * FUNC(Y0_, t0_+h1, dim)
              ZC1(3,:) = 0.5 * h1 * FUNC(VEC0, t0_+h1, dim)
          END IF 

          ZC(:,:) = ZC1(:,:)
          ZP(:,:) = ZP1(:,:)

          IF (NMETHOD == 1) NQ1 = 1
          IF (NMETHOD == 2) NQ1 = 2

          NQ    = NQ1
          NQMIN = NQ1

          time = t0_ 
!          write(*,*) 'SET_INIT_STATE, t0 = ', t0_
   
!   !      IF (h <= 0) CALL INIT_H
          EL1 = 0.
          fREJ   = .False.
          fstop  = .False.
          iQ = 1
          iH = 0
          Nh = 0
          Nr = 0  

!          write(*,*) 'INIT_STATE complete!'
          h2 = 1.
          h1 = h0
          h  = h0
          CALL RENORM(h0)
          CALL PRINT_STATE(1)

      END SUBROUTINE SET_INIT_STATE


      SUBROUTINE ACCEPT_TIME_STEP
      !DEC$ ATTRIBUTES DLLEXPORT :: ACCEPT_TIME_STEP

          USE PRECISION, ONLY : dp
          USE PARS_BDF

          IMPLICIT NONE

          CALL RENORM(h/h1)

          ZC1 = ZC
          ZP1 = ZP
          EL1 = EL
          time  = time + h1
          h2 = h1
          h1 = h
          nq1 = nq
          nh = nh + 1

          IF (time >= tstop) fstop = .True. 

      END SUBROUTINE ACCEPT_TIME_STEP      


      SUBROUTINE RESTURT
      !DEC$ ATTRIBUTES DLLEXPORT :: RESTURT
          USE PRECISION, ONLY : dp
          USE PARS_BDF
          IMPLICIT NONE

          real(dp) :: dt
          real(dp) :: dX1(dim), dX2(dim)

          nq = 1

          dt = 0.1 * h1

          dX1 = FUNC(ZC1(1,:), time, dim)
          dX2 = FUNC(ZC1(1,:), time+dt, dim)

          ZC1(3,:) = (0.5*h1**2) * (1./dt) * (dX2(:) - dX1(:))

!          CALL H_RESIZE

!          ZC1(3,:) = 0.h1*FUNC(ZC1(1,:), time+dt, dim)

      END SUBROUTINE  



      SUBROUTINE SET_JAC(JAC_)
      !DEC$ ATTRIBUTES DLLEXPORT :: SET_JAC

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          INTERFACE 
    
              SUBROUTINE JAC_(X, t, dim, PD) !RESULT(Y)
                  USE PRECISION, ONLY : dp    
                  INTEGER :: dim
                  REAL(dp) :: PD(dim,dim),X(dim)
                  REAL(dp), INTENT(IN) :: t
              END SUBROUTINE JAC_
    
          END INTERFACE

          PT_JAC       = LOC(JAC_)          

      END SUBROUTINE