!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GET_RESULT(YOUT, HOLD, HNEW, TOUT, RES)    : получение результатов расчета
!  SET_SYSTEM(dim_)                           : выделение памяти для массивов
!  INTEGRATE(h_new, FUNC, SOLVER, YOUT, res)  : интегрирование системы на один шаг
!  DEALLOC()                                  : очищение внутренней памяти
!  SET_INTEGRATOR(method)                     : задает метод интегрирования
!  SET_INIT_STATE(Y0_, t0_)                   : задает начальные условия по времени
!  SAVE_STATE                                 : сохраняет текущее состояние
!  PRINT_STATE                                : печатает текущее состояние в консоль
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GET_RESULT(YOUT, HOLD, HNEW, TOUT, RES)
      !DEC$ ATTRIBUTES DLLEXPORT :: GET_RESULT
      ! Программа для получения результатов расчетов:
      ! YOUT : вектор-решение задачи (первый элемент вектора Нордсика)
      ! HOLD : последний временной шаг
      ! HNEW : следующий временной шаг
      ! 

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          REAL(dp):: HOLD, HNEW, TOUT
          REAL(dp):: YOUT(dim)
          INTEGER :: RES

          YOUT(:) = ZC1(1,:)
          HOLD = h2
          HNEW = h
          TOUT = time
          RES = 0
          
          IF(fREJ) RES = 1

      END SUBROUTINE GET_RESULT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SET_PARAMETER(par, val)
      !DEC$ ATTRIBUTES DLLEXPORT :: SET_PARAMETER
      !DEC$ ATTRIBUTES REFERENCE :: par, val

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          CHARACTER*255 :: par, val

!          write(*,*) 'par: ', trim(par)
!          write(*,*) 'val: ', trim(val)

          par = trim(par)
          val = trim(val)

          SELECT CASE (par)

              CASE('InitStepSize') 
                  READ(val,'(F20.10)') h0

              CASE('SetStepSize') 
                  READ(val,'(F20.10)') h
                  CALL RENORM(h/h1)
                  h1 = h

              CASE('Method') 
                  IF (val == "LSODA") NMETHOD = -1
                  IF (val == "IEDS")  NMETHOD = 0
                  IF (val == "BDF")   NMETHOD = 1
                  IF (val == "SDBDF") NMETHOD = 2

              CASE('MaxOrder') 
                  READ(val,'(I10)') NQMAX

              CASE('MaxStepSize') 
                  READ(val,'(F20.10)') HMAX

              CASE('EndOfTime') 
                  READ(val,'(F20.10)') tSTOP

              CASE('MinStepSize') 
                  READ(val,'(F20.10)') HMIN

              CASE('StepSizeHold') 
                  READ(val,'(L10)') HCONST

              CASE('OrderHold') 
                  READ(val,'(L10)') QMaxHold

              CASE('ELocalControl') 
                  READ(val,'(L10)') ELocCtrl

              CASE('MaxStepsNumber') 
                  READ(val,'(I10)') NHMAX

              CASE('MaxRejectSteps') 
                  READ(val,'(I10)') NRMAX

              CASE('LogFileName')
                  if(val /= "None")then
                      fprnt = .True.
                      OPEN(1, FILE=val, STATUS='UNKNOWN',
     &                   ACTION="WRITE")                        
                  endif

              CASE('RejectFileName')
                  if(val /= "None")then
                      freject = .True.
                      OPEN(3, FILE=val, STATUS='UNKNOWN', 
     &                   ACTION="WRITE")
                  endif

              CASE('CheckFileName')
                  if(val /= "None")then
                      f_check_file = .True.
                      OPEN(2, FILE=val, STATUS='UNKNOWN', 
     &                   ACTION="WRITE")                        
                  endif

              CASE DEFAULT
                  WRITE(*,*) 'Parameter ', par, ' not found!'
                  IF (fprnt) THEN
                      WRITE(1,*) 'Parameter ', par, ' not found!'
                  END IF

           END SELECT

      END SUBROUTINE SET_PARAMETER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get(name, adr, n)
      !DEC$ ATTRIBUTES DLLEXPORT :: GET
        use precision, only : dp
        use pars_bdf
        implicit none
      
        character*255 :: name
        integer :: adr, n

        select case(name)
        case('ZC1(1)')
          adr = LOC(ZC1(1,:))
          n = SIZE(ZC1(1,:))
        case default
          write(*,*) 'Parameter ', name, ' not found!'
        end select

      end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GET_PARAMETER(pname, p)
      !DEC$ ATTRIBUTES DLLEXPORT :: GET_PARAMETER
      !DEC$ ATTRIBUTES REFERENCE :: pname, p

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          CHARACTER*255 :: pname
          REAL(dp):: p

          SELECT CASE (pname)

              CASE('h') 
                  p = h

              CASE('h1') 
                  p = h1

              CASE('t') 
                  p = time

              CASE('hr') 
                  p = hr

              CASE('frej') 
                  p = -1.
                  IF (fREJ) p = 1.

              CASE('StepsNumber') 
                  p = float(nh)

              CASE('RejectStepsNumber') 
                  p = float(nr)

              CASE DEFAULT
                  p = 0
                  write(*,*) 'Parameter ', pname, ' not found!'

          END SELECT

      END SUBROUTINE GET_PARAMETER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE DEALLOC
      !DEC$ ATTRIBUTES DLLEXPORT :: DEALLOC

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          if (allocated(vec0)) then
            DEALLOCATE(VEC0, VEC1)
            DEALLOCATE(EL, EL1)
            DEALLOCATE( ZC, ZP, ZC1, ZP1)
            DEALLOCATE(Ri, PASC, XINIT)
          endif

          IF (fprnt) CLOSE(1)
          IF (f_check_file) CLOSE(2)
          IF (freject) CLOSE(3)

      END SUBROUTINE DEALLOC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE SET_INTEGRATOR(method_name, nqmax_)
      !DEC$ ATTRIBUTES DLLEXPORT :: SET_INTEGRATOR

          USE PRECISION, ONLY : dp

          USE PARS_BDF
    
          IMPLICIT NONE
    
          CHARACTER(*) :: method_name 
          INTEGER :: nqmax_
    
          nqmax  = nqmax_

          IF (method_name == "BDF")   NMETHOD = 1
          IF (method_name == "SDBDF") NMETHOD = 2

      END SUBROUTINE SET_INTEGRATOR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE PRINT_STATE_DETAIL

      USE PARS_BDF

      IMPLICIT NONE

      WRITE(*,'(A, I5, A, ES12.3, A, ES12.3)') 
     & '  N = ', nh, '  T = ', time, '  Rel', rel

      WRITE(*,'(A, I5, A, ES12.3, A, ES12.3, A, ES12.3, A, ES12.3)') 
     & '  Next: q = ', nq1, '  h = ', h, 
     & '  Ru = ', RU, '  Rs = ', RS, '  Rd = ', RD

      WRITE(*,'(A, 5ES14.4)') 'Pred.: ', ZP1(:,1)
      WRITE(*,'(A, 6ES14.4)') 'Corr.: ', ZC1(:,1), EL(1)

      write(*,*) ''

      END SUBROUTINE PRINT_STATE_DETAIL




      SUBROUTINE PRINT_STATE(nfile)
      !DEC$ ATTRIBUTES DLLEXPORT :: PRINT_STATE

          USE PARS_BDF

          IMPLICIT NONE

          integer :: nfile
          real :: res

          IF( (fprnt   .and. (nfile == 1) ) .or.
     &        (freject .and. (nfile == 3) ) ) then

!          if (associated(PT_FUNCT)) then 
          if (PT_FUNCT /= 0)  then 
            res = FUNCT(ZC1(1,:),dim)
          else
            res = ZC1(1,1)
          endif

          WRITE(nfile,'(I5, 3ES18.6,I5,1ES18.6,1ES18.6)') 
     &      nh, time, h2, h, nq, rel, res

          END IF

!          WRITE(nfile,'(I5, 2ES18.6, I5, 4ES18.6,I10,I10,I10,4ES18.6)') 
!     &      nh, time, h, nq, rel, rd, rs, ru, ih, iq, nr, 
!     &      ZC1(1,1), ZC1(2,1), ZP(1,1), ZP(2,1)

      END SUBROUTINE PRINT_STATE




      SUBROUTINE CHECK_INPUT_DATA
      !DEC$ ATTRIBUTES DLLEXPORT :: CHECK_INPUT_DATA

          USE PARS_BDF
    
          IMPLICIT NONE

          IF (f_check_file) THEN
    
              WRITE(2,*) 'Integrate started'
              WRITE(2,*) 'Method no.: ', nmethod
              WRITE(2,*) 'Dimension: ', dim
              WRITE(2,*) 'Max order: ', nqmax
              WRITE(2,*) 'Max num. reject steps: ', nrmax
              WRITE(2,*) 'Max num. steps: ', nhmax
              WRITE(2,*) 'Max step size: ', hmax
              WRITE(2,*) 'Min step size: ', hmin
              WRITE(2,*) 'Start time :  ', time
              WRITE(2,*) 'End time :  ', tSTOP              
              WRITE(2,*) 'Init step size:  ', h0
              WRITE(2,*) 'Order hold:  ', QMaxHold
              WRITE(2,*) 'Hold step size:  ', hconst
              WRITE(2,*) 'Local error control: ', ELocCtrl
              WRITE(2,*) '------------------------------------------'
              WRITE(2,*) 'nh time h nq rel rd rs ru ih iq nr ZC1'
              WRITE(2,*) '------------------------------------------'
    
          END IF

      END SUBROUTINE CHECK_INPUT_DATA



      SUBROUTINE TESTCALL
      !DEC$ ATTRIBUTES DLLEXPORT :: TESTCALL

          USE PARS_BDF
          IMPLICIT NONE

          WRITE(*,*) 'Integrate started'
          WRITE(*,*) 'Method no.: ', nmethod
          WRITE(*,*) 'Dimension: ', dim
          WRITE(*,*) 'Max order: ', nqmax
          WRITE(*,*) 'Max num. reject steps: ', nrmax
          WRITE(*,*) 'Max num. steps: ', nhmax
          WRITE(*,*) 'Max step size: ', hmax
          WRITE(*,*) 'Min step size: ', hmin
          WRITE(*,*) 'Start time :  ', time
          WRITE(*,*) 'End time :  ', tSTOP              
          WRITE(*,*) 'Init step size:  ', h0
          WRITE(*,*) 'Order hold:  ', QMaxHold
          WRITE(*,*) 'Hold step size:  ', hconst
          WRITE(*,*) 'Local error control: ', ELocCtrl

      END SUBROUTINE TESTCALL



!      SUBROUTINE SET_FOUT(FOUT)
!      !DEC$ ATTRIBUTES DLLEXPORT :: SET_FOUT

!          USE PARS_BDF

!          CHARACTER*64 :: FOUT

!          fprnt = .True.

!          OPEN(1, FILE=FOUT, STATUS='UNKNOWN')

!      END SUBROUTINE SET_FOUT



