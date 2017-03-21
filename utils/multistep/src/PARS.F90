      MODULE PARS_BDF

          USE PRECISION, ONLY : dp

          IMPLICIT NONE

! ПАРАМЕТРЫ МОДЕЛИ

          INTEGER :: dim  ! размерность системы ODE
          INTEGER :: t0   ! начальное время


! ПАРАМЕТРЫ ИНТЕГРАТОРА

          INTEGER :: NMETHOD = 1 ! 1 - "BDF", 2 - "SDBDF"
   

! ПАРАМЕТРЫ ТЕКУЩЕГО СОСТОЯНИЯ МОДЕЛИ
  
          INTEGER :: NH = 0  ! номер временного шага
          INTEGER :: nq = 1, nq1 = 1  ! текущий порядок метода BDF
    
          REAL(dp) :: time = 0. ! текущее время
          REAL(dp) :: h  = 1.E-6 ! следующий временной шаг
          REAL(dp) :: h2 = 1.E-6  ! предыдущий временной шаг
          REAL(dp) :: h1 = 1.E-6  ! текущий временной шаг
          REAL(dp) :: h0 = 1.E-6 ! начальный временной шаг
          REAL(dp) :: RH  ! относительное изменение шага 
          REAL(dp) :: RU  ! относительное изменение шага при повышении порядка BDF
          REAL(dp) :: RS  ! относительное изменение шага при прежнем порядке BDF
          REAL(dp) :: RD  ! относительное изменение шага при понижении порядка BDF
          REAL(dp) :: REL ! относительная локальная ошибка на текущем шаге
    
          REAL(dp), ALLOCATABLE  :: XINIT(:) ! initial state

          REAL(dp), ALLOCATABLE  :: ZC(:,:) ! вектор Нордсика (корректор) на текущем шаге
          REAL(dp), ALLOCATABLE  :: ZP(:,:) ! вектор Нордсика (предиктор) на текущем шаге
    
          REAL(dp), ALLOCATABLE  :: ZC1(:,:) ! вектор Нордсика (корректор) на предыдущем шаге
          REAL(dp), ALLOCATABLE  :: ZP1(:,:) ! вектор Нордсика (предиктор) на предыдущем шаге
    
          REAL(dp), ALLOCATABLE  :: EL(:)   ! вектор локальной ошибки на текущем шаге
          REAL(dp), ALLOCATABLE  :: EL1(:)  ! вектор локальная ошибки на предыдущем шаге


! КОНТРОЛЬ ТОЧНОСТИ И КРИТЕРИИ ОСТАНОВКИ РАСЧЁТОВ

          LOGICAL :: HCONST = .False.  ! HCONST == True - временной шаг постоянный
          LOGICAL :: QMaxHold = .False.! удерживать порядок интегрирования на максимально-возможном значении
          LOGICAL :: ELocCtrl = .True. ! следить за локальной ошибкой
          LOGICAL :: fREJ    ! == True - данный шаг отброшен
          LOGICAL :: fstop = .False.  ! == True - остановка интегрирования
          LOGICAL :: fprnt = .False.   ! == True - печать промежуточных результатов
          LOGICAL :: freject = .False.   ! == True - печать отброшенных шагов
          LOGICAL :: f_check_file = .False.
    
          INTEGER :: NQMAX = 1  ! максимальный порядок метода
          INTEGER :: NQMIN = 1  ! минимальный порядок метода
          INTEGER :: NHMAX = 10000  ! максимальное число временных шагов
          INTEGER :: NRMAX = 1000  ! максимальное число отброшенных временных шагов
    
          REAL(dp):: tSTOP = 1.   ! окончание времени интегрирования
    
          REAL(dp) :: HMAX = 1.      ! максимальная величина временного шага
          REAL(dp) :: HMIN = 1.E-3   ! минимальная  величина временного шага
    
    
! ИСТОРИЯ ИЗМЕНЕНИЯ СОСТОЯНИЯ МОДЕЛИ
    
          INTEGER :: NR = 0 ! число отброшенных временных шагов

!          INTEGER :: NFUN ! число вызовов функции func
!          INTEGER :: NSOL ! число вызовов функции solve
          INTEGER :: iQ ! счетчик для изменения порядка BDF
          INTEGER :: iH ! счетчик для изменения временного шага


! ВСПОМОГАТЕЛЬНЫЕ (ВРЕМЕННЫЕ) ПЕРЕМЕННЫЕ и ПАРАМЕТРЫ

          REAL(dp) :: hr ! последний отброшенный шаг (возможно, этот параметр не нужен)
          REAL(dp) :: hreject ! последний отброшенный шаг
          REAL(dp), ALLOCATABLE  :: VEC0(:), VEC1(:)

          REAL(dp), ALLOCATABLE :: Ri(:,:)      ! коэффициенты метода BDF
          INTEGER, ALLOCATABLE :: PASC(:,:)

!          REAL(dp), DIMENSION(6,7) :: Ri = (/ 
!     & 1.,  2./3, 6./11, 24./50, 120./274,  720./1764,  
!     & 1.,  1.,   1.,     1.,      1.,      1.,       
!     & 0.,  1./3, 6./11, 35./50, 225./274, 1624./1764, 
!     & 0.,  0.,   1./11, 10./50,  85./274,  735./1764,
!     & 0.,  0.,   0.,     1./50,  15./274,  175./1764,  
!     & 0.,  0.,   0.,     0.,      1./274,   21./1764,   
!     & 0.,  0.,   0.,     0.,      0.,        1./1764/)
!

!          REAL(dp), DIMENSION(7,8) :: Ri = (/ 
!     & 1.,  2./3, 6./11, 24./50, 120./274,  720./1764,  140./363,   
!     & 1.,  1.,   1.,     1.,      1.,      1.,           1.,
!     & 0.,  1./3, 6./11, 35./50, 225./274, 1624./1764, 3283./3267, 
!     & 0.,  0.,   1./11, 10./50,  85./274,  735./1764, 6769./13068,
!     & 0.,  0.,   0.,     1./50,  15./274,  175./1764,  490./3267,  
!     & 0.,  0.,   0.,     0.,      1./274,   21./1764,  161./6534,   
!     & 0.,  0.,   0.,     0.,      0.,        1./1764,    7./3267,    
!     & 0.,  0.,   0.,     0.,      0.,        0.,         1./13068/)

!      INTEGER, DIMENSION(8,8) :: PASC                  ! матрица Паскаля


! ФУНКЦИИ, ПЕРЕДАВАЕМЫЕ В BDF

          INTERFACE 
    
              FUNCTION FUNC(X, t, dim) !RESULT(Y)
                  ! правая часть уравнения
                  USE PRECISION, ONLY : dp
                  INTEGER :: dim
                  REAL(dp), DIMENSION(dim) :: X, FUNC
                  REAL(dp), INTENT(IN) :: t
    
              END FUNCTION FUNC

              FUNCTION FUNCT(X, dim) !RESULT(Y)
                  ! правая часть уравнения
                  USE PRECISION, ONLY : dp
                  INTEGER :: dim
                  REAL(dp), DIMENSION(dim) :: X
                  REAL(dp) :: FUNCT
    
              END FUNCTION FUNCT
    
              SUBROUTINE JAC(X, t, dim, PD) !RESULT(Y)
                  ! правая часть уравнения
                  USE PRECISION, ONLY : dp
                  INTEGER :: dim
                  REAL(dp) :: PD(dim,dim), X(dim)
                  REAL(dp), INTENT(IN) :: t
    
              END SUBROUTINE JAC
    
              FUNCTION SOLVER(Y, t, alpha, h, F, dim)
                  ! решение матричной задачи 
                  USE PRECISION, ONLY : dp
                  INTEGER :: dim
                  REAL(dp), DIMENSION(dim) :: Y, F, SOLVER
                  REAL(dp), INTENT(IN) :: alpha, t, h
    
              END FUNCTION SOLVER
    
    
              FUNCTION LOCERR(dX, X0, dim)
                  ! оценка локальной ошибки
                  USE PRECISION, ONLY : dp
                  INTEGER :: dim
                  REAL(dp), INTENT(IN) :: dX(dim), X0(dim)
                  REAL(dp):: LOCERR
    
              END FUNCTION LOCERR
    
    
              SUBROUTINE CALLBACK
                  ! функция, вызываемая перед переходом на следующий шаг
              END SUBROUTINE CALLBACK

          END INTERFACE
    
          POINTER(PT_FUNC,      FUNC)
          POINTER(PT_JAC,       JAC)          
          POINTER(PT_SOLVER,    SOLVER)
          POINTER(PT_LOCERR,    LOCERR)
          POINTER(PT_CALLBACK,  CALLBACK)
          POINTER(PT_FUNCT,     FUNCT)

      END MODULE PARS_BDF


