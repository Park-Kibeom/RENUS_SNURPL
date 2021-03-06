      MODULE PARS_BDF

          USE PRECISION, ONLY : dp

          IMPLICIT NONE

! ��������� ������

          INTEGER :: dim  ! ����������� ������� ODE
          INTEGER :: t0   ! ��������� �����


! ��������� �����������

          INTEGER :: NMETHOD = 1 ! 1 - "BDF", 2 - "SDBDF"
   

! ��������� �������� ��������� ������
  
          INTEGER :: NH = 0  ! ����� ���������� ����
          INTEGER :: nq = 1, nq1 = 1  ! ������� ������� ������ BDF
    
          REAL(dp) :: time = 0. ! ������� �����
          REAL(dp) :: h  = 1.E-6 ! ��������� ��������� ���
          REAL(dp) :: h2 = 1.E-6  ! ���������� ��������� ���
          REAL(dp) :: h1 = 1.E-6  ! ������� ��������� ���
          REAL(dp) :: h0 = 1.E-6 ! ��������� ��������� ���
          REAL(dp) :: RH  ! ������������� ��������� ���� 
          REAL(dp) :: RU  ! ������������� ��������� ���� ��� ��������� ������� BDF
          REAL(dp) :: RS  ! ������������� ��������� ���� ��� ������� ������� BDF
          REAL(dp) :: RD  ! ������������� ��������� ���� ��� ��������� ������� BDF
          REAL(dp) :: REL ! ������������� ��������� ������ �� ������� ����
    
          REAL(dp), ALLOCATABLE  :: XINIT(:) ! initial state

          REAL(dp), ALLOCATABLE  :: ZC(:,:) ! ������ �������� (���������) �� ������� ����
          REAL(dp), ALLOCATABLE  :: ZP(:,:) ! ������ �������� (���������) �� ������� ����
    
          REAL(dp), ALLOCATABLE  :: ZC1(:,:) ! ������ �������� (���������) �� ���������� ����
          REAL(dp), ALLOCATABLE  :: ZP1(:,:) ! ������ �������� (���������) �� ���������� ����
    
          REAL(dp), ALLOCATABLE  :: EL(:)   ! ������ ��������� ������ �� ������� ����
          REAL(dp), ALLOCATABLE  :: EL1(:)  ! ������ ��������� ������ �� ���������� ����


! �������� �������� � �������� ��������� ���ר���

          LOGICAL :: HCONST = .False.  ! HCONST == True - ��������� ��� ����������
          LOGICAL :: QMaxHold = .False.! ���������� ������� �������������� �� �����������-��������� ��������
          LOGICAL :: ELocCtrl = .True. ! ������� �� ��������� �������
          LOGICAL :: fREJ    ! == True - ������ ��� ��������
          LOGICAL :: fstop = .False.  ! == True - ��������� ��������������
          LOGICAL :: fprnt = .False.   ! == True - ������ ������������� �����������
          LOGICAL :: freject = .False.   ! == True - ������ ����������� �����
          LOGICAL :: f_check_file = .False.
    
          INTEGER :: NQMAX = 1  ! ������������ ������� ������
          INTEGER :: NQMIN = 1  ! ����������� ������� ������
          INTEGER :: NHMAX = 10000  ! ������������ ����� ��������� �����
          INTEGER :: NRMAX = 1000  ! ������������ ����� ����������� ��������� �����
    
          REAL(dp):: tSTOP = 1.   ! ��������� ������� ��������������
    
          REAL(dp) :: HMAX = 1.      ! ������������ �������� ���������� ����
          REAL(dp) :: HMIN = 1.E-3   ! �����������  �������� ���������� ����
    
    
! ������� ��������� ��������� ������
    
          INTEGER :: NR = 0 ! ����� ����������� ��������� �����

!          INTEGER :: NFUN ! ����� ������� ������� func
!          INTEGER :: NSOL ! ����� ������� ������� solve
          INTEGER :: iQ ! ������� ��� ��������� ������� BDF
          INTEGER :: iH ! ������� ��� ��������� ���������� ����


! ��������������� (���������) ���������� � ���������

          REAL(dp) :: hr ! ��������� ����������� ��� (��������, ���� �������� �� �����)
          REAL(dp) :: hreject ! ��������� ����������� ���
          REAL(dp), ALLOCATABLE  :: VEC0(:), VEC1(:)

          REAL(dp), ALLOCATABLE :: Ri(:,:)      ! ������������ ������ BDF
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

!      INTEGER, DIMENSION(8,8) :: PASC                  ! ������� �������


! �������, ������������ � BDF

          INTERFACE 
    
              FUNCTION FUNC(X, t, dim) !RESULT(Y)
                  ! ������ ����� ���������
                  USE PRECISION, ONLY : dp
                  INTEGER :: dim
                  REAL(dp), DIMENSION(dim) :: X, FUNC
                  REAL(dp), INTENT(IN) :: t
    
              END FUNCTION FUNC

              FUNCTION FUNCT(X, dim) !RESULT(Y)
                  ! ������ ����� ���������
                  USE PRECISION, ONLY : dp
                  INTEGER :: dim
                  REAL(dp), DIMENSION(dim) :: X
                  REAL(dp) :: FUNCT
    
              END FUNCTION FUNCT
    
              SUBROUTINE JAC(X, t, dim, PD) !RESULT(Y)
                  ! ������ ����� ���������
                  USE PRECISION, ONLY : dp
                  INTEGER :: dim
                  REAL(dp) :: PD(dim,dim), X(dim)
                  REAL(dp), INTENT(IN) :: t
    
              END SUBROUTINE JAC
    
              FUNCTION SOLVER(Y, t, alpha, h, F, dim)
                  ! ������� ��������� ������ 
                  USE PRECISION, ONLY : dp
                  INTEGER :: dim
                  REAL(dp), DIMENSION(dim) :: Y, F, SOLVER
                  REAL(dp), INTENT(IN) :: alpha, t, h
    
              END FUNCTION SOLVER
    
    
              FUNCTION LOCERR(dX, X0, dim)
                  ! ������ ��������� ������
                  USE PRECISION, ONLY : dp
                  INTEGER :: dim
                  REAL(dp), INTENT(IN) :: dX(dim), X0(dim)
                  REAL(dp):: LOCERR
    
              END FUNCTION LOCERR
    
    
              SUBROUTINE CALLBACK
                  ! �������, ���������� ����� ��������� �� ��������� ���
              END SUBROUTINE CALLBACK

          END INTERFACE
    
          POINTER(PT_FUNC,      FUNC)
          POINTER(PT_JAC,       JAC)          
          POINTER(PT_SOLVER,    SOLVER)
          POINTER(PT_LOCERR,    LOCERR)
          POINTER(PT_CALLBACK,  CALLBACK)
          POINTER(PT_FUNCT,     FUNCT)

      END MODULE PARS_BDF


