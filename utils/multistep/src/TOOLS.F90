!
!   PASCAL(P, n)        : расчитываем матрицу Паскаля ранга n
!   FACTORIAL(n)        : факториал числа n
!   NORM(X, dim, p)     : p-норма вектора X
!   TRANS(X, n, x0)     : передвигает элементы массива на одну позицию
!   ANY(X, n)           : проверяет массив на наличие True
!

      SUBROUTINE PASCAL(P,n)

          USE PRECISION, ONLY : dp

          IMPLICIT NONE

          INTEGER :: n, i, j
          INTEGER :: P(n,n)

          P = 0

          DO i = 1, n
              P(1,i) = 1
          END DO

          DO i = 2, n
              DO j = 2, n

                  P(i,j) = P(i-1, j-1) + P(i, j-1) 

              END DO  
          END DO


!          DO j = 1, n
!               write(*,'(40I4)') (P(i,j), i=1,n)
!          END DO



      END SUBROUTINE PASCAL



      FUNCTION FACTORIAL(n) RESULT(fac)

          USE PRECISION, ONLY : dp

          IMPLICIT NONE

          INTEGER :: fac, n, i

          fac = 1

          DO i = 2, n
              fac = fac * i
          END DO

      RETURN
      END FUNCTION FACTORIAL


      SUBROUTINE TRANS(X, n, x0)
      ! Прибавляет элемент x0 к массиву X = (x1, x2, ..., xn),
      ! сохраняя количество его элементов
      ! Результат: X = (x2, x3,..., xn, x0)

          USE PRECISION, ONLY : dp

          IMPLICIT NONE

          INTEGER :: n, i
          LOGICAL :: X(n), x0

          DO i = 1, n-1
              X(i) = X(i+1)
          END DO

          X(n) = X0

      END SUBROUTINE



      FUNCTION ANY(X, n) RESULT(Y)
      ! Логическая функция, 
      ! проверяет массив на наличие элемента True

          USE PRECISION, ONLY : dp

          IMPLICIT NONE

          INTEGER :: i, n
          LOGICAL :: X(n), Y

    
          Y = .False.
          DO i = 1, n
              IF (X(i)) THEN 
                  Y = .True.
                  EXIT
              END IF
          END DO

      END FUNCTION ANY


      FUNCTION NORM(X, dim, p)
      !
      ! Вычисление нормы вектора 
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



      FUNCTION ROUND(x)
      !
      !
          USE PRECISION, ONLY : dp

          IMPLICIT NONE

          INTEGER :: ROUND 
          REAL(dp) :: x
          ROUND = INT(x+0.5)

      END FUNCTION ROUND



      SUBROUTINE RiMATRIX(A, dim)
      !
      !
          USE PRECISION, ONLY : dp

          IMPLICIT NONE

          INTEGER :: n, k, i, j, dim
          REAL(dp) :: A(dim, dim+1)
          REAL(dp), ALLOCATABLE :: B(:,:), X(:)

          A = 0.0

          DO k = 1, dim

              ALLOCATE(B(k, k), X(k))

              DO n = 1, k
                  X(n) = float(n)
                  B(n, 1) = 1.
                  DO j = 1, k-1
                      B(n, j+1) = (-n)**(j + 1)  
                  END DO
              END DO

              IF (k > 1) CALL LAE_SOLVE(B, k, X)

              A(k,1)   = X(1)
              A(k,2)   = 1.
              A(k,3:k+1) = X(2:)

              DEALLOCATE(B, X)

          END DO

!          DO j = 1, dim
!               write(*,'(40ES18.4)') (A(i,j), i=1, dim-1)
!          END DO

      END SUBROUTINE RiMATRIX


      SUBROUTINE LAE_SOLVE(A1, N, B1)

      USE PRECISION, ONLY : dp
      implicit none 
!=====================================================================*
! Solution of the equation A1 * x = B1 Using LU Decomposition         * 
! Last Update              Slava (c) April 1998                       *
!=====================================================================*
! Input: A1 - Matrix, N - Dimension, B1 - Right Side
      integer N
      real(dp) A1(N,N),B1(N)
! Output: A1 - LU Decomposition of the MAtrix A1, B - Solution
! Local Variables: 
      integer k, i, j
      REAL(dp) sum 

! LU Decomposition
      do j = 1,N
! first part beta_ij expression (2.3.12)
         do i= 1, j
            sum = a1(i,j)
            do k = 1, i-1
               sum = sum - a1(i,k)*a1(k,j)
            end do
            a1(i,j) = sum
         end do
! second part alfa_ij expression (2.3.13)
       do i = j+1,N
            sum = a1(i,j)
            do k = 1, j-1
               sum = sum - a1(i,k)*a1(k,j)
            end do
         a1(i,j) = sum/a1(j,j)
       end do
      end do

! Solution LU x = b

! forward substitution
      do i = 1,N
         sum = b1(i)
         do j = 1, i-1
            sum = sum - a1(i,j)*b1(j)
         end do
         b1(i) = sum
       end do

! backsubstitution

      do i = N, 1,-1
         sum = b1(i)
         do j = i+1,N
            sum = sum - a1(i,j)*b1(j)
         end do
         b1(i)=sum/a1(i,i)
      end do

      return
      END SUBROUTINE LAE_SOLVE

