!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  GET_RD
!  GET_RU
!  GET_RS
!
!  RENORM
!  NEW_HSIZE
!  INIT_H
!  RESTRICT(X, flag)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GET_RD

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          INTEGER, EXTERNAL :: FACTORIAL

          REAL(dp) :: enorm, tau, Dq, Y1(dim), Y2(dim)


          Y1(1:dim) = ZC1(1,1:dim)
          Y2(1:dim) = ZC(nq+1,1:dim)

          enorm = LOCERR(Y2(:), Y1(:), dim)

          tau = nq / REAL(FACTORIAL(nq))

          Dq = enorm / tau

          rd = 1. / (Dq**(1./nq) + 1.E-6) / 1.3

      END SUBROUTINE GET_RD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GET_RU

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          EXTERNAL FACTORIAL

          INTEGER :: FACTORIAL

          REAL(dp) :: enorm, tau, Dq, Y1(dim), Y2(dim)

          Y1(:) = EL1(:) - EL(:)
          Y2(:) = ZC1(1,:)
          enorm = LOCERR(Y1(:), Y2(:), dim)

          tau = (nq+2)/Ri(nq,nq+1) / FACTORIAL(nq)

          Dq = enorm / tau

          ru = 1. / (Dq**(1./(nq+2)) + 1.E-6) / 1.4

      END SUBROUTINE GET_RU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE GET_RS

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          INTEGER, EXTERNAL :: FACTORIAL

          REAL(dp) :: enorm, tau, Dq, Y(dim)

          Y(:) = ZC1(1,:)
          enorm = LOCERR(EL(:), Y(:), dim)

          tau = (nq+1)/Ri(nq,nq+1) / FACTORIAL(nq)

          Dq = enorm / tau

          rs = 1. / (Dq**(1./(nq+1)) + 1.E-6) / 1.2

          if (nmethod == 0) then
              rs = max(min(0.8 * 1.0/(rel ** 0.5), 2.0), 0.1)
!              write(*,*) 'rs = ', rs
!              write(*,*) 'rel = ', rel
!              write(*,*) 'enorm = ', enorm
          endif

!          write(*,*) 'rs = ', rs        

!          REL = 1./(nq+1) * Ri(nq,nq+1) * FACTORIAL(nq) * enorm

      END SUBROUTINE GET_RS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE NEW_HSIZE
      !DEC$ ATTRIBUTES DLLEXPORT :: NEW_HSIZE

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          CHARACTER :: cs

          REAL(dp), EXTERNAL :: RESTRICT, ROUND
          REAL(dp) :: r_
          INTEGER :: q_

          CALL GET_RD
          CALL GET_RS
          CALL GET_RU

          cs = ''

          IF ((NH > 0).and.(iH <= 0).and.(.not.(fREJ)) ) THEN

              IF (iQ <= 0) THEN

                  IF ((nq1 == nqmin).and.(nqmax == nqmin)) cs = 's'

                  IF ((nq1 == nqmin).and.(nqmax > nqmin)) THEN
                      IF (rs >= ru) cs = 's'
                      IF (rs <  ru) cs = 'u'
                      IF (QMaxHold) cs = 'u'          
                  END IF
    
                  IF ((nq1 == nqmax).and.(nq1 > nqmin)) THEN
                      IF (rs >= rd) cs = 's'
                      IF (rs <  rd) cs = 'd'
                      IF (QMaxHold) cs = 's'          
                  END IF
    
                  IF ((nq1 < nqmax).and.(nq1 > nqmin)) THEN
                      IF ((rd >= ru).and.(rd >= rs)) cs = 'd'
                      IF ((rs >= ru).and.(rs >= rd)) cs = 's'
                      IF ((ru >= rs).and.(ru >= rd)) cs = 'u'
                      IF (QMaxHold) cs = 'u'          
                  END IF

              ELSE

                  cs = 's'

              END IF

          END IF

                              
          SELECT CASE(cs)

              CASE ('d') 
                  q_ = - 1
                  r_ = rd
                  iQ = nq !max
    
              CASE ('s') 
                  q_ = 0
                  r_ = rs
                  iQ = iQ - 1
    
              CASE ('u') 
                  q_ = 1
                  r_ = ru
                  iQ = nq !max

              CASE DEFAULT
                  q_ = 0  ! изменение порядка метода 
                  r_ = 1. ! относительное изменение шага
                  iQ = iQ - 1

          END SELECT

!          r_ = r_**0.4 * (h2/h1)**0.6

          r_ = RESTRICT(r_, time+h1, tSTOP, h1)

          IF (HCONST) r_ = 1

          IF ( ABS(r_ - 1) < 0.01 ) THEN
              iH = iH - 1
!              r_ = 1.   ! не нужно раскомментировать, т.к. это противоречит "подстройке" временных шагов к конечной временной точке !!! 
                         ! программа будет вылетать !!
          ELSE
              iH = nq !max
          END IF

          h = h1 * r_
          nq = nq1 + q_
          rh = r_

          IF (q_ > 0) THEN
              ZC(nq+1,:) = EL(:) * Ri(nq-1, nq) / nq 
          END IF

      END SUBROUTINE NEW_HSIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE H_RESIZE
      !DEC$ ATTRIBUTES DLLEXPORT :: H_RESIZE
      
          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          EXTERNAL RESTRICT

          REAL(dp) :: RESTRICT
          REAL(dp) :: r_

          CALL GET_RS

          r_ = rs

          h = h * r_
          rh = r_
          hr = h1
          h1 = h
          CALL RENORM(r_)

      END SUBROUTINE H_RESIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE INIT_H

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

!          INTEGER :: i
!          REAL(8) :: tol, w0, w1

!          DO i = 1, dim
!              IF (ZC1(1,i) .NE. 0.) VEC0(i) = ATOL(i) / ZC1(1,i)
!          END DO

!          tol = MAX(VEC0, dim)
!          VEC1(:) = tol

!          VEC0 = LOCERR(VEC1(:), ZC1(1,:), dim)

!          DO i = 1, dim
!              w1 = w1 + ( ZC1(2,i) * VEC0(i) ) ** 2
!          END DO

!          h = (tol / (1./w0**2 + w1/dim) ) ** 0.5

!          h = MIN(h, hmax)

      END SUBROUTINE INIT_H

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE RENORM(dr)
      ! Перенормировка вектора Нордсика ZC(nq,n) и ZP(nq,n)
      ! на новый временной шаг
      ! r = h_new / h_old

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          INTEGER :: q
          REAL(dp) :: dr, a

          IF (.not.(dr == 1)) THEN
              DO q = 2, nq + 1
                  a = dr**(q-1)
                  ZC (q,:) = ZC (q,:) * a
                  ZP (q,:) = ZP (q,:) * a
                  ZC1(q,:) = ZC1(q,:) * a
                  ZP1(q,:) = ZP1(q,:) * a
              END DO
          END IF

      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION RESTRICT(x, t, t1, dt) RESULT(y)

          USE PRECISION, ONLY : dp

          USE PARS_BDF

          IMPLICIT NONE

          EXTERNAL ROUND

          INTEGER :: ROUND
          INTEGER :: num
          REAL(dp) :: x, y, t, t1, dt, tnext, a

          y = x
          IF (x > 2000.) y = 2000.

          if (y > HMAX/dt) y = HMAX/dt
          if (y < HMIN/dt) y = HMIN/dt

          tnext = t + y*dt

          IF (tnext > t1) THEN
              y = (t1 - t) / dt
          ELSE
              a = (t1 - t) / (y*dt) - 1.D-14
!              IF (a < 1.E+3) THEN
!                  num = CEILING ( a )
!                  y = (t1 - t) / dt / num
!              END IF
          END IF

          if (y < 0) then
              write(*,*) 'ERROR: negative step size, h = ', y*dt
              write(*,*) 'x, t, t1, dt = ', x, t, t1, dt
              y = 0
              fstop = .True.
          endif

!          write(*,*) 'current time, max time: ', t, t1
!          write(*,*) '1) hnew, tnext: ', y*dt, y*dt+t
!          if(y*dt+t > t1) y = (t1-t)/dt
!          write(*,*) '2) hnew, tnext: ', y*dt, y*dt+t

      RETURN
      END FUNCTION RESTRICT





