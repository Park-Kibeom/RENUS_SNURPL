      SUBROUTINE INIT_STEP_SIZE
      !DEC$ ATTRIBUTES DLLEXPORT :: INIT_STEP_SIZE

          USE PRECISION, ONLY : dp

          USE PARS_BDF
    
          IMPLICIT NONE

          REAL(dp) :: Y0(dim), Y1(dim)
          REAL(dp) :: C0, hp1, hp2, er


          Y0(:) = FUNC(ZC1(1,:), time, dim)
          er = LOCERR(Y0(:), ZC1(1,:), dim)
          hp1 = 1. / (1./(tSTOP)**2 + er**2)**0.5

          C0 = 1.
          Y0(:) = SOLVER(ZC1(1,:), time+hp1, C0, hp1, ZC1(1,:), dim)

          Y1(:) = FUNC(Y0(:), time+hp1, dim)
          er = LOCERR(Y1(:), Y0(:), dim)
          hp2 = 1. / (1./(tSTOP)**2 + er**2)**0.5

          h = min(hp1, hp2)

          CALL RENORM(h/h1)
          h1 = h

      END SUBROUTINE INIT_STEP_SIZE
