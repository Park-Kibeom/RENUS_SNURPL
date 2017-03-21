     
      subroutine LSODA_INTEGRATE
      !DEC$ ATTRIBUTES DLLEXPORT :: LSODA_INTEGRATE
      USE PRECISION, ONLY : dp
      USE PARS_BDF
      IMPLICIT NONE

!      external lsoda_func, lsoda_jac

      interface 
      subroutine lsoda_func (neq, t, y, ydot)
      integer neq
      double precision t, y(neq), ydot(neq)
      end subroutine lsoda_func


      subroutine lsoda_jac (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y(neq), pd(nrowpd,neq)
      end subroutine lsoda_jac
      end interface


      integer i, iopar, iopt, iout, istate, itask, itol
      integer leniw, lenrw, liw, lout, lrw, mband, meth, mf, miter
      integer ml, mu, neq, nerr, nfe, nfea, nje, nout, nqu, nst
      double precision atol, dtout, er, erm, ero, hu, rtol, t
      double precision rwork(697)
      integer :: iwork(45)

      itol = 1
      rtol = 1.D-4
      atol = 1.D-12
      neq = dim ! dimension of the task      

      lrw = 22 +  9*NEQ + NEQ**2
      liw = 20 + NEQ
      iopt = 1
      mf = 21

      itask = 1
      istate = 1
      iwork(5) = 5 ! order of bdf
      IWORK(6) = 10000 ! number of time steps      

      write(*,*) 'lsoda run... time = ', time, ZC1(1,:), dim, tstop

      call dlsode(lsoda_func,neq,ZC1(1,:),time,
     & tstop,itol,rtol,atol,itask,
     & istate,iopt,rwork,lrw,iwork,liw,lsoda_jac,mf)

      h = rwork(11)
      h1 = h
      nq = iwork(14)
      lenrw = iwork(17)

      write(*,*) 'lsoda complete!'


      end subroutine LSODA_INTEGRATE


      subroutine lsoda_func (neq, t, y, ydot)
      USE PARS_BDF
      USE PRECISION, ONLY : dp            
      integer neq
      real(dp) :: t, y(neq), ydot(neq)
      ydot = func(y, t, neq)
      return
      end subroutine lsoda_func


      subroutine lsoda_jac (neq, t, y, ml, mu, pd, nrowpd)
      USE PARS_BDF
      USE PRECISION, ONLY : dp      
      integer :: neq, ml, mu, nrowpd
      real(dp) :: y(neq), pd(nrowpd,neq), t          
      call jac(y,t,neq,pd)
      return
      end subroutine lsoda_jac

