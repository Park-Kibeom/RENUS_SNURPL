    module control

      logical :: ifinit = .False.

      contains
    end module control

    subroutine xenon_transient(dt,yio,yxe,absxe,xe0,io0,lio,lxe &
      ,xeout)
      !DEC$ ATTRIBUTES DLLEXPORT :: xenon_transient
      !
      !
      use precision, only : dp
      use xedriver, only : set_pars, func, solver, locerr, callback
      use control, only : ifinit

      implicit none

      character*255 :: param, value

      integer :: dim, res

      real(dp) :: t0, dt
      real(dp) :: X0(2), xe0, io0, xeout
      real(dp) :: yout(2), hold, hnew, tout
!        real(dp) :: sxe(2), ssm(2)

      real(dp) :: lxe, lio, lpm
      real(dp) :: yxe, yio, ypm
      real(dp) :: absxe

      dim = 2
      t0 = 0.D0
      X0(1) = xe0
      X0(2) = io0

      call set_pars(yio,yxe,absxe,lxe,lio)

      ! initialize integrator
      if (ifinit .eq. .False.) then
            call SET_DRIVER(dim, func, solver, locerr, callback)
            ifinit = .True.
      endif

      param = 'LogFileName'
      value = 'out/xe_bdf_out.dat'
      CALL SET_PARAMETER(param, value)

      param = 'InitStepSize'
      write(value,'(F20.10)') 1.E-6
      CALL SET_PARAMETER(param, value)   

      param = 'MaxOrder'
      write(value,'(I10)') 5
      CALL SET_PARAMETER(param, value)

!      CALL SET_PARAMETER('InitStepSize', '1.0')   
!      CALL SET_PARAMETER('Method', METHOD)  
!      CALL SET_PARAMETER('OrderHold','True')
!      CALL SET_PARAMETER('StepSizeHold','True')
!      CALL SET_PARAMETER('MaxStepsNumber', '10000')       
!      CALL SET_PARAMETER('MaxRejectSteps', '100')         
!      CALL SET_PARAMETER('MaxStepSize', '100.0')
!      CALL SET_PARAMETER('MinStepSize', '0.0')

      param = 'EndOfTime'
      write(value,'(F20.10)')dt
      CALL SET_PARAMETER(param, value) 

      CALL SET_INIT_STATE(X0, t0)                         

      CALL CHECK_INPUT_DATA()

      CALL INTEGRATE()

!      CALL GET_PARAMETER('StepsNumber', NH)                  
!      CALL GET_PARAMETER('RejectStepsNumber', NR)            

      CALL GET_RESULT(YOUT, HOLD, HNEW, TOUT, RES)    
 
!      CALL DEALLOC()

      xeout = yout(1)

    end subroutine





    subroutine XENON_ODES(dt,lio,lxe,yio,yxe, &
      frate,arate,xe0,io0,dim,xe1)
      !DEC$ ATTRIBUTES DLLEXPORT :: XENON_ODES
      !
      !
      use precision, only : dp
      use xedriver, only : set_pars, func, solver, locerr, callback
      use control, only : ifinit

      implicit none

      character*255 :: param, value

      integer :: dim, dimtot, res

      real(dp), allocatable :: xinp(:), xout(:)

      real(dp) :: t0, dt
      real(dp) :: hold, hnew, tout
!        real(dp) :: sxe(2), ssm(2)

      real(dp) :: lxe, lio!, lpm
      real(dp), dimension(dim) :: yxe, yio!, ypm
      real(dp), dimension(dim) :: xe0, io0, xe1
      real(dp), dimension(dim) :: arate, frate

      dimtot = 2*dim
      t0 = 0.D0

      ! initialize integrator
      if (ifinit .eq. .False.) then
            allocate(xinp(2*dim), xout(2*dim))
            call driver_setup(yio,yxe,arate,frate,lxe,lio)            
            call SET_DRIVER(dim, func, solver, locerr, callback)
            ifinit = .True.
      endif


      xinp(1     : dim)   = xe0(:)
      xinp(dim+1 : 2*dim) = io0(:)


      param = 'LogFileName'
      value = 'out/xe_bdf_out.dat'
      CALL SET_PARAMETER(param, value)

      param = 'InitStepSize'
      write(value,'(F20.10)') 1.E-6
      CALL SET_PARAMETER(param, value)   

      param = 'MaxOrder'
      write(value,'(I10)') 5
      CALL SET_PARAMETER(param, value)

!      CALL SET_PARAMETER('InitStepSize', '1.0')   
!      CALL SET_PARAMETER('Method', METHOD)  
!      CALL SET_PARAMETER('OrderHold','True')
!      CALL SET_PARAMETER('StepSizeHold','True')
!      CALL SET_PARAMETER('MaxStepsNumber', '10000')       
!      CALL SET_PARAMETER('MaxRejectSteps', '100')         
!      CALL SET_PARAMETER('MaxStepSize', '100.0')
!      CALL SET_PARAMETER('MinStepSize', '0.0')

      param = 'EndOfTime'
      write(value,'(F20.10)')dt
      CALL SET_PARAMETER(param, value) 

      CALL SET_INIT_STATE(X0, t0)                         

      CALL CHECK_INPUT_DATA()

      CALL INTEGRATE()

!      CALL GET_PARAMETER('StepsNumber', NH)                  
!      CALL GET_PARAMETER('RejectStepsNumber', NR)            

      CALL GET_RESULT(YOUT, HOLD, HNEW, TOUT, RES)    
 
!      CALL DEALLOC()

      xe1(:) = yout(1:dim)

    end subroutine

