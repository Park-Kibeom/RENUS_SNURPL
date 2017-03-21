    module control
      use precision, only : dp

      logical :: ifinit = .False.

    
      contains

    end module control


    subroutine XENON_ODES(dt,lio,lxe,yio,yxe, &
      frate,arate,dfrate,darate,xe0,io0,dim,xe1,io1)
      !DEC$ ATTRIBUTES DLLEXPORT,ALIAS: "xenon_odes" :: xenon_odes
      !DEC$ ATTRIBUTES STDCALL :: xenon_odes
      !DEC$ ATTRIBUTES REFERENCE :: dt, lio, lxe, dim
      !
      !
      use precision, only : dp
      use xedriver, only : driver_setup, func, solver, &
                               locerr, callback, rho_io
      use control, only : ifinit !, poisdat

      implicit none

      character*255 :: param, value

      integer :: dim, res

      real(dp), allocatable :: xinp(:), xout(:)

      real(dp) :: t0, dt
      real(dp) :: hold, hnew, tout
!        real(dp) :: sxe(2), ssm(2)

      real(dp) :: lxe, lio!, lpm
      real(dp), dimension(dim) :: yxe, yio!, ypm
      real(dp), dimension(dim) :: xe0, io0, xe1, io1
      real(dp), dimension(dim) :: arate, frate, darate, dfrate

      t0 = 0.D0

      write(*,*) 'dt = ', dt
      write(*,*) 'lio = ', lio
      write(*,*) 'lxe = ', lxe
      write(*,*) 'yxe = ', yxe(1:4)
      write(*,*) 'yio = ', yio(1:4)
!      write(*,*) 'xe0 = ', xe0(1:4)
!      write(*,*) 'dim = ', dim
      stop

      ! initialize driver and integrator
      if (ifinit .eq. .False.) then
            call driver_setup(lxe,lio,dim,yxe,yio,&
                  arate,frate,darate,dfrate,io0)
            call SET_DRIVER(dim, func, solver, locerr, callback)
            ifinit = .True.
      endif


      allocate(xinp(dim), xout(dim))

      xinp(:) = xe0(:)

      param = 'LogFileName'
      value = 'out/xe_bdf_out.dat'
      CALL SET_PARAMETER(param, value)

      param = 'InitStepSize'
      write(value,'(F20.10)') 1.
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

      param = 'MaxStepSize'
      write(value,'(F20.10)') 1000.
      CALL SET_PARAMETER(param, value)
!      CALL SET_PARAMETER('MinStepSize', '0.0')

      param = 'EndOfTime'
      write(value,'(F20.10)')dt
      CALL SET_PARAMETER(param, value) 

      CALL SET_INIT_STATE(XINP, t0)                         

      CALL CHECK_INPUT_DATA()

      CALL INTEGRATE()

!      CALL GET_PARAMETER('StepsNumber', NH)                  
!      CALL GET_PARAMETER('RejectStepsNumber', NR)            

      CALL GET_RESULT(XOUT, HOLD, HNEW, TOUT, RES)    
 
!      CALL DEALLOC()

      xe1(:) = xout(:)
      io1(:) = rho_io(dt)

      deallocate(xinp, xout)

    end subroutine

