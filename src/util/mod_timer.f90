    module timer
      use omp_lib
      integer, parameter :: size=100
      real, dimension(:) :: starttime(size)=0
      integer :: index = 0
      
      contains
        subroutine timeron()
          integer :: timevals(8)
          if(index.ge.size) return
          index = index+1
! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
!          call date_and_time(values = timevals)
!          starttime(index)=chglong(timevals)
          !call cpu_time(starttime(index))
          starttime(index) =  OMP_GET_WTIME() 
! added end
        end subroutine
         
        subroutine timeroff(totalelapsed)
          integer :: timevals(8)
          real ::  elapsed, totalelapsed
          if(index.le.0) return
! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
!          call date_and_time(values=timevals)
!          elapsed=chglong(timevals)
          !call cpu_time(elapsed)
          elapsed =  OMP_GET_WTIME() 
          !elapsed = 1.d0
! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
          totalelapsed = totalelapsed+elapsed - starttime(index)
          index = index-1 
        end subroutine
         
        subroutine reset()
          starttime=0
          index=0
        end subroutine

        function chglong(timevals)         
          integer :: timevals(8)
          real :: chglong
          chglong=timevals(5)*3600.0+timevals(6)*60.0+timevals(7)+timevals(8)*0.001
        end function
         
    end module
