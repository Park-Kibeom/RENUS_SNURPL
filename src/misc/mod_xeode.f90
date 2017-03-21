module xeode
! Solve ODE equations for Xe-135, I-135 dynamic

    PROCEDURE (), POINTER :: pFunc => NULL()  ! make true NULL-pointer on function

    contains

    subroutine set_xeodefunc(func)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS: "set_xeodefunc" :: set_xeodefunc
    !DEC$ ATTRIBUTES STDCALL :: set_xeodefunc
    !
    ! set external subroutine for xenon ODE's solution data 
        interface 
            subroutine func(deltm, rlambi, rlambxe, gami, gamxe, frate, &
            arate, dfrate, darate, rnxep, rnip, dim, rnxe, rni)

            integer :: dim                              ! dimension of ODE's system (part of xenon only!)
            real(8) :: deltm                            ! time step, s
            real(8) :: rlambi, rlambxe                  ! decay constants for I-135 and Xe-135, 1/s
            real(8), dimension(dim) :: gami, gamxe      ! yield coefficients for I-135 and Xe-135
            real(8), dimension(dim) :: frate, arate     ! rates of fission and xenon's absorbtion, 1/(cm3*s)
            real(8), dimension(dim) :: dfrate, darate   ! derivatives of fission and absorbtion's rates
            real(8), dimension(dim) :: rnxep, rnip      ! initial concentrations for Xe-135 and I-135, 1/cm3
            real(8), dimension(dim) :: rnxe, rni        ! the new concentrations for Xe-135 and I-135, 1/cm3 

            end subroutine
        end interface

        pFunc => func
        write(*,*) 'MESSAGE: Setup of external xenon ODEs solver, loc= ', loc(func)

    end subroutine set_xeodefunc


    subroutine xenon_ODES(  deltm, rlambi, rlambxe, gami, gamxe, frate, &
            arate, dfrate, darate, rnxep, rnip, dim, rnxe, rni)

        integer :: dim                              ! dimension of ODE's system (part of xenon only!)
        real(8) :: deltm                            ! time step, s
        real(8) :: rlambi, rlambxe                  ! decay constants for I-135 and Xe-135, 1/s
        real(8), dimension(dim) :: gami, gamxe      ! yield coefficients for I-135 and Xe-135
        real(8), dimension(dim) :: frate, arate     ! rates of fission and xenon's absorbtion, 1/(cm3*s)
        real(8), dimension(dim) :: dfrate, darate   ! derivatives of fission and absorbtion's rates
        real(8), dimension(dim) :: rnxep, rnip      ! initial concentrations for Xe-135 and I-135, 1/cm3
        real(8), dimension(dim) :: rnxe, rni        ! the new concentrations for Xe-135 and I-135, 1/cm3 

        if (.not.ASSOCIATED(pFunc)) pFunc => DefaultXeODESolver
    
        call pFunc(deltm, rlambi, rlambxe, gami, gamxe, frate, &
            arate, dfrate, darate, rnxep, rnip, dim, rnxe, rni)

    end subroutine xenon_ODES

    subroutine DefaultXeODESolver(  deltm, rlambi, rlambxe, gami, gamxe, frate, &
            arate, dfrate, darate, rnxep, rnip, dim, rnxe, rni)
        integer :: dim                              ! dimension of ODE's system (part of xenon only!)
        real(8) :: deltm                            ! time step, s
        real(8) :: rlambi, rlambxe                  ! decay constants for I-135 and Xe-135, 1/s
        real(8), dimension(dim) :: gami, gamxe      ! yield coefficients for I-135 and Xe-135
        real(8), dimension(dim) :: frate, arate     ! rates of fission and xenon's absorbtion, 1/(cm3*s)
        real(8), dimension(dim) :: dfrate, darate   ! derivatives of fission and absorbtion's rates
        real(8), dimension(dim) :: rnxep, rnip      ! initial concentrations for Xe-135 and I-135, 1/cm3
        real(8), dimension(dim) :: rnxe, rni        ! the new concentrations for Xe-135 and I-135, 1/cm3         
        write(*,*) 'Xenon ODEs solver by default'  
    end subroutine DefaultXeODESolver

end module xeode


