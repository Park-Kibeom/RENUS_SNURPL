module dumping
! Dump data from RENUS using external functions

    interface 
        subroutine DummyDumpFunc
        end subroutine
    end interface

    PROCEDURE (), POINTER :: pDumpFunc => NULL()  ! make true NULL-pointer 
!    POINTER(pDumpFunc, DummyDumpFunc)            ! another way to handle with pointers

    integer :: i32array(100)                         ! array for save int32 variables
    real    :: f64array(100)                         ! array for save float64 variables

    contains


    subroutine set_dumpfunc(func)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS: "set_dumpfunc" :: set_dumpfunc
    !DEC$ ATTRIBUTES STDCALL :: set_dumpfunc
    !
    ! set extern subroutine for dumping data 
        interface 
            subroutine func
            end subroutine
        end interface

        pDumpFunc => func
        write(*,*) 'MESSAGE: Set external dump function, loc= ', loc(func)

    end subroutine set_dumpfunc


    subroutine join2data(i32ptrs, i32lens, i64ptrs, i64lens, &
                         f32ptrs, f32lens, f64ptrs, f64lens)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS: "join2data" :: join2data
    !DEC$ ATTRIBUTES STDCALL :: join2data
    !DEC$ ATTRIBUTES REFERENCE :: i32ptrs, i32lens, i64ptrs, i64lens
    !DEC$ ATTRIBUTES REFERENCE :: f32ptrs, f32lens, f64ptrs, f64lens
    !
    ! pass pointers on necessary arrays
        use sfam, only : eigv   !variables: multiplication factor

        include 'global.h'      !variables: nx, ny, nxy, ng, ...
        include 'ffdm.h'        !variables: phi, psi, phif, psif, ...
        include 'thgeom.inc'    !variables: nchan, nzth, nr, nr1, ..., nr5 - number of TH-channels and axial's th-nodes + number of nodes in the pellet region
        include 'thfdbk.inc'    !variables: tcool, dcool for coolant temperatures and density
        include 'thfuel.inc'    !variables: tfuel for fuel temperatures         
        include 'pow.h'         !variables: plevel, absp, relp - power level, absolute and relative power
        include 'trancntl.inc'  !variables: time
        include 'srchppm.h'     !variables: boron concentration, ppm
        include 'xesm.h'        !variables: rnxe, rni, rnsm, rnpm for xenon/samarium concentrations

        integer, parameter :: lenptrs=100
        integer :: i32ptrs(lenptrs), i32lens(lenptrs)
        integer :: i64ptrs(lenptrs), i64lens(lenptrs)
        integer :: f32ptrs(lenptrs), f32lens(lenptrs)
        integer :: f64ptrs(lenptrs), f64lens(lenptrs)

        i32array(1) = ng
        i32array(2) = nx
        i32array(3) = ny
        i32array(4) = nxy
        i32array(5) = nz
        i32array(6) = nchan
        i32array(7) = nzth
        i32array(8) = nr

        f64array(1) = 0.
        f64array(2) = 0.

        ! variables int32
        i32ptrs(1) =  loc(i32array)
        i32lens(1) = size(i32array)

        ! variables int64

        ! variables float32

        ! variables float64
        ! array with variables
        f64ptrs(1) =  loc(f64array)
        f64lens(1) = size(f64array)

        !time, sec
        f64ptrs(2) =  loc(time)
        f64lens(2) = 1

        !level of power, rel.units
        f64ptrs(3) =  loc(plevel)
        f64lens(3) = 1       

        !neutron flux
        f64ptrs(4) =  loc(phif)
        f64lens(4) = size(phif)

        !reletive power, Wt/cm3
        f64ptrs(5) =  loc(relp)
        f64lens(5) = size(relp)

        !absolute power, Wt
        f64ptrs(6) =  loc(absp)
        f64lens(6) = size(absp)

        !coolant temperature, C
        f64ptrs(7) =  loc(tcool)
        f64lens(7) = size(tcool)

        !coolant density, g/cm3
        f64ptrs(8) =  loc(dcool)
        f64lens(8) = size(dcool)

        !fuel temperature, C
        f64ptrs(9) =  loc(tfuel)
        f64lens(9) = size(tfuel)

        !boron concentration, ppm
        f64ptrs(10) =  loc(ppm)
        f64lens(10) =  1

        !xenon-135 concentration, 1/cm3
        f64ptrs(11) =  loc(rnxe)
        f64lens(11) =  size(rnxe)

        !iodine-135 concentration, 1/cm3
        f64ptrs(12) =  loc(rni)
        f64lens(12) =  size(rni)

        !prometium concentration, 1/cm3
        f64ptrs(13) =  loc(rnpm)
        f64lens(13) =  size(rnpm)

        !samarium concentration, 1/cm3
        f64ptrs(14) =  loc(rnsm)
        f64lens(14) =  size(rnsm)

        !multiplication factor
        f64ptrs(15) =  loc(eigv)
        f64lens(15) =  1

        !Norm coefficient for power, MWt
!        f64ptrs(10)=  loc(fnorm)
!        f64lens(10)= 1  


    end subroutine join2data


    subroutine DefaultDumpFunc
        write(*,*) 'dump default'  
    end subroutine DefaultDumpFunc

    subroutine dumpfunc() 
        if (.not.ASSOCIATED(pDumpFunc)) pDumpFunc => DefaultDumpFunc
        call pDumpFunc()    
    end subroutine dumpfunc


end module dumping


