    subroutine qtrans(lendt, dtlist)
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : '_QTRANS' :: qtrans        
    !DEC$ ATTRIBUTES STDCALL :: qtrans
    !DEC$ ATTRIBUTES REFERENCE :: dtlist
    !
    ! This subroutine is designed for integrate of neutron, TH and poisons equations
    ! in quasi-state approach. On each fixed time step given from input file the 
    ! poison equations are integrated by analytically or BDF scheme. Then steady-state
    ! calculations for neutron and TH models are performed. 
    ! Scheme of the subroutine:
    !
    ! - Solve eigenvalue problem
    ! Do:
    !   - add time step
    !   - move of control rods
    !   Do:
    !      - integration of xenon equations
    !      Do:
    !          - for given 'eigv' get the boron conc. or power level 
    !          - steady-state TH simulation
    !          - update cross sections and mcfd coefficients
    !          - Solve eigenvalue problem => new 'eigv'
    !
    ! Called functions and subroutines:
    ! resetsfam
    ! runsteady
    ! runsteady_hex
    ! perturb
    ! misc_rate
    ! xenon_ODES
    ! updrelpow
    ! drivethss, drivethtr
    ! misc_fluxlevel
    ! updxsec
    ! updppm, updmwt



        use param
!        use timer 
        use sfam,     only : runsteady,resetsfam, eigv, reigv, &
                           runsteady_hex,runsteady_hex_nodal   
        use sfam,     only : psi 
        use hexfdm,   only : drive_hexfdm
        use hexfdm3,  only : drive_hexfdm3
        use sfam_cntl,only : ifcmfd
        use itrinfo,  only : iterinner2g, iterinnermg 
        use tran,     only : deltm, deltmd
        use Mod_FixedSource, only : resid,iffixed    
        use misc,     only : misc_fluxlevel, misc_rate            ! function for rates and flux level calculation    
        use dumping,  only : dumpfunc                      !function for dump data into file or using extern function 
        use xeode,    only : xenon_ODES                    !solver for xenon differencial equations

        include 'global.h'
        include 'itrcntl.h' 
        include 'ffdm.h'           !variables: phi, psi, phif, psif, ...
        include 'files.h'
        include 'geom.h'
        include 'xsec.h'
        include 'nodal.h'
        include 'times.h'
        include 'srchppm.h'        !boron concentration, ppm
        include 'thfdbk.inc'       !coolant temperatures
        include 'thfuel.inc'       !fuel temperatures 
        include 'thcntl.inc'
        include 'pow.h'
        include 'thlink.inc'       ! FREK/MARS
        include 'hexfdm.h'
        include 'geomh.h'
        include 'cntl.h'           ! 2013_05_16 . scb for restart
        include 'xesm.h'           ! xenon, samarium, iodime and prometium concentrations: rnxe, rni, rnsm, rnpm
        include 'trancntl.inc'        

!        logical :: iftrans = TRUE

        ! variables used for iteration scheme between poison and TH/neutron models:
        integer :: count              !counter for outer iterations
        integer :: count_max = 100    !max value of iterations 
!        integer :: noutmax = 100      !number of iteration for boron/power search

        real(8) :: error_max = 1.D-7  ! max variance in neutron flux between two last iterations
        real(8) :: error=1.D0         ! current values of error


!        integer :: nout

        ! variables for steady-state calculations:
        logical :: ssconvchk
        logical :: ifupdxsc=FALSE
        logical :: ifupdls=TRUE
        logical :: flag2g
        logical, save :: first=TRUE ! this flag leads to high number (ncmfd=10000) of eigenvalue iterations

        integer :: noutbegmg=0, ninbegmg=0
        integer :: noutbeg2g=0, ninbeg2g=0
        integer :: dim      
        integer :: num_dt=0 !current number of time step

        integer :: lendt         ! number of time steps
        real(8), optional :: dtlist(lendt) ! list of time steps
        real    :: eigvfirst
        real    :: eigvdif
        real    :: eigvnodal
        common / convcheck / eigvfirst
        common / convcheck2 / eigvnodal

!        write(*,*) 'lendt: ',  lendt, PRESENT(dtlist)
!        write(*,*) 'dtlist: ', loc(dtlist)
!        write(*,*) 'dtlist: ', dtlist(1:lendt)
     

        if (PRESENT(dtlist)) tend = sum(dtlist)
        time = 0.D0 ! start to integrate from this initial time

        deltm  = deltm0 !current  time step size 
        deltmd = deltm0 !previous time step size

        flag2g = ng2.eq.ng ! check on two groups

        call resetsfam() !update nodal and mcfd coefficients
!        iftran = FALSE

!        call exfbdata(TRUE,nthdim,fbdata,bankpos,iok)     
        tmmax = maxval(tcool(:,:))
        tcoolavg = 0.
        tfmax = 0.
        dcoolavg = 0.

        ! iterations by time
        do while (TRUE)

            ! check the time step on excess of total time
            num_dt = num_dt + 1
            if (PRESENT(dtlist)) deltm = dtlist(num_dt)
            if (time+deltm .gt. tend) deltm = tend - time
            time = time + deltm  ! increasing of time value by one time step

            ! moving of control rods
            call perturb(deltm)
            
            ! simple iterations between xenon and steady-state equations
            phifp_(:,:,:) = 0.
            phifp(:,:,:)  = phif(:,:,:)
            tcoolp(:,:)   = tcool(:,:)
            tfuelp(:,:,:) = tfuel(:,:,:)

            do count=0, count_max

                !transient calculations of xenon equations             
                if (flagxesm) then
                    call misc_rate(phif,  xsff,   frate)  ! fission rate
                    call misc_rate(phif,  xsxeaf, arate)  ! rate of absortion by xenon
                    call misc_rate(phifp, xsff,   dfrate) ! on previous time step
                    call misc_rate(phifp, xsxeaf, darate) ! on previous time step

                    ! normolize of arrays
                    dfrate(:,:) = fnorm/deltm * (frate - dfrate)
                    darate(:,:) = fnorm/deltm * (arate - darate)
                    frate(:,:)  = fnorm * frate(:,:)
                    arate(:,:)  = fnorm * arate(:,:)
                    dim = size(frate)

                    ! integrate of xenon and iodine equations
                    call xenon_ODES(     &
                        deltm,           & ! time step, s
                        rlambi, rlambxe, & ! decay constants for I-135 and Xe-135, 1/s
                        gami, gamxe,     & ! yield coefficients for I-135 and Xe-135
                        frate, arate,    & ! rates of fission and xenon's absorbtion, 1/(cm3*s)
                        dfrate, darate,  & ! derivatives of fission and absorbtion's rates
                        rnxep, rnip,     & ! initial concentrations for Xe-135 and I-135, 1/cm3
                        dim,             & ! dimension of ODE's system (part of xenon only!)
                        rnxe, rni)         ! the new concentrations for Xe-135 and I-135, 1/cm3 
                endif


            ! iterations for total power or boron concentration
            first = TRUE
            do nout=0, noutmax

                srchmwt = .true.  
                srchppm = .false.
                if(srchppm .and. nout.ne.1) call updppm(eigv,1.0) ! search boron concentration, ppm  
                if(srchmwt) call updmwt(eigv,1.0,nout) ! search power level, rel.values  

                call updrelpow(FALSE)        !update power (ifupdplevel=FALSE)
                
                ! TH simulation
                if (fdbk) then
                    tcool(:,:)   = tcoolp(:,:)
                    tfuel(:,:,:) = tfuelp(:,:,:)
                    thetaf = 1.              ! setup full implicit euler scheme for fuel and coolant integration
                    thetac = 1.
                    call drivethss           !steady-state TH simulation
!                    call drivethtr(deltm)        !transient TH simulation
                endif

                call updrodfrac()            !update rod frations in nodes
                call updxsec(fdbk, FALSE)    !update cross sections (iftrans = FALSE)
                call resetsfam()             !update nodal and mcfd coefficients

                !solve of eigenvalue problem 
                if(rect) then
                    call runsteady(flag2g,noutbegmg,ninbegmg,&
                        noutbeg2g,ninbeg2g,epsl2,erreig,errl2)
                else
                    call runsteady_hex(flag2g,noutbegmg,ninbegmg,&
                        noutbeg2g,ninbeg2g,epsl2,erreig,errl2)
                endif

                !write(*,*) time, count, nout

                !check convergence for steady-state (eigenvalue) iterations
                if(ssconvchk(eigv,errl2)) exit

            enddo ! end of eigenvalue iterations

                iterinner2g=noutbeg2g
                iterinnermg=noutbegmg

                fnorm = misc_fluxlevel(phif)
                phifp_(:,:,:) = phifp(:,:,:) 
                phifp(:,:,:)  = phif(:,:,:)

                error = sum( (phifp - phifp_)**2 )**0.5 / sum(phifp)

                if ((count>2).and.(error<error_max)) exit !when difference between fluxes on two last iterations will be less than error_max
                if(count .eq. count_max-1) then 
                    write(*,*) 'ERROR: Xenon dynamic did not converge, &
                    error=', error, '> ', error_max
                    stop
                endif

            enddo ! end of iterations between neutrons and xenon equations

        if(flagxesm) call updxesm()              !save values of xenon and iodine concentrations
        call editout( .true. , time )   ! 2014_10_06 . scb
        tmmax = maxval(tcool(:,:))
        call dumpfunc()                  !save results in file

!        call updrelpow(TRUE)
!        call updxsec(fdbk, TRUE)    !update cross sections (iftrans = TRUE)

       
        if (time .ge. tend) exit
!        write(*,*) 'Time, s: ', time, &
!        'Power, %: ', 100*plevel,'b.c., ppm: ', ppm

        enddo ! end iterations by time steps

    end subroutine
