subroutine renusinput(finame)
!DEC$ ATTRIBUTES DLLEXPORT :: renusinput
!DEC$ ATTRIBUTES STDCALL :: renusinput
!DEC$ ATTRIBUTES REFERENCE :: finame
! finame : name of file with input data for renus

    use pydini

    character(255) :: finame

    call timeron() ! init 1st timer 
    call timeron() ! init 2d timer
     
    call default() !initialize default variables

    localfn=trim(finame)
    open(io5,file=localfn, status='old') ! open input file
      
    call initbd()      ! initialize common blocks that are not specified by block data
    call scaninput()   ! read cards from input file

    ! aaaaaa, what is it ???
    open(io8,file=trim(caseid)//'.out',status='unknown')      ! 2014_08_08 . scb
    if(flagout(1))  open(io9,file=trim(caseid)//'.outl',status='unknown')   ! 2014_08_08 . scb
    irstin=io19
    irstout=io20

    call allocpdm0() ! allocate memory required for input processing

    call readinput() ! read input file and cards

end subroutine renusinput


subroutine renusinit
!DEC$ ATTRIBUTES DLLEXPORT :: renusinit

    use pydini

    call init()      ! initialize variables for: xs, geom, TH, BDF, ... 

    call setcntl(ninmax,ifcmfd2g)  ! 2012_11_08 . scb
  
    call setgeom( ng,&
        nx,ny,nz,nxy,&
        ndir,&
        symopt,isymang,                       &
        isymloc,nxs,nxe,nys,nye,nxsf,nxef,                          &
                  kfbeg,kfend,jfbeg,jfend,nodel,                              &
                  hmesh,volnode,volcore,volfuel,                              &
                  (/albxl(1),albxr(1),albyl(1),albyr(1),albzb(1),albzt(1)/),  &
                  ltola,ktoka,iassytyp,iffuela,                               &
                  ifrecsp3,                                                   &   ! 2013_07_05 . scb
                  nxa,nya,nza,rhmesh,ltoi,ltoj)        ! 2013_07_15 . scb

    if(.not.rect) then 
      call setgeomhex( ng,nz,nxy,nassy,nxpnt,ncorn,nxsfc,                    &
                        nsurf,hf2f,ndivhs,kfbeg,kfend,                        &                   
                        volnode,hz,albr(1),albzb(1),albzt(1),                 &
                        wtass,neignd,neigjin,neigpt,neigz,                    &
                        neigsfc,neigsnd,neigsndz,                             &
                        codpnt,pbdv,wtdhat,neigsfcz,ineigcond,                &
                        ipntr,ilubnd,iastopnt,ifhexfdm,ifhexsp3,              &
                        ndiv,neigtri,neigtria)

          if(ifhexfdm) then
            if(ifhexsp3) then
              call malloc_hexfdm3
            else
              call malloc_hexfdm
            endif
          else
            if(ifhexsp3) then
              call malloctpen_sp3
              call mallochex2g_sp3
            endif
            call malloctpen
            call mallochex2g
          endif
    endif

         
       call setxsec(ng,mgb(2),nxy,nz,xstf,xsaf,xsdf,xsnff,      &
                     xskpf,xschif,xsff,xbeta,xsadf,xscdf,xssf,xssfs,   &
                     xssfe,xss2nf,xsdf2,xstrf,xbeta2,xsmax)   
        ! 2013_07_15 . scb :: xbeta2 added
        ! 2013_07_19 . scb :: xsmax added
        ! 2013_10_02 . scb :: xsff added

        call dmalloc(curil,ndirmax,nxy,nz,ng)
        call dmalloc(curir,ndirmax,nxy,nz,ng)
        call dmalloc(curol,ndirmax,nxy,nz,ng)
        call dmalloc(curor,ndirmax,nxy,nz,ng)

        call mallocsfam(ng,nxy,nz,phif,phisfc,jnet,curol,curor,curil,curir)  

  ! 2013_07_16 . scb      
        if(ifrecsp3) then
          call malloccmfd(ng,nxy,nz,nzp1,nsurf)
          call initsp3senm(ng)
        endif
  ! added end      

        if(.not.ifhexfdm) then
          call initsfam(TRUE)
          if(.not.rect) then   
            if(ifhexsp3) then
              call inittpen_sp3
            else
              call inittpen
            endif 
          endif
        endif

      call timeroff(tinit)
end subroutine renusinit


subroutine renustate
!DEC$ ATTRIBUTES DLLEXPORT :: renustate
!DEC$ ATTRIBUTES STDCALL :: renustate
    use pydini
    use dumping,  only : dumpfunc                      !function for dump data into file or using extern function 

! standard code
    call timeron()    ! 2015_06_22 . scb . tsteady
    call runss(ncmfdmgtot,ninmgtot,ncmfd2gtot,nin2gtot)
    call timeroff(tsteady)
    call updxsec(fdbk, FALSE)
                
    eigvsave=eigv 
    eigv0 = eigv   ! 2014_09_04 . scb      

    call dumpfunc() ! save state in file

end subroutine renustate



subroutine renustransient
!DEC$ ATTRIBUTES DLLEXPORT :: renustransient
!DEC$ ATTRIBUTES STDCALL :: renustransient
    use pydini

    if(ifreact) call cal_beta
    call normalize  
    first=FALSE
    if(printff) call writeff
    if(pinpower) call driveppr(FALSE) ! pin power calculation

    call timeron()  

    if(transient) call runtr

    call timeroff(ttransient)

end subroutine renustransient
