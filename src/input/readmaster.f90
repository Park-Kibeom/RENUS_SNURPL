! 2013_10_15 . scb
  subroutine readmaster(indev)
    
    use input
    use allocs
  
    include 'global.h'
    include 'cards.h'
    include 'files.h'
    include 'geom.h'    
  
    integer :: indev
    integer :: iline, irow, jrow
    real(8) :: rdum
    
    logical,save :: firstjobtyp=.true.
    
    do while (.true.)
      read(indev,'(a512)',end=1000) oneline
      if(probe.eq.'#' .or. probe.eq.'/' .or. oneline.eq.BLANK .or.ifnumeric(oneline)) cycle
      
      read(oneline,*) blockname
      call toupper(blockname)
      
      if(blockname.eq.'END') exit
      
      select case(blockname) 
        case('%JOB_TYP') 
          iline=0
          nline=1
          if(firstjobtyp) nline=3
          firstjobtyp=.false.
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.nline) exit
          enddo
		    case('%JOB_TIT')     
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            caseid=oneline
            
            exit
          enddo             
        case('%JOB_VER')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            exit
          enddo
        case('%JOB_IDE')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            exit
          enddo
        case('%JOB_HEX')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            exit
          enddo
        case('%JOB_MDL')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.2) exit
          enddo
        case('%GEN_MTH')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.4) exit
          enddo
        case('%GEN_LMT')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.2) exit
          enddo
        case('%GEN_THD')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.5) exit
          enddo
        case('%GEN_DCH')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.3) exit
          enddo
        case('%GEN_PIN')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            exit
          enddo
        case('%EDT_OUT')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.3) exit
          enddo
        case('%GEN_DIM')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.1) then
              read(oneline,*)  nxa,nya,nza,nbatch,ncomp
            elseif(iline.eq.2) then
              read(oneline,*)  ndim,ngeo,nsym,ndivxy,ndivz
            elseif(iline.eq.3) then
              read(oneline,*)  ng
              exit
            endif            
          enddo
          
          nx=nxa/ndivxy
          ny=nya/ndivxy
          nz=nza/ndivz
          
          nxy=nx*ny
          
          nassytyp=nbatch
                      
          isymang=360/ngeo
        case('%GEN_GEO')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.4) exit
          enddo
        case('%GEN_SYM')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.2) exit
          enddo
        case('%LPD_BCH')          
          allocate(cassymap(nx,ny))
          cassymap=blank
          allocate(cdum(nx))
          cdum=blank
          
          iline=0
          nrow=nx/2
          if(nrow*2 .ne. nx) nrow=nrow+1
          irow=nrow-1
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            irow=irow+1
            read(oneline,*) cassymap(1:irow,iline)
            
            if(iline.eq.nrow) exit
          enddo
          
          jrow=0  ! for saving dummy assy type
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            irow=irow-1
            jrow=jrow+1
            read(oneline,*) cdum(jrow),cassymap(1:irow,iline)
            
            if(iline.eq.ny) exit
          enddo
          
        case('%LPD_B&C')
          allocate(cassyname(nbatch))
          cassyname=blank
          call dmalloc(iassycomp,nz,nbatch)
          
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            read(oneline,*) cassyname(iline),iassycomp(1:nz,iline)
            
            if(iline.eq.nbatch) exit
          enddo
        case('%LPD_C&X')
          iline=0
          nfcomp=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            read(oneline,*) icompnum(iline), nmxset(iline), indblock(iline)
            
            if(icompnum(iline).gt.0) nfcomp=nfcomp+1
            
            if(iline.eq.ncomp) exit
          enddo
        case('%LPD_HFF')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            read(oneline,*)
            
            if(iline.eq.nfcomp) exit
          enddo
        case('%EDT_PIN')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.2) exit
          enddo
        case('%MSC_MEM')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            exit
          enddo
        case('%COB_INP')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.3) exit
          enddo
        case('%ROD_CFG')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            read(oneline,*) ncrgr              
            
            exit
          enddo
          
          allocate(crname(ncrgr))
          crname=blank
          call dmalloc(mattip,ncrgr)
          call dmalloc(matabs,ncrgr)
          call dmalloc(matfol,ncrgr)
          call dmalloc(ifgtp,ncrgr)
          call dmalloc(lentip,ncrgr)
          call dmalloc(lenabs,ncrgr)
          call dmalloc(crups,ncrgr)
          call dmalloc(crlos,ncrgr)
          call dmalloc(crpos,ncrgr)
          call dmalloc(crstp,ncrgr)
          
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            read(oneline,*) crname(iline),mattip(iline),matabs(iline),matfol(iline),lentip(iline),lenabs(iline), &
                            crups(iline),crlos(iline),crpos(iline),crstp(iline),ifgtp(iline)
            
            if(iline.eq.ncrgr) exit
          enddo          
        case('%ROD_COR')
          
          
          
          
          
        case('%EXE_STD')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            exit
          enddo
        case('%EXE_DEP')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.2) exit
          enddo
        case('%EDT_OPT')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            exit
          enddo
        case('%EXE_ROD')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.ncrgr) exit
          enddo
        case('%TRN_OPT')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.4) exit
          enddo
        case('%TRN_EPS')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.2) exit
          enddo
        case('%TRN_DEF')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            read(oneline,*) nprec
            
            exit
          enddo
          
          call dmalloc(chid,nprec,ng)
          call dmalloc(betij,nprec)
          call dmalloc(lambd,nprec)
          call dmalloc(rvel,ng)
          
          
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            read(oneline,*) chid(iline,1:ng)
            
            if(iline.eq.nprec) exit
          enddo
          
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.1) then
              read(oneline,*) betij(1:nprec)
            elseif(iline.eq.2) then
              read(oneline,*) lambd(1:nprec)
            elseif(iline.eq.3) then
              read(oneline,*) rvel(1:ng)
              exit
            endif            
          enddo
        case('%TRN_PRN')
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle

            read(oneline,*) rdum
            if(rdum.lt.0) exit
          enddo
        case('%TRN_SHD')
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            if(iline.eq.3) then
              read(oneline,*) nstuck
              exit
            endif            
          enddo
          allocate(crstuck(nstuck))
          nstuck=blank
          
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            read(oneline,*) crstuck(iline)
            
            if(iline.eq.nstuck) exit
          enddo          
        case('%EXE_TRN')     
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            read(oneline,*) nconf
            
            exit
          enddo
          call dmalloc(timechg,nconf)
          call dmalloc(ppm,nconf)
          call dmalloc(irodchg,nconf)
          allocate(crtype(ncrgr,nconf))
          crtype=blank
          call dmalloc(crstep,ncrgr,nconf)
          
          iline=0
          do while(.true.)
            read(indev,'(a512)',end=1000) oneline
            if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
            iline=iline+1
            
            read(oneline,*) timechg(iline),ppm(iline),irodchg(iline)
            
            if(irodchg(iline).gt.0) then
              jline=0
              do while(.true.)
                read(indev,'(a512)',end=1000) oneline
                if(probe.eq.'#' .or. oneline.eq.BLANK) cycle
            
                jline=jline+1
            
                read(oneline,*) crtype(jline,iline),crstep(jline,iline)
            
                if(jline.eq.irodchg(iline)) exit
              enddo
            endif                          
            
            if(iline.eq.nconf) exit
          enddo
        case default
          print *, blockname
          go to 200
        end select 
        
 200    continue
      enddo                     
         
 1000 continue
         
  end subroutine
! added end  