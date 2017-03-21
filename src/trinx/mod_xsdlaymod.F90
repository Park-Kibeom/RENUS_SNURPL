! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB  
        module XSDlayMod
           use param
           use BasicParam
           use ErrorHandle
           use allocs

           implicit none
           

           integer, parameter :: ADLA_BAS = 1                                            ! A.DLA original format
           integer, parameter :: ADLA_EXT = 2                                            ! A.DLA extend format by hanty

           type XSDlay
              

              integer(NBI) :: no_group      
              integer(NBI) :: no_isotope
              integer(NBI) :: no_family                                           ! number of delayed neutron precursor group
              integer(NBI) :: no_chid                                             ! number of delayed emission spectrum

              integer(NBI), pointer, dimension(:) :: inkfam                       ! ??
              integer(NBI), pointer, dimension(:) :: skip_block                   ! number of records to be skipped to read data for isotope i. 
              integer(NBI), pointer, dimension(:,:) :: no_infam                         ! ??
               

         
              character(ISTN), pointer, dimension(:) :: iso_name

              real(NBF), pointer, dimension(:,:) :: chid                          ! delayed neutron emission spectrum
              real(NBF), pointer, dimension(:,:,:) :: chid_ext                    ! delayed neutron emission spectrum per isotope
              real(NBF), pointer, dimension(:) :: xs_dct                          ! decay constant
              real(NBF), pointer, dimension(:,:) :: xs_dct_ext                    ! decay constant per isotope 
              real(NBF), pointer, dimension(:,:,:) :: dnu                         ! delayed neutron yield
              real(NBF), pointer, dimension(:) :: max_energy                      ! energy maximum bound 
              real(NBF) :: min_energy

           end type XSDlay
           
           !type(XSDlay) :: dlayxs



           contains

           subroutine init_xsdlay( dlayxs )
              type(XSDlay), intent(inout) :: dlayxs
              integer(NBI) :: i, m
              
              character(ISTN), parameter :: BLANK = '        '

              allocate( dlayxs%iso_name( dlayxs%no_isotope ) )
              do i=1, dlayxs%no_isotope
                 dlayxs%iso_name(i) = '         '
              enddo
              call dmalloc( dlayxs%xs_dct,   dlayxs%no_family  )
              call dmalloc( dlayxs%xs_dct_ext, dlayxs%no_family, dlayxs%no_isotope )
              call dmalloc( dlayxs%inkfam,   dlayxs%no_isotope )
              call dmalloc( dlayxs%skip_block, dlayxs%no_isotope )

              dlayxs%no_chid = dlayxs%no_group*dlayxs%no_family
              call dmalloc( dlayxs%chid, dlayxs%no_group, dlayxs%no_family )
              call dmalloc( dlayxs%chid_ext, dlayxs%no_group, dlayxs%no_family, dlayxs%no_isotope )
              call dmalloc( dlayxs%max_energy, dlayxs%no_group )

           end subroutine init_xsdlay

           subroutine read_xsdlay( dlayxs, dly_flag, filedla )
              type(XSDlay), intent(inout) :: dlayxs
              integer(NBI) :: dly_flag
					character*72 :: filedla

              if( dly_flag .eq. ADLA_BAS )  then
	             call read_xsdlay_bas( dlayxs, filedla )
              else if( dly_flag .eq. ADLA_EXT ) then
                 call read_xsdlay_ext( dlayxs, filedla )
                 call trans_xsdlay_form( dlayxs )
!				print *, sum(dlayxs%dnu(:,1,1))
              else 
                 call write_error3( 'DLAYXS Flag Error', dly_flag )
              endif        

           end subroutine read_xsdlay

           subroutine read_xsdlay_bas( dlayxs, filedla )

			  character*72 :: filedla
              type(XSDlay), intent(inout) :: dlayxs

              logical        :: check = .TRUE.
              character(MCL) :: hname, bname
              character(MCL) :: ctmp1, ctmp2, ctmp3, ctmp4
              integer(NBI)   :: itmp, maxn
              integer(NBI)   :: niso, ng, nfam, max_fam, ninfam, nsnudel, max_del_size

              integer(NBI)   :: g, i, j, k, l
              real(NBF)      :: rtmp
           
        !  file open
              open(nadla,file=filedla)

        ! 0V
        !      read(nadla,*) hname, maxn, ctmp2, itmp                                                                    
        !      read(nadla,*) (itmp, i=1,maxn)
              read(nadla,*) bname                                                                                

              call check_word( bname, '0V', check )
              if( check .eq. .FALSE. )  call write_error( bname )

        ! 1D    
              read(nadla,*) bname, dlayxs%no_group, dlayxs%no_isotope, dlayxs%no_family         

              call check_word( bname, '1D', check )
              if( check .eq. .FALSE. )  call write_error( bname )

              call init_xsdlay(dlayxs)
              niso = dlayxs%no_isotope
              nfam = dlayxs%no_family
              ng   = dlayxs%no_group

        !  2D    
!              read(nadla,'(1X,A2,1X,9(1X,A6)/(10(1X,A6)))') bname, (dlayxs%iso_name(i),i=1,niso)
			read(nadla,*) bname, (dlayxs%iso_name(i),i=1,niso)
              call check_word( bname, '2D', check )
              if( check .eq. .FALSE. )  call write_error( bname )

              read(nadla,*)  (dlayxs%xs_dct(j),j=1,nfam),   ((dlayxs%chid(g,j),g=1,ng),j=1,nfam),       &
     	                     (dlayxs%max_energy(g),g=1,ng), dlayxs%min_energy 
                                                                  
              read(nadla,*)  (dlayxs%inkfam(i),i=1,niso), (dlayxs%skip_block(i),i=1,niso)

        !  get max nkfam
              max_fam = 0
              do i=1, niso
                 if( max_fam .lt. dlayxs%inkfam(i) ) max_fam =  dlayxs%inkfam(i)
              enddo
        ! 
              max_del_size = max_fam*ng+1
        !!
              call dmalloc(dlayxs%no_infam,max_fam,niso)
              call dmalloc(dlayxs%dnu,ng,max_fam,niso)
        !!
        !  3D

              do i=1,niso
                 ninfam  = dlayxs%inkfam(i)
                 read(nadla,*) bname, ((dlayxs%dnu(g,j,i),g=1,ng),j=1,ninfam) 
                 call check_word( bname, '3D', check )
                 if( check .eq. .FALSE. )  call write_error( bname )

                 read(nadla,*) (dlayxs%no_infam(j,i),j=1,ninfam)
              enddo

           end subroutine read_xsdlay_bas


           subroutine read_xsdlay_ext( dlayxs, filedla )

			  character*72 :: filedla
              type(XSDlay), intent(inout) :: dlayxs

              logical        :: check = .TRUE.
              character(MCL) :: hname, bname
              character(MCL) :: ctmp1, ctmp2, ctmp3, ctmp4
              integer(NBI)   :: itmp, maxn
              integer(NBI)   :: niso, ng, nfam, max_fam, ninfam, nsnudel, max_del_size

              integer(NBI)   :: g, i, j, k, l
              real(NBF)      :: rtmp
           
        !  file open
              open(nadla,file=filedla)

        ! 0V
        !      read(nadla,*) hname, maxn, ctmp2, itmp                                                                    
        !      read(nadla,*) (itmp, i=1,maxn)
              read(nadla,*) bname                                                                                

              call check_word( bname, '0V', check )
              if( check .eq. .FALSE. )  call write_error( bname )

        ! 1D    
              read(nadla,*) bname, dlayxs%no_group, dlayxs%no_isotope, dlayxs%no_family         

              call check_word( bname, '1D', check )
              if( check .eq. .FALSE. )  call write_error( bname )

              call init_xsdlay(dlayxs)
              niso = dlayxs%no_isotope
              nfam = dlayxs%no_family
              ng   = dlayxs%no_group

        !  2D    
              read(nadla,'(1X,A2,1X,9(1X,A6)/(10(1X,A6)))') bname, (dlayxs%iso_name(i),i=1,niso)
              call check_word( bname, '2D', check )
              if( check .eq. .FALSE. )  call write_error( bname )

              read(nadla,*)  ((dlayxs%xs_dct_ext(j,i),j=1,nfam),i=1,niso),                 &                   
                            (((dlayxs%chid_ext(g,j,i),g=1,ng),j=1,nfam),i=1,niso),         &                   
                            (dlayxs%max_energy(g),g=1,ng), dlayxs%min_energy 
                                                                  
              read(nadla,*)  (dlayxs%inkfam(i),i=1,niso), (dlayxs%skip_block(i),i=1,niso)

        !  get max nkfam
              max_fam = 0
              do i=1, niso
                 if( max_fam .lt. dlayxs%inkfam(i) ) max_fam =  dlayxs%inkfam(i)
              enddo
        ! 
              max_del_size = max_fam*ng+1
        !!
              call dmalloc(dlayxs%no_infam,max_fam,niso)
              call dmalloc(dlayxs%dnu,ng,max_fam,niso)
        !!
        !  3D

              do i=1,niso
                 ninfam  = dlayxs%inkfam(i)
                 read(nadla,*) bname, ((dlayxs%dnu(g,j,i),g=1,ng),j=1,ninfam) 
                 call check_word( bname, '3D', check )
                 if( check .eq. .FALSE. )  call write_error( bname )

                 read(nadla,*) (dlayxs%no_infam(j,i),j=1,ninfam)
              enddo

           end subroutine read_xsdlay_ext


           subroutine trans_xsdlay_form( dlayxs )
              type(XSDlay), intent(inout) :: dlayxs


           end subroutine trans_xsdlay_form

        end module XSDlayMod