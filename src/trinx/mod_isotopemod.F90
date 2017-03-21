! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB       
        module IsotopeMod
           use param
           use BasicParam
           use XSConstMod
           use ErrorHandle
           use allocs

           implicit none
           integer, parameter :: TOTAL_SCAT = 0, ELASTIC_SCAT = 100, INELASTIC_SCAT = 200, N2N_SCAT = 300 

!#define MOD_ISO
           type Isotope
             
              character(ISTN) :: name
              real(NBF) :: mass 
	          real(NBF) :: kappa                                ! total energy yield per fission                 
	          real(NBF) :: energy_yld_capt                      ! total energy yield per capture
	          real(NBF) :: temperature               
	          real(NBF) :: potent_scat_res                      ! potential scattering in resonance range  [barn/atom]   
	          real(NBF) :: density
        	  
	          real(NBF), pointer, dimension(:,:) :: pl_xs_tr    ! PL weighted transport xs
	          real(NBF), pointer, dimension(:,:) :: pl_xs_tt    ! PL weighted total xs
              real(NBF), pointer, dimension(:,:) :: dir_xs_tr   ! coordinate direction i transport xs

              real(NBF), pointer, dimension(:) :: xs_ngm        ! (n,gamma)  xs
	          real(NBF), pointer, dimension(:) :: xs_nf         ! (n,f)      xs
              real(NBF), pointer, dimension(:) :: xs_na         ! (n,alpha)  xs 
              real(NBF), pointer, dimension(:) :: xs_np         ! (n,p)      xs
	          real(NBF), pointer, dimension(:) :: xs_n2n        ! (n,n2n)    xs
              real(NBF), pointer, dimension(:) :: xs_nd         ! (n,d)      xs 
              real(NBF), pointer, dimension(:) :: xs_nt         ! (n,t)      xs
	          real(NBF), pointer, dimension(:) :: tnu           ! total neutron yield/fission 
	          real(NBF), pointer, dimension(:) :: chiv          ! chi(g) vector : flag_chi = 1

              real(NBF), pointer, dimension(:,:) :: xs_scat
              real(NBF), pointer, dimension(:,:) :: xs_n2ns     ! n2n scattering xs
              real(NBF), pointer, dimension(:,:) :: scatmat     

	          real(NBF), pointer, dimension(:,:) :: chim         ! chi(k,g) matrix : flag_chi >1 
	                                                             ! fraction of neutrons emitted into group g from any group, using spectrum k

              integer(NBI) :: no_group
              integer(NBI) :: class                              ! 0 : undefined,       1 : fissile,    2 : fertile,  3 : other actinide, 
	                                                             ! 4 : fission product, 5 : structure,  6 : coolant,  7 : control
	          integer(NBI) :: flag_chi                           ! 0 : use file-wide chi,   1 : isotope chi vector,  >1 : isotope chi matrix
              integer(NBI) :: flag_xs_nf                         ! (n,f)     XS flag ,      0 : no data,    1 : present
	          integer(NBI) :: flag_xs_na                         ! (n,alpha) XS flag
	          integer(NBI) :: flag_xs_np                         ! (n,p)     XS flag
	          integer(NBI) :: flag_xs_n2n                        ! (n,2n)    XS flag
	          integer(NBI) :: flag_xs_nd                         ! (n,d)     XS flag
	          integer(NBI) :: flag_xs_nt                         ! (n,t)     XS flag
	          integer(NBI) :: no_mnt_xs_tt                       ! number of moments of total xs provided in xs record
	          integer(NBI) :: no_mnt_xs_tr                       ! number of moments of transport xs provided in xs record
	          integer(NBI) :: no_dir_xs_tr                       ! number of coordinate directions of transport xs provided in xs record      &
	                                                             ! 0 : no coordinate dependent transport xs.

	          integer(NBI), pointer, dimension(:)   :: id_scat             ! scattering matrix type identification for scattering block N.
                                                                           ! 000 : total,  100 : elastic,  200 : inelastic,  300 : (n,2n) scattering
	          integer(NBI), pointer, dimension(:)   :: no_scat_order       ! number of scattering orders in block N.
	          integer(NBI), pointer, dimension(:,:) :: no_grp_scat         ! number of groups that scatter into group g.
	          integer(NBI), pointer, dimension(:,:) :: pos_ingrp_scat      ! position of in-group scattering xs in scattering data for group g       

              integer(NBI), pointer, dimension(:)   :: iso_spec            ! spectrum index used to calculate emission spectrum from fission in group g 
           end type Isotope

           contains

           subroutine init_isotope( isotp, const )
              type(XSConst), intent(in)  :: const
	          type(Isotope), pointer, dimension(:) :: isotp
	          integer(NBI) :: i
	          integer(NBI) :: niso, ng, nscmax
          
              niso   = const%no_isotope
	          ng     = const%no_group
	          nscmax = const%max_no_block_sct
            
	          allocate( isotp( niso ) )

              do i=1, niso 
                 isotp(i)%no_group = ng
		         call dmalloc( isotp(i)%id_scat, nscmax )
		         call dmalloc( isotp(i)%no_scat_order, nscmax )
		         call dmalloc( isotp(i)%no_grp_scat, ng, nscmax )
		         call dmalloc( isotp(i)%pos_ingrp_scat, ng, nscmax )

                 isotp(i)%flag_chi    = FLAG_OFF
                 isotp(i)%flag_xs_nf  = FLAG_OFF
                 isotp(i)%flag_xs_na  = FLAG_OFF
                 isotp(i)%flag_xs_np  = FLAG_OFF
                 isotp(i)%flag_xs_nd  = FLAG_OFF
                 isotp(i)%flag_xs_nt  = FLAG_OFF
                 isotp(i)%flag_xs_n2n = FLAG_OFF
              enddo
            
           end subroutine init_isotope


           subroutine allocate_isotope( isotp, const )
              type(XSConst), intent(in)  :: const
	          type(Isotope), intent(out) :: isotp
              integer(NBI) :: g, i, j, k
	          integer(NBI) :: ng, nscmax, nchi, ntrn, ntot, nstrpd
              character(ISTN) :: string

	          ng     = const%no_group
	          nscmax = const%max_no_block_sct

              call trimw6( string, isotp%name )
              isotp%name = string
             
              nchi = isotp%flag_chi
              ntrn = isotp%no_mnt_xs_tr
              ntot = isotp%no_mnt_xs_tt
              nstrpd = isotp%no_dir_xs_tr

              call dmalloc( isotp%pl_xs_tt, ng, ntot )
	          call dmalloc( isotp%pl_xs_tr, ng, ntrn )
	          call dmalloc( isotp%xs_ngm, ng )
	          call dmalloc( isotp%xs_nf, ng )
              call dmalloc( isotp%xs_na, ng )
	          call dmalloc( isotp%xs_np, ng )
	          call dmalloc( isotp%xs_n2n, ng )
	          call dmalloc( isotp%xs_nd, ng )
	          call dmalloc( isotp%xs_nt, ng )
              call dmalloc( isotp%xs_scat, ng, ng )
              call dmalloc( isotp%xs_n2ns, ng, ng )
             
	          call dmalloc( isotp%tnu, ng )
	          call dmalloc( isotp%chiv, ng )
              call dmalloc( isotp%dir_xs_tr, ng, nstrpd )
              call dmalloc( isotp%chim, nchi, ng )
           end subroutine allocate_isotope


        !!
           subroutine delete_isotope( isotp, const )
              type(XSConst), intent(in)  :: const
	          type(Isotope), pointer, dimension(:) :: isotp
 	          integer(NBI) :: niso, i
	          integer(NBI) :: alloc_stat

              niso   = const%no_isotope
              do i=1, niso 
                 deallocate( isotp(i)%pl_xs_tt, stat=alloc_stat )
		         deallocate( isotp(i)%pl_xs_tr, stat=alloc_stat )
		         deallocate( isotp(i)%xs_ngm, stat=alloc_stat )
		         deallocate( isotp(i)%xs_nf,  stat=alloc_stat )
		         deallocate( isotp(i)%xs_na,  stat=alloc_stat )
		         deallocate( isotp(i)%xs_np,  stat=alloc_stat )
		         deallocate( isotp(i)%xs_n2n, stat=alloc_stat )
		         deallocate( isotp(i)%xs_nd,  stat=alloc_stat )
		         deallocate( isotp(i)%xs_nt,  stat=alloc_stat )
             
		         deallocate( isotp(i)%tnu, stat=alloc_stat )
		         deallocate( isotp(i)%chiv, stat=alloc_stat )
                 deallocate( isotp(i)%dir_xs_tr, stat=alloc_stat )

                 deallocate( isotp(i)%chim, stat=alloc_stat )
		         deallocate( isotp(i)%id_scat, stat=alloc_stat )
		         deallocate( isotp(i)%no_scat_order, stat=alloc_stat )
		         deallocate( isotp(i)%no_grp_scat, stat=alloc_stat )
		         deallocate( isotp(i)%pos_ingrp_scat, stat=alloc_stat )
              enddo

              deallocate( isotp, stat=alloc_stat )
              
           end subroutine delete_isotope

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !    Subroutine to read isotxs format and storage to isotope and xs_const structure
        !    Reference to 'isotope.f90' and 'xs_const.f90' about the deifinition of the variables          
        !

           subroutine read_isotxs( isotp, const, fileiso)

              type(XSConst), intent(out) :: const
              type(Isotope), pointer, dimension(:) :: isotp

              logical        :: check = .TRUE.
              character(MCL) :: hname, bname
              character(MCL) :: ctmp1, ctmp2, ctmp3, ctmp4
              integer(NBI)   :: itmp, maxn
              integer(NBI)   :: niso, ng, nscmax, nchi, ltrn, lordn, ltot, kmax, nstrpd
              integer(NBI)   :: nchig, nnfg, nnag, nnpg, nn2ng, nndg, nntg
        ! debug
              integer(NBI), pointer, dimension(:) :: nkmax
              integer(NBI)   :: nkk
        !
              integer(NBI)   :: g, i, j, k, l
              integer(NBI)   :: sblk, eblk
              integer(NBI)   :: naiso
              real(NBF)      :: rtmp
              
              character*72 :: fileiso
           
        !  file open
              naiso=1001
!              open(naiso,file=fileiso, status='old', err=1000)
			open(naiso,file=fileiso, status='old')

        ! 0V
        !      read(naiso,*) hname, maxn, ctmp2, itmp                                                                    
        !!!      read(naiso,*) (itmp, i=1,maxn)
        !      read(naiso,*) itmp
              read(naiso,*) bname                                                                                

              call check_word( bname, '0V', check )
              if( check .eq. .FALSE. )  call write_error( bname )

        ! 1D    
              read(naiso,*) bname,  &                                                                      
                            const%no_group,         const%no_isotope,     const%max_upsct,  &          
                            const%min_dnsct,        const%max_sct_order,  const%flag_chi,   &            
     			             const%max_no_block_sct, const%no_subblock   
   
              call check_word( bname, '1D', check )
              if( check .eq. .FALSE. )  call write_error( bname )

              call init_xsconst( const )
              call init_isotope( isotp, const )

        !  2D    
              read(naiso,*) bname
              call check_word( bname, '2D', check )
              if( check .eq. .FALSE. )  call write_error( bname )

              niso = const%no_isotope
              ng   = const%no_group
        !! isotope name & identification
              read(naiso,'(A8,9(1X,A6)/(10(1X,A6)))') ctmp1, (isotp(i)%name,i=1,niso)
              if( const%flag_chi .eq. FLAG_ON ) then
#ifdef MOD_ISO
				read(naiso,'(5E14.7)') (const%wide_chiv(g),g=1,ng) ! mh
#else
                read(naiso,'(6E12.5)') (const%wide_chiv(g),g=1,ng)
#endif
              endif
#ifdef MOD_ISO
              read(naiso,*)  (const%velocity(g),g=1,ng), (const%max_energy(g),g=1,ng),    
     +	                      const%min_energy  
#else
!              read(naiso,'((6E12.5/6E12.5),(/6E12.5),E12.5/)')  (const%velocity(g),g=1,ng), (const%max_energy(g),g=1,ng),   
              read(naiso,'(6E12.5)')  (const%velocity(g),g=1,ng), (const%max_energy(g),g=1,ng), const%min_energy   			                              
#endif
              read(naiso,'(12I6)') (const%skip_block(i),i=1,niso)

        !  3D
        !! file-wide fission spectrum matrix
              if( const%flag_chi .gt. FLAG_ON ) then
        !         read(naiso,'(6E12.5)') (const%wide_chi_vec(g),g=1,ng)
              endif

              sblk = 0                ! count block               
              do i=1,niso
! FIX_ME
!                 if( sblk .ne. const%skip_block(i) )  call write_error( 'unexpected block' )
        !  4D
                 sblk = sblk + 1
                 read(naiso,*) bname
                 call check_word( bname, '4D', check )
                 if( check .eq. .FALSE. )  call write_error( bname )
        !! isotope constants
                 nscmax = const%max_no_block_sct
         
#ifdef MOD_ISO
	             read(naiso,'(5E14.7)')  isotp(i)%mass,          isotp(i)%kappa,            isotp(i)%energy_yld_capt, &
                                        isotp(i)%temperature,   isotp(i)%potent_scat_res
                 read(naiso,'(E14.7)')   isotp(i)%density                                     ! mh
#else
                 read(naiso,'(6E12.5)')  isotp(i)%mass,          isotp(i)%kappa,            isotp(i)%energy_yld_capt, &
                                         isotp(i)%temperature,   isotp(i)%potent_scat_res,  isotp(i)%density
#endif


                 read(naiso,'(12I6)' )   isotp(i)%class,         isotp(i)%flag_chi,         isotp(i)%flag_xs_nf,      &
                                         isotp(i)%flag_xs_na,    isotp(i)%flag_xs_np,       isotp(i)%flag_xs_n2n,     &
                                         isotp(i)%flag_xs_nd,    isotp(i)%flag_xs_nt,       isotp(i)%no_mnt_xs_tt,    &
                                         isotp(i)%no_mnt_xs_tr,  isotp(i)%no_dir_xs_tr,                               &
     				                         (isotp(i)%id_scat(j),      j=1,nscmax),                                   &  
     						                 (isotp(i)%no_scat_order(j),j=1,nscmax),                                   & 
     						                 ((isotp(i)%no_grp_scat(g,j),g=1,ng),   j=1,nscmax),                       &
     						                 ((isotp(i)%pos_ingrp_scat(g,j),g=1,ng),j=1,nscmax)

        ! allocate
                 call allocate_isotope( isotp(i), const )

        !  5D
                 sblk = sblk + 1
        !! check the length of array
                 nchig = 0
    	         nnfg  = 0
                 nnag  = 0
                 nnpg  = 0
	             nn2ng = 0
                 nndg  = 0
                 nntg  = 0      

                 if( isotp(i)%flag_chi    .eq. FLAG_ON  ) nchig = ng
	             if( isotp(i)%flag_xs_nf  .gt. FLAG_OFF ) nnfg  = ng
	             if( isotp(i)%flag_xs_na  .gt. FLAG_OFF ) nnag  = ng
	             if( isotp(i)%flag_xs_np  .gt. FLAG_OFF ) nnpg  = ng
	             if( isotp(i)%flag_xs_n2n .gt. FLAG_OFF ) nn2ng = ng
	             if( isotp(i)%flag_xs_nd  .gt. FLAG_OFF ) nndg  = ng
                 if( isotp(i)%flag_xs_nt  .gt. FLAG_OFF ) nntg  = ng

	             nstrpd = isotp(i)%no_dir_xs_tr
	             ltrn = isotp(i)%no_mnt_xs_tr
	             ltot = isotp(i)%no_mnt_xs_tt

        !! read isotope xs
#ifdef MOD_ISO
			   read(naiso,'(1X,A2,1X,4E14.7/(5E14.7))') bname,  
#else   
               read(naiso,'(1X,A2,1X,5E12.5/(6E12.5))') bname, &
#endif                                                               
                                                        ((isotp(i)%pl_xs_tr(g,l),g=1,ng),l=1,ltrn),                                &
     	                                                   ((isotp(i)%pl_xs_tt(g,l),g=1,ng),l=1,ltot),                               &
     	                                                   (isotp(i)%xs_ngm(g),g=1,ng),           (isotp(i)%xs_nf(g),g=1,nnfg),      &
     	                                                   (isotp(i)%tnu(g),g=1,nnfg),            (isotp(i)%chiv(g),g=1,nchig),      &
     	                                                   (isotp(i)%xs_na(g),g=1,nnag),          (isotp(i)%xs_np(g),g=1,nnpg),      &
     	                                                   (isotp(i)%xs_n2n(g),g=1,nn2ng),        (isotp(i)%xs_nd(g),g=1,nndg),      &
     	                                                   (isotp(i)%xs_nt(g),g=1,nntg),                                             &
     							                           ((isotp(i)%dir_xs_tr(g,l),g=1,ng),l=1,nstrpd) 

                 call check_word( bname, '5D', check )
                 if( check .eq. .FALSE. )  call write_error( bname )
          
             
        !6D
        !! isotope chi data
                 nchi = isotp(i)%flag_chi
                 if( nchi .gt. FLAG_ON ) then
                    sblk = sblk + 1
                    read(naiso,'(1X,A2,1X,5E12.5/(6E12.5))') bname, ((isotp(i)%chim(k,g),k=1,nchi),g=1,ng)
                    call check_word( bname, '6D', check )
                    if( check .eq. .FALSE. )  call write_error( bname )

                    read(naiso,'(12I6)')  (isotp(i)%iso_spec(g),g=1,ng) 
                 endif

        !7D 
                 call dmalloc(nkmax,nscmax)
                 nkk = 0
                 do k=1, nscmax
                    lordn = isotp(i)%no_scat_order(k)
		            if( lordn .gt. FLAG_OFF ) then
                       sblk = sblk + 1                                

   		               kmax = 0
                       do g=1,ng
	                      kmax = kmax + isotp(i)%no_grp_scat(g,k)
                       enddo
                       nkmax(k) = kmax
                       if( nkk .lt. kmax ) nkk = kmax
                    endif
                 enddo

        !! scattering sub-block
                 call dmalloc( isotp(i)%scatmat, nkk, nscmax )

                 do k=1, nscmax
                    lordn = isotp(i)%no_scat_order(k)
		            if( lordn .gt. FLAG_OFF ) then
#ifdef MOD_ISO
					 read(naiso,'(1X,A2,1X,4E14.7/(5E14.7))') bname, (isotp(i)%scatmat(j,k),j=1,nkmax(k))  ! mh
#else
                       read(naiso,'(1X,A2,1X,5E12.5/(6E12.5))') bname, (isotp(i)%scatmat(j,k),j=1,nkmax(k))
#endif
                       call check_word( bname, '7D', check )
                       if( check .eq. .FALSE. )  call write_error( bname )
                    endif
                 enddo

        !! generate scattering matrix    
                 call make_scat_matrix( isotp(i) )
              enddo

              close(naiso)         
!1000          print *, 'error occured while reading isotxs'

           end subroutine read_isotxs

        !
        !  transform scatmat(ID,NBLK) -> xs_scat(NG,NG)
        !  get sg->tg  ::  xs_scat(sg,tg)
        !
           subroutine make_scat_matrix( isotp )
              type(Isotope)   :: isotp
              integer(NBI)    :: nb, ng, sg, tg, id, nscmat
              integer(NBI)    :: elnb, ienb, n2nb, tsnb
              integer(NBI)    :: igs, lmg, mxg
              real(NBF)       :: fact, ratio, sum_2ns

              ng = isotp%no_group
              nscmat = size(isotp%id_scat)

              do tg=1,ng
                 do sg=1,ng
                    isotp%xs_scat(sg,tg) = 0.0
                 enddo
              enddo

        ! get no of block for scattering type
              n2nb = NO_ID
              elnb = NO_ID
              ienb = NO_ID
              tsnb = NO_ID

              do nb=1,nscmat     
                 if( isotp%no_scat_order(nb) .eq. FLAG_ON ) then
                    select case( isotp%id_scat(nb) )
                       case( ELASTIC_SCAT   )
                           elnb = nb                                   ! no_block of elastic scattering 
                       case( INELASTIC_SCAT )
                           ienb = nb                                   ! no_block of inelastic scattering 
                       case( N2N_SCAT       )
                           n2nb = nb                                   ! no_block of n2n scattering 
                       case( TOTAL_SCAT     ) 
                           tsnb = nb                                   ! no_block of total scattering
                    end select
                 endif
              enddo

        ! make xs_scat matrix
              do nb=1,nscmat
                 if( isotp%no_scat_order(nb) .eq. FLAG_ON ) then
                    id = 1
                    do tg=1,ng
                       igs = isotp%pos_ingrp_scat(tg,nb)
                       mxg = tg + igs - 1                                   ! maximum energy group ( reference to ISOTXS Manual ) 
                       lmg = tg - isotp%no_grp_scat(tg,nb) + igs                                   ! minimum energy group 

                       do sg=mxg,lmg,-1
                          if( (nb .eq. elnb) .or. (nb .eq. ienb) ) then
                             isotp%xs_scat(sg,tg) = isotp%xs_scat(sg,tg) + isotp%scatmat(id,nb) ! elastic and inelastic scattering
                          else if( nb .eq. n2nb ) then
                             isotp%xs_n2ns(sg,tg) = isotp%scatmat(id,nb)                           ! n2n scattering
                          else 
                             isotp%xs_scat(sg,tg) = isotp%scatmat(id,nb)                           ! total scattering 
                          endif
                          id = id + 1
                       enddo
                    enddo
                 endif
              enddo

        ! process n2n scattering
              if( n2nb .ne. NO_ID ) then
                 do tg=1,ng
                    sum_2ns = 0.0
                    fact    = 1.0

                    do sg=1,ng
                       sum_2ns = sum_2ns + isotp%xs_n2ns(sg,tg)          
                    enddo
                    if( (isotp%xs_n2n(tg) .gt. SMALL) .and. (sum_2ns .gt. SMALL) ) then
                       ratio = isotp%xs_n2n(tg)/sum_2ns
                       if( ratio .gt. 0.75 ) fact = 2.0                                ! ratio,0.75, and fact reference to STEP
                    endif
        !!! fact(tg) ???
                    do sg=1,ng
        !               isotp%xs_scat(sg,tg) = isotp%xs_scat(sg,tg) + fact*isotp%xs_n2ns(sg,tg)
                       isotp%xs_scat(sg,tg) = isotp%xs_scat(sg,tg) + 2*isotp%xs_n2ns(sg,tg)
                    enddo
                 enddo
              endif

              do tg=1,ng
                 isotp%xs_na(tg) =  isotp%xs_na(tg) - isotp%xs_n2n(tg) ! mod_n2n
              enddo



           end subroutine make_scat_matrix

        ! debug
           subroutine write_isotxs( isotp, const )

              type(XSConst), intent(in) :: const
              type(Isotope), intent(in), dimension(:) :: isotp

              character(MCL) :: hname
              character(MCL) :: ctmp1, ctmp2, ctmp3, ctmp4
              integer(NBI)   :: itmp, maxn
              integer(NBI)   :: i, j, k
              real(NBF), dimension(MCL) :: rtmp 
         
              open(nxslib,file='xslib.out',status='unknown')

              write(nxslib,*) const%no_group,         const%no_isotope,     const%max_upsct,          &
                            const%min_dnsct,        const%max_sct_order,  const%flag_chi,             &
     	      		              const%max_no_block_sct, const%no_subblock   

              write(*,*)  'complete to write isotxs!'
              close(nxslib)               

           end subroutine write_isotxs

        end module IsotopeMod