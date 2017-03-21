! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB  
        module XSConstMod
                 use param
                 use BasicParam
		         use allocs
!		         implicit none

                 integer, parameter :: NO_FILE_WIDE_SPECTRUM=0, FILE_WIDE_CHI_SPECTRUM=1, FILE_WIDE_CHI_MATRIX=2
                 

           type XSConst
                 
        ! Reference : ISOTXS format manual written in 11/30/76
        ! XSConst   : global constants over total isotopes
                 
                 integer(NBI) :: no_group
		         integer(NBI) :: no_isotope
		         integer(NBI) :: max_upsct
		         integer(NBI) :: min_dnsct 
		         integer(NBI) :: max_sct_order                        ! maximum scattering order
		         integer(NBI) :: flag_chi                             ! 0 : No file-wide spectrum,   1 : file-wide chi vector,   2 : file-wide chi matrix
		         integer(NBI) :: max_no_block_sct                     ! maximum number of blocks of scattering data 
		         integer(NBI) :: no_subblock                          ! subblocking control for scatter matrices.       &
 		                                                              ! the scattering data are subblocked into no_subblock records per scattering block.

		         integer(NBI), pointer, dimension(:) :: skip_block    ! number of records to be skipped to read data for isotope i.

		         real(NBF) ::  min_energy

		         real(NBF), pointer, dimension(:)   :: max_energy
		         real(NBF), pointer, dimension(:)   :: wide_chiv
		         real(NBF), pointer, dimension(:,:) :: wide_chim
		         real(NBF), pointer, dimension(:)   :: velocity


           end type XSConst

           contains

           subroutine init_xsconst( const )
              type(XSConst), intent(inout) :: const
              integer(NBI) :: i, g, ng, nf, niso

              ng = const%no_group
	          nf = const%flag_chi
	          niso = const%no_isotope

              call dmalloc( const%velocity, ng )
	          call dmalloc( const%max_energy, ng )
              call dmalloc( const%skip_block, niso )

              do g=1,ng
	             const%velocity(g)   = 0.0
		         const%max_energy(g) = 0.0
	          enddo
	          do i=1,niso
	             const%skip_block(i) = 0
	          enddo

              select case ( nf )
    	          case ( 0 ) 
	                  continue
	              case ( 1 )
                      call dmalloc( const%wide_chiv, ng )
			          do g=1,ng
			             const%wide_chiv(g) = 0.0
                      enddo
                  case ( 2: ) 
                      call dmalloc( const%wide_chim, nf, ng )
                      do g=1,ng
			             do i=1,nf
			                const%wide_chim(i,g) = 0.0
                         enddo
			          enddo
	          end select

           end subroutine init_xsconst

        end module XSConstMod