! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB  
    module MaterialMod
	use param
	use BasicParam
	use IsotopeMod
	use ErrorHandle
	use XSDlayMod
	use XSConstMod
	use OptionMod
	use allocs
	use trinx_cntl

	implicit none

	integer, parameter :: NTHVAR = 4                                            ! th-variable 1~4
	integer, parameter :: TBASE  = 1
	integer, parameter :: TCURR  = 5                                            ! current temp 5 
	integer, parameter :: TOTNT  = 5                                            ! total number of temp variable array : 4+1
	real(NBF), parameter    :: CONV   = 3.20435306e-11 ! 3.2e-11

! for UNCARDS
	integer, parameter :: NTHVAR_UNC = 3


	type Material
	character(ISTN) :: name
	integer(NBI)    :: no_isotope
	integer(NBI)    :: no_group

	type(Isotope),   pointer, dimension(:,:) :: isotope
	character(ISTN), pointer, dimension(:,:) :: iso_name

	real(NBF), pointer, dimension(:)     :: iso_density                       ! atoms/cc * 1.e-24
	real(NBF), pointer, dimension(:)     :: iso_density0 
	real(NBF), pointer, dimension(:,:)   :: cap_coef
	real(NBF), pointer, dimension(:,:)   :: fis_coef
	real(NBF), pointer, dimension(:,:)   :: nfis_coef
	real(NBF), pointer, dimension(:,:)   :: kfis_coef
	real(NBF), pointer, dimension(:,:)   :: tran_coef
	real(NBF), pointer, dimension(:,:)   :: chip_coef

	real(NBF), pointer, dimension(:,:,:) :: sct_coef


	real(NBF), pointer, dimension(:)     :: fde_den    
	real(NBF), pointer, dimension(:,:)   :: xs_nf
	real(NBF), pointer, dimension(:,:)   :: xs_f
	real(NBF), pointer, dimension(:,:)   :: xs_a
	real(NBF), pointer, dimension(:,:)   :: xs_c
	real(NBF), pointer, dimension(:,:)   :: xs_r
	real(NBF), pointer, dimension(:,:)   :: xs_tr
	real(NBF), pointer, dimension(:,:)   :: xs_kf
	real(NBF), pointer, dimension(:,:,:) :: xs_s
	real(NBF), pointer, dimension(:,:)   :: chip
	real(NBF), pointer, dimension(:,:)   :: chit
	real(NBF), pointer, dimension(:,:)   :: chid

	real(NBF), pointer, dimension(:,:)   :: xs_nf0
	real(NBF), pointer, dimension(:,:)   :: xs_a0
	real(NBF), pointer, dimension(:,:)   :: xs_tr0
	real(NBF), pointer, dimension(:,:)   :: xs_kf0
	real(NBF), pointer, dimension(:,:,:) :: xs_s0

	real(NBF), pointer, dimension(:,:)   :: tnu
	real(NBF), pointer, dimension(:)     :: kappa

	real(NBF), pointer, dimension(:)   :: beta 
	real(NBF), pointer, dimension(:)   :: dct

	real(NBF), pointer, dimension(:)     :: tfuel_var
	integer(NBI), pointer, dimension(:)  :: no_temp_var
	integer(NBI) :: max_temp_var
	real(NBF) :: tcool_base
	real(NBF) :: tcool_var

	logical :: update
	end type Material


	type Composition

	Character(ISTN) :: name
	integer(NBI)    :: no_material
	integer(NBI)    :: no_group

	type(Material),  pointer, dimension(:) :: material
	character(ISTN), pointer, dimension(:) :: mat_name

	real(NBF), pointer, dimension(:) :: vol_fraction
	real(NBF), pointer, dimension(:) :: vol_fraction0

	real(NBF), pointer, dimension(:)   :: xs_nf
	real(NBF), pointer, dimension(:)   :: xs_f
	real(NBF), pointer, dimension(:)   :: xs_a
	real(NBF), pointer, dimension(:)   :: xs_c
	real(NBF), pointer, dimension(:)   :: xs_r
	real(NBF), pointer, dimension(:)   :: xs_tr
	real(NBF), pointer, dimension(:,:) :: xs_s
	real(NBF), pointer, dimension(:)   :: xs_kf

	real(NBF), pointer, dimension(:)   :: xs_nf0
	real(NBF), pointer, dimension(:)   :: xs_f0
	real(NBF), pointer, dimension(:)   :: xs_a0
	real(NBF), pointer, dimension(:)   :: xs_c0
	real(NBF), pointer, dimension(:)   :: xs_r0
	real(NBF), pointer, dimension(:)   :: xs_tr0
	real(NBF), pointer, dimension(:,:) :: xs_s0
	real(NBF), pointer, dimension(:)   :: xs_kf0

	real(NBF), pointer, dimension(:)   :: chip
	real(NBF), pointer, dimension(:)   :: chit
	real(NBF), pointer, dimension(:)   :: chid

	real(NBF), pointer, dimension(:)   :: chip0
	real(NBF), pointer, dimension(:)   :: chit0
	real(NBF), pointer, dimension(:)   :: chid0

	real(NBF), pointer, dimension(:)   :: tnu
	real(NBF) :: kappa

	real(NBF), pointer, dimension(:)   :: beta 
	real(NBF), pointer, dimension(:)   :: dct

! for UNCARDS.
	real(NBF), pointer, dimension(:,:)   :: xsu_nf
	real(NBF), pointer, dimension(:,:)   :: xsu_f
	real(NBF), pointer, dimension(:,:)   :: xsu_r
	real(NBF), pointer, dimension(:,:)   :: xsu_a
	real(NBF), pointer, dimension(:,:)   :: xsu_tr
	real(NBF), pointer, dimension(:,:,:) :: xsu_s
	real(NBF), pointer, dimension(:,:)   :: xsu_kf

!       real(NBF) :: kappa

	end type Composition
       
       
      contains

! read the material information from dif3d input.
!      copy the block '13' of A.NIP3 into a new file, "A.INP"
!

!  Allocate memory
	subroutine init_mat( mat, niso_mat, ng )
	type(Material) :: mat

	integer(NBI), intent(in) :: niso_mat
	integer(NBI) :: ng
	integer(NBI) :: i, m, t
	integer(NBI) :: nmat, mxm

	character(ISTN), parameter :: BLANK = '      '

	mxm = niso_mat
	mat%no_isotope = mxm
	mat%no_group   = ng
	allocate( mat%iso_name(NTHVAR,mxm) )
	allocate( mat%isotope(NTHVAR,mxm)  )
	call dmalloc( mat%no_temp_var, mxm )
	call dmalloc( mat%tfuel_var, NTHVAR )
	call dmalloc( mat%cap_coef,NTHVAR-1,ng )
	call dmalloc( mat%fis_coef,NTHVAR-1,ng )
	call dmalloc( mat%chip_coef,NTHVAR-1,ng )
	call dmalloc( mat%sct_coef,NTHVAR-1,ng,ng )
	call dmalloc( mat%nfis_coef,NTHVAR-1,ng )
	call dmalloc( mat%tran_coef,NTHVAR-1,ng )
	call dmalloc( mat%kfis_coef,NTHVAR-1,ng )
	call dmalloc( mat%fde_den,       mxm )
	call dmalloc( mat%iso_density,   mxm )
	call dmalloc( mat%iso_density0,  mxm )

	call dmalloc( mat%xs_f,   TOTNT, ng )
	call dmalloc( mat%xs_nf,  TOTNT, ng )
	call dmalloc( mat%xs_a,   TOTNT, ng )
	call dmalloc( mat%xs_c,   TOTNT, ng )
	call dmalloc( mat%xs_r,   TOTNT, ng )
	call dmalloc( mat%xs_s,TOTNT,ng, ng )
	call dmalloc( mat%xs_tr,  TOTNT, ng )
	call dmalloc( mat%xs_kf,  TOTNT, ng )
	call dmalloc( mat%tnu,    TOTNT, ng )
	call dmalloc( mat%chip,   TOTNT, ng )
	call dmalloc( mat%chit,   TOTNT, ng )
	call dmalloc( mat%chid,   TOTNT, ng )
	call dmalloc( mat%kappa,  TOTNT     )

	call dmalloc( mat%beta,   NPRECT )
	call dmalloc( mat%dct,    NPRECT )

	do i=1,mxm
	do t=1,NTHVAR
	mat%iso_name(t,i) = BLANK
	enddo
	mat%iso_density(i) = 0.0
	mat%no_temp_var(i) = 1
	enddo
	end subroutine init_mat

	subroutine init_material( mat, comp, niso_mat, nmat_comp, ng )
	type(Material),    pointer, dimension(:) :: mat
	type(Composition), pointer, dimension(:) :: comp

	integer(NBI), intent(in), dimension(:) :: niso_mat
	integer(NBI), intent(in), dimension(:) :: nmat_comp

	integer(NBI) :: ng
	integer(NBI) :: i, m, l, t
	integer(NBI) :: nmat, ncomp, mxm, mxc

	character(ISTN), parameter :: BLANK = '      '

	nmat  = size(niso_mat)
	ncomp = size(nmat_comp)   

	allocate( mat(nmat) )
	allocate( comp(ncomp) )

	do m=1,nmat
		mxm = niso_mat(m)
		mat(m)%no_isotope = mxm
		mat(m)%no_group   = ng
		allocate( mat(m)%iso_name(NTHVAR,mxm) )
		allocate( mat(m)%isotope(NTHVAR,mxm)  )
		call dmalloc( mat(m)%no_temp_var, mxm )
		call dmalloc( mat(m)%tfuel_var, NTHVAR )
		call dmalloc( mat(m)%cap_coef,NTHVAR-1,ng )
		call dmalloc( mat(m)%fis_coef,NTHVAR-1,ng )
		call dmalloc( mat(m)%chip_coef,NTHVAR-1,ng )

		call dmalloc( mat(m)%sct_coef,NTHVAR-1,ng,ng )
		call dmalloc( mat(m)%nfis_coef,NTHVAR-1,ng )
		call dmalloc( mat(m)%tran_coef,NTHVAR-1,ng )
		call dmalloc( mat(m)%kfis_coef,NTHVAR-1,ng )
		call dmalloc( mat(m)%fde_den,       mxm )
		call dmalloc( mat(m)%iso_density,   mxm )
		call dmalloc( mat(m)%iso_density0,  mxm)

		call dmalloc( mat(m)%xs_f,   TOTNT, ng )
		call dmalloc( mat(m)%xs_nf,  TOTNT, ng )
		call dmalloc( mat(m)%xs_a,   TOTNT, ng )
		call dmalloc( mat(m)%xs_c,   TOTNT, ng )
		call dmalloc( mat(m)%xs_r,   TOTNT, ng )
		call dmalloc( mat(m)%xs_s,TOTNT,ng, ng )
		call dmalloc( mat(m)%xs_tr,  TOTNT, ng )
		call dmalloc( mat(m)%xs_kf,  TOTNT, ng )
		call dmalloc( mat(m)%tnu,    TOTNT, ng )
		call dmalloc( mat(m)%chip,   TOTNT, ng )
		call dmalloc( mat(m)%chit,   TOTNT, ng )
		call dmalloc( mat(m)%chid,   TOTNT, ng )
		call dmalloc( mat(m)%kappa,  TOTNT     )

		call dmalloc( mat(m)%beta,   NPRECT )
		call dmalloc( mat(m)%dct,    NPRECT )
		do i=1,mxm
			do t=1,NTHVAR
				mat(m)%iso_name(t,i) = BLANK
			enddo
				mat(m)%iso_density(i) = 0.0
				mat(m)%no_temp_var(i) = 1
		enddo
		mat(m)%update = .FALSE.
	enddo 


	do l=1,ncomp
		mxc = nmat_comp(l)
		comp(l)%no_material = mxc
		comp(l)%no_group    = ng
		comp(l)%kappa       = 0.0
		allocate( comp(l)%mat_name(mxc) )
		allocate( comp(l)%material(mxc) )
		call dmalloc( comp(l)%vol_fraction,mxc )
		call dmalloc( comp(l)%vol_fraction0,mxc )

		call dmalloc( comp(l)%xs_f,   ng )
		call dmalloc( comp(l)%xs_nf,  ng )
		call dmalloc( comp(l)%xs_a,   ng )
		call dmalloc( comp(l)%xs_c,   ng )
		call dmalloc( comp(l)%xs_r,   ng )
		call dmalloc( comp(l)%xs_s,ng,ng )
		call dmalloc( comp(l)%xs_tr,  ng )
		call dmalloc( comp(l)%xs_kf,  ng )

		call dmalloc( comp(l)%xs_f0,   ng )
		call dmalloc( comp(l)%xs_nf0,  ng )
		call dmalloc( comp(l)%xs_a0,   ng )
		call dmalloc( comp(l)%xs_c0,   ng )
		call dmalloc( comp(l)%xs_r0,   ng )
		call dmalloc( comp(l)%xs_s0,ng,ng )
		call dmalloc( comp(l)%xs_tr0,  ng )
		call dmalloc( comp(l)%xs_kf0,  ng )

		call dmalloc( comp(l)%chip,   ng )
		call dmalloc( comp(l)%chit,   ng )
		call dmalloc( comp(l)%chid,   ng )

		call dmalloc( comp(l)%chip0,   ng )
		call dmalloc( comp(l)%chit0,   ng )
		call dmalloc( comp(l)%chid0,   ng )

		call dmalloc( comp(l)%tnu,    ng )
		call dmalloc( comp(l)%beta,   NPRECT )
		call dmalloc( comp(l)%dct,    NPRECT )
! for UNCARDS
		call dmalloc( comp(l)%xsu_f,   ng, NTHVAR_UNC )
		call dmalloc( comp(l)%xsu_nf,  ng, NTHVAR_UNC )
		call dmalloc( comp(l)%xsu_a,   ng, NTHVAR_UNC )
		call dmalloc( comp(l)%xsu_r,   ng, NTHVAR_UNC )
		call dmalloc( comp(l)%xsu_s,ng,ng, NTHVAR_UNC )
		call dmalloc( comp(l)%xsu_tr,  ng, NTHVAR_UNC )
		call dmalloc( comp(l)%xsu_kf,  ng, NTHVAR_UNC )
!
		do i=1,nmat_comp(l)
			comp(l)%mat_name(i) = BLANK
			comp(l)%vol_fraction(i) = 0.0
		enddo
	enddo

	end subroutine init_material

	subroutine reset_macxs_mat( mat )

	type(Material) :: mat
	integer(NBI)   :: t, g, ng, sg

	ng = mat%no_group

	do t=1, TOTNT
		mat%kappa(t) = 0.0
	enddo
	do g=1,ng
		do t=1, TOTNT
			mat%xs_f(t,g)  = 0.0
			mat%xs_nf(t,g) = 0.0
			mat%xs_kf(t,g) = 0.0
			mat%xs_a(t,g)  = 0.0
			mat%xs_c(t,g)  = 0.0
			mat%xs_r(t,g)  = 0.0
			mat%tnu(t,g)   = 0.0
			mat%chip(t,g)   = 0.0
			mat%chit(t,g)   = 0.0
			mat%chid(t,g)   = 0.0
			mat%xs_tr(t,g) = 0.0
		enddo
		do sg=1,ng
			do t=1, TOTNT
			mat%xs_s(t,sg,g) = 0.0
			enddo
		enddo
	enddo
	do g=1,NPRECT
	mat%beta(g)=0.0
	mat%dct(g)=0.0
	enddo
	end subroutine reset_macxs_mat

	subroutine reset_macxs_comp( comp )
	type(Composition) :: comp
	integer(NBI)   :: t, g, ng, sg

	ng = comp%no_group
	do g=1,ng
		comp%xs_f(g)  = 0.0
		comp%xs_nf(g) = 0.0
		comp%xs_kf(g) = 0.0
		comp%xs_a(g)  = 0.0
		comp%xs_r(g)  = 0.0
		comp%xs_c(g)  = 0.0
		comp%tnu(g)   = 0.0
		comp%chip(g)  = 0.0
		comp%chit(g)  = 0.0
		comp%chid(g)  = 0.0
		comp%xs_tr(g) = 0.0
		do sg=1,ng
			comp%xs_s(sg,g) = 0.0
		enddo
	enddo

	do g=1,NPRECT
		comp%beta(g)=0.0
		comp%dct(g)=0.0
	enddo
! for UNCARDS
	do g=1,ng
		do t=1,NTHVAR_UNC
			comp%xsu_f(g,t)  = 0.0
			comp%xsu_nf(g,t) = 0.0
			comp%xsu_kf(g,t) = 0.0
			comp%xsu_r(g,t)  = 0.0
			comp%xsu_a(g,t)  = 0.0
			comp%xsu_tr(g,t) = 0.0
			do sg=1,ng
			comp%xsu_s(sg,g,t) = 0.0
			enddo
		enddo
	enddo

	end subroutine reset_macxs_comp


	subroutine read_material( mat, comp, xs_const, filemat )

	integer, parameter :: NMPL   =  3                             ! number of materials per line
	integer, parameter :: MBLK   = 13
	integer, parameter :: CBLK   = 14
	integer, parameter :: MTFBLK = 15
	integer, parameter :: MTCBLK = 16
	integer, parameter :: TCBLK  = 17
!      integer, parameter :: TFBLK  = 18
	character*72 :: filemat

	type(XSConst)      :: xs_const
	type(Material),    pointer, dimension(:) :: mat
	type(Composition), pointer, dimension(:) :: comp


	logical         :: check
	character(MCL)  :: hname, ctmp1, isname
	character(ISTN) :: mname, pmname
	integer(NBI)    :: iblk, ierr
	integer(NBI)    :: nmat, ncomp, nmax, nmline, ncline, ng, niso, nmiso
	integer(NBI)    :: msum, lsum
	integer(NBI)    :: i, t, mid, lid, lm, lc, mci, lci, m
	real(NBF)       :: rtmp, tcbase

	character(ISTN), dimension(NMPL)    :: iname
	character(ISTN), dimension(NTHVAR)  :: itname

	real(NBF),       dimension(NMPL)    :: frac
	real(NBF),       dimension(NTHVAR)  :: thvar
	integer(NBI), pointer, dimension(:) :: niso_mat,  nmat_comp
	integer(NBI), pointer, dimension(:) :: niso_line, nmat_line
	integer(NBI) :: nainp

	real(NBF) :: fuelv, cladv
	logical :: check_mat	

	character(ISTN), parameter :: BLANK = '      '     

	real(NBF) :: dum(100), fraction
	!  file open
	nainp = 1002
	!open(nainp,file='A.INP')
	open(nainp,file=filemat)
	read(nainp,*) hname                             ! A.INP , TITLE

	nmat   = 0
	ncomp  = 0
	nmline = 0
	ncline = 0
	mname  = BLANK

	!  count number of materials and compositions
	read(nainp,*) hname                             ! A.INP

	do while( .not. EOF(nainp) )
		pmname = mname
		read(nainp,'(I2,10X,A6)',iostat=ierr) iblk, mname
!		read(nainp,'(I,A)',iostat=ierr) iblk, mname ! mod_trinx
		if( ierr .ne. SUCCESS ) continue
	!  block 13
		if( iblk .eq. MBLK  ) then
			if( mname .ne. pmname ) then
				nmat = nmat + 1
			endif
			nmline = nmline + 1
	!  block 14
		else if( iblk .eq. CBLK ) then
			if( mname .ne. pmname ) then
				ncomp = ncomp + 1   
			endif
			ncline = ncline + 1
		else
			continue
		endif         
	enddo

	call dmalloc(niso_mat,  nmat )
	call dmalloc(nmat_comp, ncomp)
	call dmalloc(niso_line, nmline)
	call dmalloc(nmat_line, ncline)

	!  count isotope per materal
	rewind(nainp)
	read(nainp,*) hname                   ! A.INP

	mid  = 0
	lid  = 0
	lm   = 0
	lc   = 0
	mname = BLANK

	do while( .not. EOF(nainp) )
		pmname = mname
		do i=1, NMPL
		iname(i) = BLANK
		enddo

		read(nainp,'(I2,10X,A6,3(A6,E12.5))',iostat=ierr) iblk, mname, ( (iname(i),frac(i)),i=1,NMPL )
!		read(nainp,'(I2,10X,A10,3(A10,E12.5))',iostat=ierr) iblk, mname, ( (iname(i),frac(i)),i=1,NMPL ) ! mod_trinx
		if( ierr .ne. SUCCESS ) continue

		if( iblk .eq. MBLK  ) then
			lm = lm + 1
			if( mname .ne. pmname ) mid = mid + 1
			do i=1,NMPL
			if(iname(i) .ne. BLANK) then
			niso_mat(mid) = niso_mat(mid) + 1
			niso_line(lm) = niso_line(lm) + 1 
			endif
			enddo
		else if( iblk .eq. CBLK ) then
			lc = lc + 1
			if( mname .ne. pmname ) lid = lid + 1
			do i=1,NMPL
			if(iname(i) .ne. BLANK) then
			nmat_comp(lid) = nmat_comp(lid) + 1 
			nmat_line(lc)  = nmat_line(lc)  + 1
			endif
			enddo
		else 
			continue
		endif
	enddo

	! allocate memory to material and composition array

	call init_material( mat, comp, niso_mat, nmat_comp, xs_const%no_group )

	rewind(nainp)
	! read materials and comp
	mid  = 0
	lid  = 0
	lm   = 0
	lc   = 0
	read(nainp,*) hname                 ! A.INP

	do while( .not. EOF(nainp) )
		pmname = mname
		do i=1, NMPL
		iname(i) = BLANK
		enddo

		read(nainp,'(I2,10X,A6,3(A6,E12.5))',iostat=ierr) iblk, mname, ( (iname(i),frac(i)),i=1,NMPL ) 
!		read(nainp,'(I2,10X,A10,3(A10,E12.5))',iostat=ierr) iblk, mname, ( (iname(i),frac(i)),i=1,NMPL )  ! mod_trinx
		if( ierr .ne. SUCCESS ) continue

		if( iblk .eq. MBLK  ) then
			lm = lm + 1                                         ! line index
			if( mname .ne. pmname ) then
			mid = mid + 1                                    ! material index
			mci = 1                                          ! initialize isotope counter
			mat(mid)%name = mname
			endif

			do i=1,niso_line(lm)      
			mat(mid)%iso_name(TBASE,mci) = iname(i)
			mat(mid)%iso_density(mci) = frac(i) 
			mat(mid)%iso_density0(mci) = frac(i)    
			mci = mci + 1                                     ! count isotope
			enddo
		else if( iblk .eq. CBLK ) then
			lc = lc + 1
			if( mname .ne. pmname ) then
			lid = lid + 1
			lci = 1
			comp(lid)%name = mname
			endif

			do i=1,nmat_line(lc)
			comp(lid)%mat_name(lci) = iname(i)
			comp(lid)%vol_fraction(lci) = frac(i)
			comp(lid)%vol_fraction0(lci) = frac(i)
			lci = lci + 1
			enddo
		else 
			continue
		endif
	enddo

	rewind(nainp)

	! read material block 15 for fuel temperature feedback
	read(nainp,*) hname                 ! A.INP
	do while( .not. EOF(nainp) )
		do t=1, NTHVAR
		itname(t) = BLANK
		enddo

		read(nainp,*,iostat=ierr) iblk, mname, ( itname(t),t=1,NTHVAR )      
		if( ierr .ne. SUCCESS ) continue

		if( iblk .eq. MTFBLK  ) then
			check_mat = .false.
			do mid=1, nmat
				if( mname .eq. mat(mid)%name ) then
					check_mat = .true.
					niso = mat(mid)%no_isotope
					do i=1, niso
						if( itname(1) .eq. mat(mid)%iso_name(TBASE,i) ) then
						do t=2,NTHVAR
						mat(mid)%iso_name(t,i) = itname(t)
						enddo
						mat(mid)%no_temp_var(i) = NTHVAR
						endif
					enddo
				endif
			enddo

			lid = 1
			do while( ( check_mat .eq. .false. ) .and. ( lid .le. ncomp ) )
				do m=1, nmat_comp(lid)
					if( comp(lid)%mat_name(m) .eq. mname ) then
					check_mat = .true.
					nmiso = 1
					do t=1, NTHVAR
					comp(lid)%material(m)%iso_name(t,nmiso) = itname(t)        ! if material is isotope, no_iso is one.
					enddo
					endif
				enddo  
				lid = lid + 1
			enddo
			if( check_mat .eq. .false. )  then
				call write_error2('Not found material',mname)
			endif
		endif
	enddo

	rewind(nainp)
	! read material block 16 for coolant temperature feedback
	read(nainp,*) hname                 ! A.INP
	do while( .not. EOF(nainp) )
		read(nainp,*,iostat=ierr) iblk, mname, isname, rtmp     
		if( ierr .ne. SUCCESS ) continue

		if( iblk .eq. MTCBLK  ) then
			do mid=1, nmat
				if( mname .eq. mat(mid)%name ) then
					niso = mat(mid)%no_isotope
					do i=1, niso
					if( itname(1) .eq. mat(mid)%iso_name(TBASE,i) ) then
					mat(mid)%fde_den = rtmp
					mat(mid)%no_temp_var(i) = 1                                     ! coolant
					endif
					enddo
				endif
			enddo
		endif
	enddo

	! read material block 17 for fuel temperature variance
	rewind(nainp)
	read(nainp,*) hname                   ! A.INP

	do while( .not. EOF(nainp) )
		read(nainp,*,iostat=ierr) iblk, tcbase     
		if( ierr .ne. SUCCESS ) continue

		if( iblk .eq. TCBLK  ) then
			do mid=1, nmat
			mat(mid)%tcool_base = tcbase
			enddo
		endif
	enddo

	close(nainp)

	end subroutine read_material



!  Preprocess for materials
	subroutine make_material( mat, isotps )

	type(Isotope),  pointer, dimension(:) :: isotps
	type(Material), pointer, dimension(:) :: mat

	character(ISTN) :: iname
	integer(NBI) :: i, j, l, m, nmat, ntiso, nmiso, t, nth 

	ntiso = size(isotps)
	nmat  = size(mat)

	do m=1, nmat
		nmiso = mat(m)%no_isotope
		do i=1, ntiso
			iname = isotps(i)%name
			do j=1, nmiso
				nth = mat(m)%no_temp_var(j)
				do t=1, nth
					if( mat(m)%iso_name(t,j) .eq. iname ) then
						mat(m)%isotope(t,j) = isotps(i)
					endif
				enddo  
			enddo
		enddo
	enddo

	! set temperature
	do m=1, nmat
		nmiso = mat(m)%no_isotope
		do j=1, nmiso
			nth = mat(m)%no_temp_var(j)
			if( nth .eq. NTHVAR ) then
				do t=1, NTHVAR
				mat(m)%tfuel_var(t) = mat(m)%isotope(t,j)%temperature
				enddo
			else
				mat(m)%tcool_var = mat(m)%isotope(TBASE,j)%temperature
			endif
		enddo
	enddo
	end subroutine make_material


  !  Preprocess for compositions
	subroutine make_composition( comp, mat, isotps, opt )

	type(Option) :: opt
	type(Material),    pointer, dimension(:) :: mat
	type(Composition), pointer, dimension(:) :: comp
	type(Isotope),     pointer, dimension(:) :: isotps

	logical :: check_mat
	character(ISTN) :: mname
	integer(NBI) :: i, j, l, t, m, ncomp, nmat, nmiso, niso, nmat_comp 

	nmat  = size(mat)
	ncomp = size(comp)
	niso  = size(isotps)

	do l=1, ncomp
		nmat_comp = comp(l)%no_material
		do m=1, nmat_comp
			mname = comp(l)%mat_name(m)
			check_mat = .false.
			do j=1, nmat
				if( mat(j)%name .eq. mname ) then
					comp(l)%material(m) = mat(j)
					check_mat = .true.			
				endif  
			enddo
			i = 1
			do while( ( check_mat .eq. .false. ) .and. ( i .le. niso ) )
				if( isotps(i)%name .eq. mname ) then
					nmiso = 1                                                          !! if material is isotope, nmiso is one.
					call init_mat( comp(l)%material(m), nmiso , comp(l)%no_group )
					comp(l)%material(m)%name = mname
					comp(l)%material(m)%isotope(TBASE,nmiso) = isotps(i)
					comp(l)%material(m)%iso_density = 1.0
					check_mat = .true.
				endif
				i = i + 1
			enddo

			if( check_mat .eq. .false. ) then
				call write_error2('Not found material',mname)
			endif
		enddo
	enddo

	end subroutine make_composition


	subroutine make_macxs( comp, const, dlayxs, opt, ifinit)

	logical,intent(in) :: ifinit

	type(XSConst) :: const
	type(XSDlay)  :: dlayxs
	type(Composition) :: comp
	type(Option)  :: opt

	integer(NBI) :: nmat, niso, ncomp
	integer(NBI) :: l,m,ifa

	nmat = comp%no_material
	if(ifinit) then
		do m=1,nmat
			if( opt%shape_option .eq. QUARTIC ) then
				call reset_macxs_mat(comp%material(m))
				call make_temp_coef(comp%material(m))
				call generate_macxs_mat_temp(comp%material(m), const, dlayxs, opt, TRUE )
			else if( opt%shape_option .eq. FLAT ) then
				call generate_macxs_mat(comp%material(m), const, dlayxs, opt )
			endif
		enddo
		call generate_macxs_comp(comp,TRUE) 
	else
		do m=1,nmat
			if( opt%shape_option .eq. QUARTIC ) then
				do ifa=1,nfuelassm
				if(comp%mat_name(m).eq.fuelmat(ifa,1)) then 
					call generate_macxs_mat_temp(comp%material(m), const, dlayxs, opt, FALSE )	
				endif
				enddo	
			else if( opt%shape_option .eq. FLAT ) then
				call generate_macxs_mat(comp%material(m), const, dlayxs, opt )
			endif
		enddo
		call generate_macxs_comp(comp,FALSE) 	   
	endif ! ifinit
	

	return
	end subroutine make_macxs

! Make temperature coefficients for capture and fission xs 
	subroutine make_temp_coef( mat )
	type(Material) :: mat

	integer(NBI) :: i, t, niso, nvar, mnvar, ig
	integer(NBI) :: ng, sg, g
	real(NBF)    :: tfbase, tcbase, pert, det, rdet
	real(NBF)    :: a11, a12, a13, a21, a22, a23, a31, a32, a33
	real(NBF), dimension(NTHVAR-1) :: TC1, TC2, TC3, TC4, DXSF, DXSC, DXSNF, DXSKF, DXSTR, DXSCHIP      ! coeff(3)
	real(NBF), pointer, dimension(:,:) :: DXSSCT    ! coeff(3,NG)

	real(NBF)    :: sum_mac_chi, sum_mac_xs_nf, sum_iso_xs_nf

	ng   = mat%no_group
	niso = mat%no_isotope 

	do g=1,ng
		do t=1, NTHVAR
			mat%xs_f(t,g)  = 0.0
			mat%xs_nf(t,g)  = 0.0
			mat%xs_kf(t,g)  = 0.0
			mat%xs_tr(t,g)  = 0.0

			mat%xs_a(t,g)  = 0.0
			mat%xs_c(t,g)  = 0.0

			do sg=1,ng
			mat%xs_s(t,sg,g) = 0.0
			enddo
		enddo
	enddo

	do t=1, NTHVAR-1
		TC1(t)   = 0.0
		TC2(t)   = 0.0
		TC3(t)   = 0.0
		TC4(t)   = 0.0
		DXSF(t)  = 0.0
		DXSC(t)  = 0.0
		DXSTR(t) = 0.0
		DXSNF(t) = 0.0
		DXSKF(t) = 0.0
		DXSCHIP(t) = 0.0 
	enddo

	if(.not.associated(DXSSCT)) call dmalloc( DXSSCT, NTHVAR-1, ng)

	mnvar = 0
	do i=1, niso
		nvar = mat%no_temp_var(i)
		if( nvar .gt. mnvar ) mnvar = nvar 
	enddo

	mat%max_temp_var = mnvar

	if( mnvar .eq. NTHVAR ) then
		do g=1, ng
		do t=1, NTHVAR ! t=1 : BASE temperature, t=2~4 : perturbated temperature
			do i=1, niso							
				mat%xs_f(t,g)  = mat%xs_f(t,g)  +  mat%isotope(t,i)%xs_nf(g)*mat%iso_density(i)
				mat%xs_nf(t,g) = mat%xs_nf(t,g) +  mat%isotope(t,i)%xs_nf(g)*mat%isotope(t,i)%tnu(g)*mat%iso_density(i)
				mat%xs_kf(t,g) = mat%xs_kf(t,g) +  mat%isotope(t,i)%xs_nf(g)*mat%isotope(t,i)%kappa*mat%iso_density(i)
				mat%xs_tr(t,g) = mat%xs_tr(t,g) +  mat%isotope(t,i)%pl_xs_tr(g,1)*mat%iso_density(i)

				mat%xs_c(t,g)  = mat%xs_c(t,g)  + (mat%isotope(t,i)%xs_na(g) &
     				+  mat%isotope(t,i)%xs_np(g) + mat%isotope(t,i)%xs_nt(g)  &
     				+  mat%isotope(t,i)%xs_nd(g) + mat%isotope(t,i)%xs_ngm(g))*mat%iso_density(i)
				mat%xs_a(t,g)  = mat%xs_c(t,g)  + mat%xs_f(t,g) 
				do sg=1,ng
				mat%xs_s(t,sg,g) = mat%xs_s(t,sg,g) + mat%isotope(t,i)%xs_scat(sg,g)*mat%iso_density(i)
				enddo
			enddo

			sum_mac_chi   = 0.0
			sum_mac_xs_nf = 0.0
			do i=1,niso
				sum_iso_xs_nf = 0.0
				do ig=1,ng
					sum_iso_xs_nf = sum_iso_xs_nf + mat%isotope(t,i)%xs_nf(ig)*mat%isotope(t,i)%tnu(ig)
				enddo
				sum_mac_chi = sum_mac_chi + mat%isotope(t,i)%chiv(g)*mat%iso_density(i)*sum_iso_xs_nf
				sum_mac_xs_nf = sum_mac_xs_nf + mat%iso_density(i)*sum_iso_xs_nf
			enddo
			mat%chip(t,g) = sum_mac_chi/sum_mac_xs_nf
		enddo
		enddo
	else
		do i=1, niso
			do g=1, ng
! t = TBASE
				mat%xs_f(TBASE,g)  = mat%xs_f(TBASE,g)  +  mat%isotope(TBASE,i)%xs_nf(g)*mat%iso_density(i)
				mat%xs_nf(TBASE,g) = mat%xs_nf(TBASE,g) &
     				+  mat%isotope(TBASE,i)%xs_nf(g)*mat%isotope(TBASE,i)%tnu(g)*mat%iso_density(i)
				mat%xs_kf(TBASE,g) = mat%xs_kf(TBASE,g) &
     				+  mat%isotope(TBASE,i)%xs_nf(g)*mat%isotope(TBASE,i)%kappa*mat%iso_density(i)
				mat%xs_tr(TBASE,g) = mat%xs_tr(TBASE,g) +  mat%isotope(TBASE,i)%pl_xs_tr(g,1)*mat%iso_density(i)

				mat%xs_c(TBASE,g)  = mat%xs_c(TBASE,g)  + (mat%isotope(TBASE,i)%xs_na(g) &
     				+  mat%isotope(TBASE,i)%xs_np(g) + mat%isotope(TBASE,i)%xs_nt(g)  &
     				+  mat%isotope(TBASE,i)%xs_nd(g) + mat%isotope(TBASE,i)%xs_ngm(g))*mat%iso_density(i)
				mat%xs_a(TBASE,g)  = mat%xs_c(TBASE,g)  + mat%xs_f(TBASE,g)

				do sg=1,ng
					mat%xs_s(TBASE,sg,g) = mat%xs_s(TBASE,sg,g) + mat%isotope(TBASE,i)%xs_scat(sg,g)*mat%iso_density(i)
				enddo

				mat%xs_f(TCURR,g) = mat%xs_f(TBASE,g)
				mat%xs_nf(TCURR,g)= mat%xs_nf(TBASE,g)
				mat%xs_kf(TCURR,g)= mat%xs_kf(TBASE,g)
				mat%xs_tr(TCURR,g)= mat%xs_tr(TBASE,g)
				mat%xs_c(TCURR,g) = mat%xs_c(TBASE,g)				
				mat%xs_a(TCURR,g) = mat%xs_a(TBASE,g)
				do sg = 1, ng
					mat%xs_s(TCURR,sg,g) = mat%xs_s(TBASE,sg,g)
				enddo

! t = 2 for coolant 
#ifdef T2COEFF
				mat%xs_f(2,g)  = mat%xs_f(2,g)  +  mat%isotope(TBASE,i)%xs_nf(g)*mat%fde_den(i)
				mat%xs_tr(2,g) = mat%xs_tr(2,g) +  mat%isotope(TBASE,i)%pl_xs_tr(g,1)*mat%fde_den(i)
				mat%xs_nf(2,g) = mat%xs_nf(2,g) +  mat%isotope(TBASE,i)%xs_nf(g)*mat%isotope(TBASE,i)%tnu(g)*mat%fde_den(i)
				mat%xs_kf(2,g) = mat%xs_kf(2,g) +  mat%isotope(TBASE,i)%xs_nf(g)*mat%isotope(TBASE,i)%kappa*mat%fde_den(i)

				mat%xs_c(2,g)  = mat%xs_c(2,g)  + (mat%isotope(TBASE,i)%xs_na(g) &
     				+ mat%isotope(TBASE,i)%xs_np(g) + mat%isotope(TBASE,i)%xs_nt(g)     &
     				+  mat%isotope(TBASE,i)%xs_nd(g) + mat%isotope(TBASE,i)%xs_ngm(g))*mat%fde_den(i)
				mat%xs_a(2,g)  = mat%xs_c(2,g)  + mat%xs_f(2,g)    
				do sg=1,ng
				mat%xs_s(2,sg,g) = mat%xs_s(2,sg,g) + mat%isotope(TBASE,i)%xs_scat(sg,g)*mat%fde_den(i)
				enddo
#endif
			enddo
		enddo
	endif

	if( mnvar .eq. NTHVAR ) then
		tfbase = mat%tfuel_var(TBASE) ! In fuel, TBASE is t=1.
		do g=1, ng
			do t=2, NTHVAR  ! t=1 : BASE temperature, t=2~4 : perturbated temperature
				pert = mat%tfuel_var(t)
				TC1(t-1) = sqrt(pert) - sqrt(tfbase)
				TC2(t-1) = log(pert/tfbase)
				TC3(t-1) = 1/pert - 1/tfbase
				DXSF(t-1)   = mat%xs_f(t,g)  - mat%xs_f(TBASE,g)
				DXSC(t-1)   = mat%xs_c(t,g)  - mat%xs_c(TBASE,g)
				DXSNF(t-1)  = mat%xs_nf(t,g) - mat%xs_nf(TBASE,g)
				DXSKF(t-1)  = mat%xs_kf(t,g) - mat%xs_kf(TBASE,g)
				DXSTR(t-1)  = mat%xs_tr(t,g) - mat%xs_tr(TBASE,g)

				DXSCHIP(t-1)= mat%chip(t,g)  - mat%chip(TBASE,g) ! XSCHI
				do sg=1,ng
					DXSSCT(t-1,sg) = mat%xs_s(t,sg,g)  - mat%xs_s(TBASE,sg,g)
				enddo
			enddo
			det = - TC1(3)*TC2(2)*TC3(1) + TC1(2)*TC2(3)*TC3(1) + TC1(3)*TC2(1)*TC3(2)      &               
     			  - TC1(1)*TC2(3)*TC3(2) - TC1(2)*TC2(1)*TC3(3) + TC1(1)*TC2(2)*TC3(3)
			rdet=1./det

			a11 = (TC2(2)*TC3(3) - TC2(3)*TC3(2))*rdet
			a12 = (TC2(3)*TC3(1) - TC2(1)*TC3(3))*rdet
			a13 = (TC2(1)*TC3(2) - TC2(2)*TC3(1))*rdet
			a21 = (TC1(3)*TC3(2) - TC1(2)*TC3(3))*rdet
			a22 = (TC1(1)*TC3(3) - TC1(3)*TC3(1))*rdet
			a23 = (TC1(2)*TC3(1) - TC1(1)*TC3(2))*rdet
			a31 = (TC1(2)*TC2(3) - TC1(3)*TC2(2))*rdet
			a32 = (TC1(3)*TC2(1) - TC1(1)*TC2(3))*rdet
			a33 = (TC1(1)*TC2(2) - TC1(2)*TC2(1))*rdet

			mat%fis_coef(1,g)  = a11*DXSF(1)  + a12*DXSF(2)  + a13*DXSF(3)
			mat%fis_coef(2,g)  = a21*DXSF(1)  + a22*DXSF(2)  + a23*DXSF(3)
			mat%fis_coef(3,g)  = a31*DXSF(1)  + a32*DXSF(2)  + a33*DXSF(3)

			mat%nfis_coef(1,g) = a11*DXSNF(1) + a12*DXSNF(2) + a13*DXSNF(3)
			mat%nfis_coef(2,g) = a21*DXSNF(1) + a22*DXSNF(2) + a23*DXSNF(3)
			mat%nfis_coef(3,g) = a31*DXSNF(1) + a32*DXSNF(2) + a33*DXSNF(3)

			mat%kfis_coef(1,g) = a11*DXSKF(1) + a12*DXSKF(2) + a13*DXSKF(3)
			mat%kfis_coef(2,g) = a21*DXSKF(1) + a22*DXSKF(2) + a23*DXSKF(3)
			mat%kfis_coef(3,g) = a31*DXSKF(1) + a32*DXSKF(2) + a33*DXSKF(3)

			mat%tran_coef(1,g) = a11*DXSTR(1) + a12*DXSTR(2) + a13*DXSTR(3)
			mat%tran_coef(2,g) = a21*DXSTR(1) + a22*DXSTR(2) + a23*DXSTR(3)
			mat%tran_coef(3,g) = a31*DXSTR(1) + a32*DXSTR(2) + a33*DXSTR(3)

			mat%cap_coef(1,g)  = a11*DXSC(1)  + a12*DXSC(2)  + a13*DXSC(3)
			mat%cap_coef(2,g)  = a21*DXSC(1)  + a22*DXSC(2)  + a23*DXSC(3)
			mat%cap_coef(3,g)  = a31*DXSC(1)  + a32*DXSC(2)  + a33*DXSC(3)

			mat%chip_coef(1,g)  = a11*DXSCHIP(1)  + a12*DXSCHIP(2)  + a13*DXSCHIP(3)
			mat%chip_coef(2,g)  = a21*DXSCHIP(1)  + a22*DXSCHIP(2)  + a23*DXSCHIP(3)
			mat%chip_coef(3,g)  = a31*DXSCHIP(1)  + a32*DXSCHIP(2)  + a33*DXSCHIP(3)

			do sg=1,ng
			mat%sct_coef(1,sg,g)  = a11*DXSSCT(1,sg)  + a12*DXSSCT(2,sg)  + a13*DXSSCT(3,sg)
			mat%sct_coef(2,sg,g)  = a21*DXSSCT(1,sg)  + a22*DXSSCT(2,sg)  + a23*DXSSCT(3,sg)
			mat%sct_coef(3,sg,g)  = a31*DXSSCT(1,sg)  + a32*DXSSCT(2,sg)  + a33*DXSSCT(3,sg)
			enddo
		enddo
	else 
		do g=1, ng
#ifdef T2COEFF
			mat%fis_coef(1,g)  = mat%xs_f(2,g)                                  !  D = SUM ( sig_x*fde_den )
			mat%nfis_coef(1,g) = mat%xs_nf(2,g)
			mat%kfis_coef(1,g) = mat%xs_kf(2,g)
			mat%cap_coef(1,g)  = mat%xs_c(2,g)                                  !  fde_den : first derivative density      
			do sg=1, ng  
			mat%sct_coef(1,sg,g) = mat%xs_s(2,sg,g)                      !  fde_den : first derivative density         
			enddo 
#endif
		enddo
	endif

	return
	end subroutine make_temp_coef

!------------------------------------------------------------------------------------------------------- 
!  define reaction
!    absorption : (n,g)+(n,a)+(n,p)+(n,t)+(n,d)+(n,f)
!    capture    : (n,g)+(n,a)+(n,p)+(n,t)+(n,d)
!    fission    : (n,f)
!    scattering : (n,n)+(n,n')+(n,2ns)  :: elastic or inelastic

!!!      SCAPT(M,N,K)=SCAPT(M,N,K)+(SNGAM(N)+SNALF(N)+SNP(N)+SND(N)+SNT(N))*DEN                              
!!!      SFIST(M,N,K)=SFIST(M,N,K)+(SFIS(N))*DEN           
!!   SCATTERING XS      
!       IF TOTAL SCATTERING MATRIX IS PRESENT, USE IT TO COMPUTE REMOVAL  
!       CROSS SECTION                                                     
!                                                                        
!       SET FACN2N DEPENDING ON NORMALIZATION OF (N,2N)SCATTERING.        
!       FACN2N=1 IF PRODUCTION NORMALIZED.  FACN2N=2 IF REACTION          
!       NORMALIZED                                                        

!!   REMOVAL XS 
!       IF TOTAL SCATTERING IS NOT PRESENT AND SN2N IS PRESENT, USE       
!       IT TO SUBTRACT OUT THE N2N SELF-SCATTERING                        
!                                                                       
!       IF SN2N IS NOT PRESENT, MAKE SAME CORRECTION USING N2N SCATTERING 
!       MATRIX (N2NSCT(N,N))                                              
!                                                                       
!       IF TOTAL SCATTERING MATRIX IS PRESENT, ASSUME THAT N2N SCATTERING 
!       MATRIX IS NOT,  SO USE N2N REACTION CROSS SECTION (SN2N)          
!       TO SUBTRACT OUT 1*N2NSCT(L,N), L.GT.N                             
!       THIS IS EQUIVALENT TO ASSUMING THAT THERE IS NO N2N SELFSCATTERING

!!!      SIGTR(M,N)=SIGTR(M,N)+STRPL(N,1)*DEN                              

!  macro xs for material dependent on temperature
!------------------------------------------------------------------------------------------------------- 
	subroutine generate_macxs_mat_temp( mat, const, dlayxs, opt, ifinit )

	type(XSConst)  :: const
	type(Material) :: mat
	type(XSDlay)   :: dlayxs
	type(Option)   :: opt
	logical,intent(in) :: ifinit

	integer(NBI)   :: i,  t, ig, k, sg, g, mi, m
	integer(NBI)   :: nkfam, niso, ng, nvar, ndiso
	real(NBF)      :: sum_chi, fact, sum_2ns, max_kf
	real(NBF)      :: tfuel, tcool, tfbase, tcbase
	real(NBF)      :: sqrt_temp, log_temp, inv_temp, dif_temp
	real(NBF)      :: sum_mac_xs_nf, sum_iso_xs_nf, sum_mac_chi, sum_dnu
	real(NBF)      :: sum_beta, sum_chid, sum_tnu, sum_chip, sum_chit

#define xsdlay
#ifdef xsdlay
	real(NBf), pointer :: beta(:,:), dct(:), nuf(:), chidnuf(:,:), chid(:)
	integer(NBI) :: eff, gd, TEST
	real(NBF) :: sumu, suml
#else
	real(NBF) :: beta			
#endif
	real(NBF), pointer, dimension(:) :: dnu
	logical, save :: ifalloc=TRUE
	character*4 :: cdum1, cdum2

	niso  = mat%no_isotope
	ng    = mat%no_group

	tfuel = opt%fuel_temp
	tcool = opt%cool_temp

#ifdef TEMP
	if(mat%xs_nf(1,1).gt.0) then
	if(tfuel.gt.mat%tfuel_var(NTHVAR)) then
		write(mesg,'(a50)')  'Wanring ! : A fuel temperature is out of range'
		call message(true,true,mesg)
		tfuel=mat%tfuel_var(NTHVAR)
	elseif(tfuel.lt.mat%tfuel_var(TBASE)) then
		write(mesg,'(a50)')  'Wanring ! : A fuel temperature is out of range'
		call message(true,true,mesg)
		tfuel=mat%tfuel_var(TBASE)
	endif
	endif
#endif

	   
!	call reset_macxs_mat( mat ) 
!	call make_temp_coef( mat )

	tfbase = mat%tfuel_var(TBASE) 
	tcbase = mat%tcool_var
	sqrt_temp = sqrt(tfuel) - sqrt(tfbase)
	log_temp  = log(tfuel/tfbase)
	inv_temp  = 1/tfuel - 1/tfbase
	if( mat%max_temp_var .eq. NTHVAR ) then
		do g=1,ng

			mat%xs_f(TCURR,g)  = mat%xs_f(TBASE,g)  + mat%fis_coef(1,g)*sqrt_temp  + mat%fis_coef(2,g)*log_temp  &
     			+ mat%fis_coef(3,g)*inv_temp
			mat%xs_nf(TCURR,g) = mat%xs_nf(TBASE,g) + mat%nfis_coef(1,g)*sqrt_temp + mat%nfis_coef(2,g)*log_temp &
     			+ mat%nfis_coef(3,g)*inv_temp
			mat%xs_kf(TCURR,g) = mat%xs_kf(TBASE,g) + mat%kfis_coef(1,g)*sqrt_temp + mat%kfis_coef(2,g)*log_temp &
     			+ mat%kfis_coef(3,g)*inv_temp
			mat%xs_tr(TCURR,g) = mat%xs_tr(TBASE,g) + mat%tran_coef(1,g)*sqrt_temp + mat%tran_coef(2,g)*log_temp &
     			+ mat%tran_coef(3,g)*inv_temp

			mat%xs_c(TCURR,g)  = mat%xs_c(TBASE,g)  + mat%cap_coef(1,g)*sqrt_temp  + mat%cap_coef(2,g)*log_temp  &
     			+ mat%cap_coef(3,g)*inv_temp

			mat%xs_a(TCURR,g)  = mat%xs_c(TCURR,g)  + mat%xs_f(TCURR,g)

			mat%chip(TCURR,g) = mat%chip(TBASE,g) + mat%chip_coef(1,g)*sqrt_temp + mat%chip_coef(2,g)*log_temp &
     			+ mat%chip_coef(3,g)*inv_temp

			if(mat%xs_nf(TCURR,g).gt.0.0) mat%update=.TRUE.

			do sg=1,ng
			mat%xs_s(TCURR,sg,g) = mat%xs_s(TBASE,sg,g) + mat%sct_coef(1,sg,g)*sqrt_temp + mat%sct_coef(2,sg,g)*log_temp &
     			+ mat%sct_coef(3,sg,g)*inv_temp
			enddo
		enddo
	else
		do g=1,ng
			dif_temp  = tcool - tcbase
			mat%xs_f(TCURR,g)  = mat%xs_f(TBASE,g)  + mat%fis_coef(1,g)*dif_temp
			mat%xs_nf(TCURR,g) = mat%xs_nf(TBASE,g) + mat%nfis_coef(1,g)*dif_temp
			mat%xs_kf(TCURR,g) = mat%xs_kf(TBASE,g) + mat%kfis_coef(1,g)*dif_temp
			mat%xs_tr(TCURR,g) = mat%xs_tr(TBASE,g) + mat%tran_coef(1,g)*dif_temp

			mat%xs_c(TCURR,g)  = mat%xs_c(TBASE,g)  + mat%cap_coef(1,g)*dif_temp
			mat%xs_a(TCURR,g)  = mat%xs_c(TCURR,g)  + mat%xs_f(TCURR,g)
			do sg=1,ng
			mat%xs_s(TCURR,sg,g)  = mat%xs_s(TBASE,sg,g) + mat%sct_coef(1,sg,g)*dif_temp
			enddo
		enddo
	endif

	if(.not.ifinit) goto 100
!	goto 100

! tnu, dnu
#ifdef xsdlay
	mat%tnu(TBASE,:)=0.0
	do g=1,ng
		mat%tnu(TBASE,g)=mat%tnu(TBASE,g)+mat%xs_nf(TBASE,g)
!		mat%tnu(TBASE,g)=mat%tnu(TBASE,g)+mat%xs_nf(TCURR,g) 		 
	enddo              
#else							 
	do i=1,niso
		nvar = mat%no_temp_var(i)
		do g=1,ng
		do t=1, nvar
			mat%tnu(t,g)   = mat%tnu(t,g)   + mat%isotope(t,i)%tnu(g)*mat%iso_density(i)
		enddo
		enddo
	enddo
#endif


	TEST=1
! chi_total
	if( opt%dly_flag .eq. 2 ) then
		nkfam = dlayxs%no_family
		ndiso = dlayxs%no_isotope
		
#ifdef xsdlay
		eff=ndiso+1

		if(ifalloc) then
			call dmalloc(dnu,nkfam)	
			call dmalloc(beta,nkfam,eff)
			call dmalloc(dct,nkfam)
			call dmalloc(nuf,ndiso)
			call dmalloc(chidnuf,ng,ndiso)
			call dmalloc(chid,ng)
			ifalloc=FALSE
		endif

		sum_dnu=0.0
		do k=1,nkfam
			dnu(k) = 0.0
			do g=1,ng
			do i=1,ndiso
			do mi=1, niso
				cdum1=dlayxs%iso_name(i)
				cdum2=mat%iso_name(TBASE,mi)
				if( cdum1 .eq. cdum2 ) then
!				if( dlayxs%iso_name(i) .eq. mat%iso_name(TBASE,mi) ) then ! mod_xsdlay
				dnu(k) = dnu(k) + dlayxs%dnu(g,k,i)*mat%iso_density(mi)* mat%isotope(TEST ,mi)%xs_nf(g) 
				endif
			enddo
			enddo
			enddo
			sum_dnu=sum_dnu+dnu(k)
		enddo

		sum_tnu=0.0
		do g=1,ng
			sum_tnu=sum_tnu+mat%xs_nf(TCURR ,g)
!			sum_tnu=sum_tnu+mat%tnu(TCURR,g)
		enddo

		sum_beta=0.0
		do k=1,nkfam
			if( sum_tnu .lt. NEAR_ZERO ) then
				mat%beta(k) = 0.0
			else
				mat%beta(k) = dnu(k)/(sum_tnu+sum_dnu)      ! if dlayxs data exists, tnu is not total nu but prompt nu.
			endif
			sum_beta = sum_beta + mat%beta(k)
!			sum_chid = sum_chid + dlayxs%chid(g,k)*beta
!			mat%chit(TBASE,g) = (1.0-sum_beta)*mat%chip(TBASE,g) + sum_chid
		enddo
	!            
		nuf=0.0
		chidnuf=0.0
		do i=1,ndiso
			sum_dnu=0.0
			sum_tnu=0.0
			do k=1,nkfam
				dnu(k) = 0.0
				do g=1,ng
				do mi=1, niso
				cdum1=dlayxs%iso_name(i)
				cdum2=mat%iso_name(TBASE,mi)
				if( cdum1 .eq. cdum2 ) then
!				if( dlayxs%iso_name(i) .eq. mat%iso_name(TBASE,mi) ) then ! mod_xsdlay
					dnu(k) = dnu(k)+dlayxs%dnu(g,k,i)
					if(k.eq.1) then
						sum_tnu = sum_tnu+mat%isotope(TEST,mi)%tnu(g)
						nuf(i) = nuf(i)+mat%isotope(TEST,mi)%tnu(g)*mat%isotope(TEST,mi)%xs_nf(g)
					endif
				endif
				enddo
				enddo
				do gd=1,ng
					chidnuf(gd,i)=chidnuf(gd,i)+dlayxs%chid_ext(gd,k,i)*nuf(i)
				enddo
				sum_dnu=sum_dnu+dnu(k)
			enddo
			do k=1,nkfam
				if( sum_tnu .lt. NEAR_ZERO ) then
					beta(k,i) = 0.0
				else
					beta(k,i)=dnu(k)/(sum_dnu+sum_tnu)     ! if dlayxs data exists, tnu is not total nu but prompt nu.
				endif
			enddo
		enddo


		do g=1,ng
			sumu=0.0
			suml=0.0
			do i=1,ndiso
				sumu=sumu+chidnuf(g,i)
				suml=suml+nuf(i)
			enddo
			if(suml.lt.NEAR_ZERO) then
				mat%chid(TBASE,g)=0.0
			else
				mat%chid(TBASE,g)=sumu/suml
			endif
		enddo

		do k=1,nkfam
			sumu=0.0
			suml=0.0
			do i=1,ndiso
				sumu=sumu+dlayxs%xs_dct_ext(k,i)*beta(k,i)*nuf(i)
				suml=suml+beta(k,i)*nuf(i)
			enddo

			if( sumu .lt. NEAR_ZERO ) then
				mat%dct(k) = 0.0
			else
				mat%dct(k)=sumu/suml     
			endif	
		enddo
#else
		do g=1,ng
			sum_beta = 0.0
			sum_chid = 0.0
			sum_dnu = 0.0

			do k=1,nkfam
				dnu(k) = 0.0
				do i=1,ndiso
				do mi=1, niso
					cdum1=dlayxs%iso_name(i)
					cdum2=mat%iso_name(TBASE,mi)
					if( cdum1 .eq. cdum2 ) then
!					if( dlayxs%iso_name(i) .eq. mat%iso_name(TBASE,mi) ) then ! mod_xsdlay
						dnu(k) = dnu(k) + dlayxs%dnu(g,k,i)*mat%iso_density(mi)           !!!!!! ndiso ???
					endif
				enddo
				enddo
				sum_dnu = sum_dnu + dnu(k)
			enddo

			do k=1,nkfam
				if( mat%tnu(TBASE,g) .lt. NEAR_ZERO ) then
					beta = 0.0
				else
					beta = dnu(k)/(mat%tnu(TBASE,g)+sum_dnu)      ! if dlayxs data exists, tnu is not total nu but prompt nu.
				endif
				sum_beta = sum_beta + beta
				sum_chid = sum_chid + dlayxs%chid(g,k)*beta
			enddo
			mat%chit(TBASE,g) = (1.0-sum_beta)*mat%chip(TBASE,g) + sum_chid
		enddo
#endif
	elseif ( opt%dly_flag .eq. FlAG_OFF ) then
		do g=1,ng
			mat%chit(TBASE,g) = mat%chip(TBASE,g)
		enddo  
	endif

! normalize chi
	sum_chip = 0.0
	sum_chit = 0.0
	sum_chid = 0.0
	do g=1,ng
		sum_chip = sum_chip + mat%chip(TBASE,g)
		sum_chit = sum_chit + mat%chit(TBASE,g)
		sum_chid = sum_chid + mat%chid(TBASE,g)
	enddo
	do g=1,ng
		if( sum_chip .lt. NEAR_ZERO )  sum_chip = NEAR_INF
		mat%chip(TBASE,g) = mat%chip(TBASE,g)/sum_chip
		mat%chip(TCURR,g) = mat%chip(TBASE,g)

		if( sum_chit .lt. NEAR_ZERO )  sum_chit = NEAR_INF
		mat%chit(TBASE,g) = mat%chit(TBASE,g)/sum_chit
		mat%chit(TCURR,g) = mat%chit(TBASE,g)

		if( sum_chid .lt. NEAR_ZERO )  sum_chid = NEAR_INF
		mat%chid(TBASE,g) = mat%chid(TBASE,g)/sum_chid
		mat%chid(TCURR,g) = mat%chid(TBASE,g)
	enddo

100   continue
! correction sigma_kappa_f
	max_kf = 0.0
	do g=1,ng
		if( mat%xs_kf(TCURR,g) .gt. max_kf ) max_kf = mat%xs_kf(TCURR,g)
	enddo

	if( max_kf .lt. NEAR_ZERO ) then
		do g=1,ng
			mat%xs_kf(TCURR,g) = mat%xs_f(TCURR,g)*CONV
		enddo
	endif

	return
	end subroutine generate_macxs_mat_temp

  !
  !  macro xs according to given base temperature
       subroutine generate_macxs_mat( mat, const, dlayxs, opt )
          type(XSConst)  :: const
          type(Material) :: mat
	    type(XSDlay)   :: dlayxs
	    type(Option)   :: opt

          integer(NBI)   :: i, niso, ng, sg, g, nvar, t, ig, k
          real(NBF)      :: sum_chi, fact, sum_2ns, max_kf
          real(NBF)      :: temp, sum_mac_chi, sum_mac_xs_nf, sum_iso_xs_nf

          niso = mat%no_isotope
          ng   = mat%no_group

          do t=1, TCURR
             mat%kappa(t) = 0.0
          enddo
          do g=1,ng
             do t=1, TCURR
                mat%xs_f(t,g)  = 0.0
                mat%xs_nf(t,g) = 0.0
                mat%xs_kf(t,g) = 0.0
                mat%xs_a(t,g)  = 0.0
                mat%xs_c(t,g)  = 0.0
                mat%xs_r(t,g)  = 0.0
                mat%tnu(t,g)   = 0.0
                mat%chip(t,g)  = 0.0
                mat%xs_tr(t,g) = 0.0
             enddo
             do sg=1,ng
                do t=1, TCURR
                   mat%xs_s(t,sg,g) = 0.0
                enddo
             enddo
          enddo

          do i=1,niso
             nvar = mat%no_temp_var(i)
             do g=1,ng
                mat%xs_f(TCURR,g)  = mat%xs_f(TCURR,g)  +  mat%isotope(TBASE,i)%xs_nf(g)*mat%iso_density(i)
                mat%xs_nf(TCURR,g) = mat%xs_nf(TCURR,g) &
               +  mat%isotope(TBASE,i)%xs_nf(g)*mat%isotope(TBASE,i)%tnu(g)*mat%iso_density(i)
                mat%xs_kf(TCURR,g) = mat%xs_kf(TCURR,g) &
               +  mat%isotope(TBASE,i)%xs_nf(g)*mat%isotope(TBASE,i)%kappa*mat%iso_density(i)
                mat%xs_tr(TCURR,g) = mat%xs_tr(TCURR,g) +  mat%isotope(TBASE,i)%pl_xs_tr(g,1)*mat%iso_density(i)
                mat%xs_c(TCURR,g)  = mat%xs_c(TCURR,g)  + (mat%isotope(TBASE,i)%xs_na(g) + mat%isotope(TBASE,i)%xs_np(g)  &
               + mat%isotope(TBASE,i)%xs_nt(g)  &
               +  mat%isotope(TBASE,i)%xs_nd(g) + mat%isotope(TBASE,i)%xs_ngm(g) )*mat%iso_density(i)
                mat%xs_a(TCURR,g)  = mat%xs_c(TCURR,g)  + mat%xs_f(TCURR,g)

                do sg=1,ng
                   mat%xs_s(TCURR,sg,g) = mat%xs_s(TCURR,sg,g) + mat%isotope(TBASE,i)%xs_scat(sg,g)*mat%iso_density(i)
                enddo
             enddo
          enddo
  ! chi
          do g=1,ng
             sum_mac_chi   = 0.0
             sum_mac_xs_nf = 0.0
             do i=1,niso
                sum_iso_xs_nf = 0.0
                do ig=1,ng
                   sum_iso_xs_nf = sum_iso_xs_nf + mat%isotope(TBASE,i)%xs_nf(ig)*mat%isotope(TBASE,i)%tnu(ig)
                enddo
                if( mat%isotope(TBASE,i)%flag_chi .eq. FLAG_ON ) then
                   sum_mac_chi   = sum_mac_chi + mat%isotope(TBASE,i)%chiv(g)*mat%iso_density(i)*sum_iso_xs_nf
                   sum_mac_xs_nf = sum_mac_xs_nf + mat%iso_density(i)*sum_iso_xs_nf
  !               else if( mat%isotope(TBASE,i)%flag_chi .gt. FLAG_ON ) then
  !                  do k=1,FLAG_ON
  !                     mat%chi(TBASE,g) = mat%chi(TBASE,g)                                         ! isotope chi matrix
  !                  enddo
                else if( mat%isotope(TBASE,i)%flag_chi .eq. FLAG_OFF ) then
                   if( const%flag_chi .eq. FLAG_ON ) then
                      sum_mac_chi = sum_mac_chi + const%wide_chiv(g)*mat%iso_density(i)*sum_iso_xs_nf
                      sum_mac_xs_nf = sum_mac_xs_nf + mat%iso_density(i)*sum_iso_xs_nf
                   else if( const%flag_chi .gt. FLAG_ON ) then
                      do k=1,FLAG_ON
                         sum_mac_chi = sum_mac_chi + const%wide_chim(k,g)*mat%iso_density(i)*sum_iso_xs_nf  ! wide chi matrix
                      enddo
                      sum_mac_xs_nf = sum_mac_xs_nf + mat%iso_density(i)*sum_iso_xs_nf
                   endif
                endif
             enddo
             if( sum_mac_xs_nf .lt. NEAR_ZERO ) then
                mat%chip(TBASE,g) = 0.0
             else
                mat%chip(TBASE,g) = sum_mac_chi/sum_mac_xs_nf
             endif
          enddo

          sum_chi = 0.0
          do g=1,ng
             sum_chi = sum_chi + mat%chip(TBASE,g)
          enddo
          do g=1,ng
             if( sum_chi .lt. NEAR_ZERO )  sum_chi = NEAR_INF
             mat%chip(TBASE,g) = mat%chip(TBASE,g)/sum_chi
             mat%chip(TCURR,g) = mat%chip(TBASE,g)
             mat%chit(TCURR,g) = mat%chip(TBASE,g)
          enddo

  ! correction sigma_kappa_f
          max_kf = 0.0
          do g=1,ng
  !         do t=1, nvar
             if( mat%xs_kf(TCURR,g) .gt. max_kf ) max_kf = mat%xs_kf(TCURR,g)
          enddo

          if( max_kf .lt. NEAR_ZERO ) then
             do g=1,ng
                mat%xs_kf(TCURR,g) = mat%xs_f(TCURR,g)*CONV
             enddo
          endif

       end subroutine generate_macxs_mat

  

	subroutine generate_macxs_comp( comp, ifinit )

	type(Composition) :: comp  
	integer(NBI) :: l, nmat, ng, sg, g, ig
	real(NBF) :: sum_chip, sum_chit, sum_chid, sum_nf
	logical :: ifinit

	nmat = comp%no_material
	ng   = comp%no_group

	call reset_macxs_comp( comp )

	do g=1,ng
	do l=1,nmat
! modification : adding the sum of xs_nf
		sum_nf = 0.0
		do ig = 1, ng
			sum_nf = sum_nf + comp%material(l)%xs_nf(TCURR,ig)
		enddo

		comp%xs_f(g)  = comp%xs_f(g)  + comp%material(l)%xs_f(TCURR,g)*comp%vol_fraction(l)   ! 2014_10_24 . scb
		comp%xs_nf(g) = comp%xs_nf(g) + comp%material(l)%xs_nf(TCURR,g)*comp%vol_fraction(l)
		comp%xs_kf(g) = comp%xs_kf(g) + comp%material(l)%xs_kf(TCURR,g)*comp%vol_fraction(l)

!		comp%xs_c(g)  = comp%xs_c(g)  + comp%material(l)%xs_c(TCURR,g)*comp%vol_fraction(l)
		comp%xs_a(g)  = comp%xs_a(g)  + comp%material(l)%xs_a(TCURR,g)*comp%vol_fraction(l)
!		comp%tnu(g)   = comp%tnu(g)   + comp%material(l)%tnu(TCURR,g)*comp%vol_fraction(l)
		comp%xs_tr(g) = comp%xs_tr(g) + comp%material(l)%xs_tr(TCURR,g)*comp%vol_fraction(l)

		comp%chip(g)  = comp%chip(g)  + comp%material(l)%chip(TCURR,g)*comp%vol_fraction(l)*sum_nf
		comp%chit(g)  = comp%chit(g)  + comp%material(l)%chit(TCURR,g)*comp%vol_fraction(l)*sum_nf
		comp%chid(g)  = comp%chid(g)  + comp%material(l)%chid(TCURR,g)*comp%vol_fraction(l)*sum_nf

		do sg=1,ng
			comp%xs_s(sg,g) = comp%xs_s(sg,g) + comp%material(l)%xs_s(TCURR,sg,g)*comp%vol_fraction(l)
		enddo
	enddo
	if(ifinit) then
		comp%xs_f0(g)=comp%xs_f(g)
		comp%xs_nf0(g)=comp%xs_nf(g)
		comp%xs_kf0(g)=comp%xs_kf(g)
		comp%xs_c0(g)=comp%xs_c(g)
		comp%xs_a0(g)=comp%xs_a(g)
		comp%xs_tr0(g)=comp%xs_tr(g)

		comp%chip0(g)=comp%chip(g)
		comp%chit0(g)=comp%chit(g) 
		comp%chid0(g)=comp%chid(g)
	endif
	enddo


	do g=1,NPRECT
	do l=1,nmat
		comp%beta(g)  = comp%beta(g)  + comp%material(l)%beta(g)
		comp%dct(g)  = comp%dct(g)  + comp%material(l)%dct(g)
	enddo
	enddo

!#ifdef NOMALIZE
	sum_chip = 0.0
	sum_chit = 0.0
	sum_chid = 0.0
	do g=1,ng
		sum_chip = sum_chip + comp%chip(g)
		sum_chit = sum_chit + comp%chit(g)
		sum_chid = sum_chid + comp%chid(g)
	enddo
	do g=1,ng
		if( sum_chip .lt. NEAR_ZERO )  sum_chip = NEAR_INF
		comp%chip(g) = comp%chip(g)/sum_chip

		if( sum_chit .lt. NEAR_ZERO )  sum_chit = NEAR_INF
		comp%chit(g) = comp%chit(g)/sum_chit

		if( sum_chid .lt. NEAR_ZERO )  sum_chid = NEAR_INF
		comp%chid(g) = comp%chid(g)/sum_chid
	enddo
!#endif

!     call generate_macxs_comp_uncards(comp)
	return
	end subroutine generate_macxs_comp


       subroutine generate_macxs_comp_uncards(comp)
          type(Composition) :: comp 
          integer(NBI) :: l, nmat, ng, sg, g, tg, t, nt
          real(NBF) :: sum_xs_s

          nmat = comp%no_material
          ng   = comp%no_group

          do l=1,nmat
             nt = comp%material(l)%max_temp_var
             do g=1,ng
                if( nt .eq. NTHVAR ) then
   	               do t=1, NTHVAR_UNC
                      comp%xsu_f(g,t)  = comp%xsu_f(g,t)  + comp%material(l)%xs_f(t,g)*comp%vol_fraction(l)
                      comp%xsu_nf(g,t) = comp%xsu_nf(g,t) + comp%material(l)%xs_nf(t,g)*comp%vol_fraction(l)
                      comp%xsu_kf(g,t) = comp%xsu_kf(g,t) + comp%material(l)%xs_kf(t,g)*comp%vol_fraction(l)

                      comp%xsu_a(g,t)  = comp%xsu_a(g,t)  + comp%material(l)%xs_a(t,g)*comp%vol_fraction(l)
                      comp%xsu_tr(g,t) = comp%xsu_tr(g,t) + comp%material(l)%xs_tr(t,g)*comp%vol_fraction(l)
                      do sg=1,ng
                         comp%xsu_s(sg,g,t) = comp%xsu_s(sg,g,t) + comp%material(l)%xs_s(t,sg,g)*comp%vol_fraction(l)
                      enddo
			       enddo
                else
   	               do t=1, NTHVAR_UNC
                      comp%xsu_f(g,t)  = comp%xsu_f(g,t)  + comp%material(l)%xs_f(TBASE,g)*comp%vol_fraction(l)
                      comp%xsu_nf(g,t) = comp%xsu_nf(g,t) + comp%material(l)%xs_nf(TBASE,g)*comp%vol_fraction(l)
                      comp%xsu_kf(g,t) = comp%xsu_kf(g,t) + comp%material(l)%xs_kf(TBASE,g)*comp%vol_fraction(l)

                      comp%xsu_a(g,t)  = comp%xsu_a(g,t)  + comp%material(l)%xs_a(TBASE,g)*comp%vol_fraction(l)
                      comp%xsu_tr(g,t) = comp%xsu_tr(g,t) + comp%material(l)%xs_tr(TBASE,g)*comp%vol_fraction(l)
                      do sg=1,ng
                         comp%xsu_s(sg,g,t) = comp%xsu_s(sg,g,t) + comp%material(l)%xs_s(TBASE,sg,g)*comp%vol_fraction(l)
                      enddo
                   enddo
                endif
             enddo
          enddo

          do g=1,ng
		     do t=1, NTHVAR_UNC
                sum_xs_s = 0.0
                do tg=1,ng
			       if( tg .ne. g ) then
                      sum_xs_s = sum_xs_s + comp%xsu_s(g,tg,t)
			       endif
                enddo
                comp%xsu_r(g,t)  = comp%xsu_a(g,t) + sum_xs_s
             enddo
          enddo

       end subroutine generate_macxs_comp_uncards

	  end module MaterialMod



  !      DO J=1,NTFMAT
  !         DO I=1,NGROUP
  !            BETJ = 0.0
  !	         DO K=1,NFAM
  !               BETJ = BETJ + DNU(K,J,I)
  !	         END DO
  !            TNU  = PNU(J,I) + BETJ
  !            if(tnu.eq.0) tnu=1.e10 !temp hj02jul05

  !	         BETJ = BETJ / TNU
  !            CHIT(I,J) = (1.D+00 - BETJ)*CHIP(J,I)

  !	         DO K=1,NFAM
  !               BETI = DNU(K,J,I) / TNU
  !               CHIT(I,J) = CHIT(I,J) + CHID(I,K)*BETI
  !	         END DO
  !         END DO
  !      END DO



  !      do i=1,niso
  !         do g=1,ng
  !            if( mat%isotope(TBASE,i)%flag_chi .eq. FLAG_ON ) then
  !               mat%chip(TCURR,g) = mat%chip(TCURR,g) + mat%isotope(TBASE,i)%chiv(g)*mat%iso_density(i)
  !!               else if( mat%isotope(TBASE,i)%flag_chi .gt. FLAG_ON ) then
  !!                  do k=1,FLAG_ON
  !!                     mat%chip(TCURR,g) = mat%chip(TCURR,g) + mat%isotope(TBASE,i)%chim(k,g)*mat%iso_density(i)              ! isotope chi matrix
  !!                  enddo
  !            else if( mat%isotope(TBASE,i)%flag_chi .eq. FLAG_OFF ) then
  !               if( const%flag_chi .eq. FLAG_ON ) then
  !                  mat%chip(TCURR,g) = mat%chip(TCURR,g) + const%wide_chiv(g)*mat%iso_density(i)
  !!                  else if( const%flag_chi .gt. FLAG_ON ) then
  !!                     do k=1,FLAG_ON
  !!                        mat%chip(TCURR,g) = mat%chip(TCURR,g) + const%%wide_chim(k,g)*mat%iso_density(i)              ! wide chi matrix
  !!                     enddo
  !               endif
  !            endif
  !         enddo
  !      enddo
