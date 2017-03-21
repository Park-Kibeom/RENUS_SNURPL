! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB  
! program for generating multigroup cross-section from the ISOTXS format file.
! 2007. 5 
!
	subroutine trinx0

	use param
	use BasicParam
	use XSConstMod
	use IsotopeMod
	use MaterialMod
      use XSDlayMod
	use OptionMod
	use trinx_cntl

    include 'global.h'
	include 'thexpan.inc'
	include 'thlink.inc'
	include 'files.h'
	include 'xsec.h'

	integer :: icomp
	integer :: atype
	integer :: l
	character*6 :: compname 
	real(NBF) :: rodfrac, delt
	logical :: ifinit
	logical, save :: first=.TRUE.

	type(Option)  :: opt
	type(XSConst) :: const
	type(XSDlay)  :: dlayxs
	type(Isotope),     pointer, dimension (:) :: isotxs
	type(Material),    pointer, dimension (:) :: mat
	type(Composition), pointer, dimension (:) :: comp
      
      common / scbtrinx / opt, const, dlayxs, isotxs, mat, comp
      ! 2013_04_01 . scb

	call read_input(opt,filemat)
	if(opt%iso_flag .eq. FLAG_ON ) call read_isotxs(isotxs,const,fileiso)		
	if(opt%dly_flag .ne. FLAG_OFF ) call read_xsdlay(dlayxs,opt%dly_flag,filedla)		
 
	call read_material(mat,comp,const,filemat)    
	call make_material(mat,isotxs)
	call make_composition(comp,mat,isotxs,opt)
	
	ncomp = size(comp)
	do g=1,ng
		velof(g,:)=const%velocity(g)
	enddo

	return


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          ENTRY          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	entry trinx(ifinit,icomp,compname,atype,rodfrac)

	nmat = size(mat)
	ncomp = size(comp)


	if(ifreact) then
		if(first) then
			opt%fuel_temp = 700.
			opt%cool_temp = 660.
			first=.FALSE.
		else
			opt%fuel_temp = base_ft + del_dop
			opt%cool_temp = base_ct
		endif
	else
		opt%fuel_temp = tfueln
		opt%cool_temp = tcooln
	endif

	do l=1,ncomp
		if(compname.eq.comp(l)%name) then
			if(.not.ifinit) call updmat(comp(l),atype,rodfrac) 
			call make_macxs(comp(l),const,dlayxs,opt,ifinit)
			call linkxs(icomp,compname,comp(l))
		endif
	enddo
        

!	call write_xs(comp,const,opt,mxc,compname)
!     call delete_isotope(isotxs,const)
!	close(nerror)

	return

	end subroutine 
