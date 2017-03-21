! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
module geomhex

	use const
	implicit none

! input arguement
	integer               :: ng                   ! number of group
	integer               :: nz, nassy            ! number of nodes in x, y and z direction and total nodes
	integer               :: nxpnt
	integer               :: ncorn
	integer               :: nxsfc
	integer               :: nsurf
	integer               :: nxy

	integer :: kfbeg,kfend
	integer :: ndivhs
	real :: hf2f, hexarea
 ! hex_fdm
  integer :: ndiv, nsub
  real :: hdiv
  logical :: ifhexfdm, ifhexsp3

	integer,pointer,dimension(:,:)  :: neigjin, neignd, neigpt, neigz,  &
	                                   neigsfc, neigsnd, neigsndz

	real,pointer,dimension(:) :: wtass, codpnt
  real,pointer,dimension(:,:) :: pbdv, wtdhat

   ! albedo
	real :: alxr, alzr, alzl
  real :: reflrat, reflratzb, reflratzt
  logical :: ifbcref, ifbcrefzb, ifbcrefzt

	! related to node size 
	real,pointer,dimension(:,:)     :: volnode
	real,pointer,dimension(:)       :: hz

	integer,pointer,dimension(:,:)  :: neigsfcz, ipntr, iastopnt
	integer,pointer,dimension(:,:,:) :: ineigcond
	integer,pointer,dimension(:)    :: ilubnd

  integer,pointer :: neigtri(:,:,:), neigtria(:,:,:)
  real :: triarea
  real,pointer :: volnodet(:,:,:)
! 
! processing param.
	real :: rt3, rsqrt3, hside, twort3o9h
	integer :: ntph

	contains

	subroutine setgeomhex( ngl,nzl,nxyl,nassyl,nxpntl,ncornl,nxsfcl,    &
                           nsurfl,                                 &
	                        hf2fl,ndivhsl,                          &
						         kfbegl,kfendl,                          &
						         volnodel,hzl,                        &
						         alxrfl,alzlfl,alzrfl,                   &
						         wtassl,neigndl,neigjinl,neigptl,neigzl, &
						         neigsfcl,neigsndl,neigsndzl,            &
						         codpntl,                                &
						         pbdvl,                                  &
						         wtdhatl,                                &
						         neigsfczl,                              &
						         ineigcondl,                             &
						         ipntrl,                                 &
						         ilubndl,                                &
						         iastopntl,                              &
                           ifhexfdml,ifhexsp3l,                    &
                           ndivl,neigtril,neigtrial                &
						        )
										

		use allocs

		integer               :: ngl                   ! number of group
		integer               :: nzl, nassyl            ! number of nodes in x, y and z direction and total nodes
		integer               :: nxpntl
		integer               :: ncornl
		integer               :: nxsfcl
		integer               :: nsurfl
		integer               :: nxyl
		integer               :: ndivl
    logical               :: ifhexfdml, ifhexsp3l



		integer :: kfbegl,kfendl
		integer :: ndivhsl
		real :: hf2fl
		integer,pointer,dimension(:,:)  :: neigjinl, neigndl, neigptl, neigzl,  &
		                                   neigsfcl, neigsndl, neigsndzl

		real,pointer,dimension(:) :: wtassl, codpntl
		real,pointer,dimension(:,:) :: pbdvl, wtdhatl

		! albedo
		real :: alxrfl, alzrfl, alzlfl

		! related to node size 
		real,pointer,dimension(:,:)     :: volnodel
		real,pointer,dimension(:)       :: hzl

		integer,pointer,dimension(:,:)  :: neigsfczl, ipntrl, iastopntl
		integer,pointer,dimension(:,:,:) :: ineigcondl
		integer,pointer,dimension(:)    :: ilubndl

    integer,pointer :: neigtril(:,:,:), neigtrial(:,:,:)

    integer :: k, l, n
    real :: hzn

		ng=ngl; nz=nzl; nxy=nxyl
		nxpnt=nxpntl; ncorn=ncornl;
		ndivhs=ndivhsl; nassy=nassyl;
		hf2f=hf2fl; nxsfc=nxsfcl;
		kfbeg=kfbegl; kfend=kfendl;
		nsurf=nsurfl
    ndiv=ndivl 
    ifhexfdm=ifhexfdml 
    ifhexsp3=ifhexsp3l

		alxr=alxrfl
		alzr=alzrfl
		alzl=alzlfl

		neigjin =>  neigjinl
		neignd  =>  neigndl
		neigpt  =>  neigptl
		neigz   =>  neigzl
		wtass   =>  wtassl
		codpnt  =>  codpntl
		pbdv    =>  pbdvl
		neigsfc =>  neigsfcl
		neigsnd =>  neigsndl
		neigsndz=>  neigsndzl
		wtdhat  =>  wtdhatl

		volnode =>  volnodel
		hz      =>  hzl
    neigtri =>  neigtril
    neigtria=>  neigtrial

		neigsfcz  => neigsfczl
		ineigcond => ineigcondl
		ipntr     => ipntrl
		ilubnd    => ilubndl
		iastopnt  => iastopntl


		reflrat=(1-2*alxr)/(1+2*alxr)
		reflratzb=(1-2*alzl)/(1+2*alzl)
		reflratzt=(1-2*alzr)/(1+2*alzr)

    ifbcref=FALSE
    ifbcrefzb=FALSE
    ifbcrefzt=FALSE
    if(alxr.eq.0) ifbcref=TRUE
    if(alzl.eq.0) ifbcrefzb=TRUE
    if(alzr.eq.0) ifbcrefzt=TRUE

      
		rt3=sqrt(3.)
		rsqrt3=1./rt3
		hside=hf2f*rsqrt3

		ntph=6*ndivhs*ndivhs
		hexarea=1.5d0*hf2f*hside

		twort3o9h=2*rt3/9./hside

! hex_fdm
    if(ifhexfdm) then
        nsub=6*ndiv*ndiv
        hdiv=hside/ndiv
        triarea=0.25*rt3*hdiv*hdiv
        call dmalloc(volnodet,nsub,nassy,nz)
        do k=1,nz
          if(nz.eq.1) then
              hzn=hz(k)
          else
              hzn=hz(k)
          endif
          do l=1,nassy
              do n=1,nsub             
                volnodet(n,l,k)=triarea*hzn
              enddo
          enddo
        enddo
    endif

	end subroutine
    
end module
