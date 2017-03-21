! file name = dimpar.h - array dimension parameters 
! dimension parameters
      logical rect,hex,symedge
      logical,pointer,dimension(:) :: iffuela,iffuella
      integer,pointer,dimension(:) :: iassytyp,  &
        nxsa,nxea,nrowxa,nxs,nxe,nrowx,latoia,latoja,ltoi,ltoj,ltola,  &
        itoia,jtoja,ktoka,nysa,nyea,nrowya,nys,nye,nrowy,  &
        nmeshx,nmeshy,nmeshz,kfmb,kfme,  &
        nodema,nodepa,jsfcda,nodem,nodep,jsfcd,  &
        nxsfa,nxefa,nxsf,nxef

      integer,pointer,dimension(:,:) :: nodela,nodel,  &    !(i,j)
                                        neibz,  & !(2,nza)
                                        neibr,neibra,  & !(nrdir2,nxy)
                                        lsfca,lsfc       !(nrdir2,nxy)
      integer,pointer :: icompnumz(:,:),icomprod(:)  !(k,iat),(nrodtype)
      character(4) symopt
      real,pointer,dimension(:) :: hxa,hya,hza,heleva,helev
      real,pointer,dimension(:) :: albr,albrl,albrr,albzb,albzt,albxl,albxr,albyl,albyr
      real,pointer,dimension(:,:) :: volnodea,volnode
!
!      type(cmtofm),pointer,dimension(:) :: latol  ! 2013_05_14 . scb
!

! mesh size
      real, pointer, dimension(:,:,:) :: hmesh,rhmesh

      common /geoml/rect,hex,iffuela,iffuella
      common /geoma/symopt
      common /geomi/isymang,isymloc,kfbeg,kfend,kfbega,kfenda,  &
          iassytyp,icompnumz,  &
          nxsa,nxea,nrowxa,nxs,nxe,nrowx,latoia,latoja,ltoi,ltoj,ltola,  &
          itoia,jtoja,ktoka,nysa,nyea,nrowya,nys,nye,nrowy,  &
          nmeshx,nmeshy,nmeshz,kfmb,kfme,  &
          nodema,nodepa,jsfcda,nodem,nodep,jsfcd,  &
          nodela,nodel,  &
          neibz,  &
          neibr,neibra,  &
          lsfca,lsfc,nsurfx,nsurfxa,  &
          nxsfa,nxefa,nxsf,nxef,   &
		  jfbeg,jfend         ! added in ARTOS ver. 0.2 ( for zeroleak option ). 2012_07_03 by SCB

      common /geoms/domheight,coreheight,volcore,volfuel
      common /geomf/hmesh,rhmesh,hxa,hya,hza,heleva,helev,  &
                 volnodea,volnode,  &
                 albr,albrl,albrr,albzb,albzt,albxl,albxr,albyl,albyr

!      common /ffdmdt/latol   ! 2013_05_14 . scb


!
! r/b map
      integer,pointer,dimension(:,:) :: maprb
      common /geomrb/ maprb
!
! corner points
      common/geomcorn/ ncorn


! rod configuration
      integer :: nrodtyp
      integer,pointer,dimension(:) :: irodtyp
      real,pointer,dimension(:) :: rodfullpos,rodstep,rodstepsize,rodfrac(:,:), &
	                               rodstepsize0   ! added in ARTOS ver. 0.2 ( for thermal expansion ). 2012_07_03 by SCB

      real,pointer,dimension(:) :: rodstep0   ! 2014_12_22 . scb

! 2012_08_23 . scb
 	  real,pointer,dimension(:) :: pscrmbeg, rodstepd
	  common / geomrodfa2 / pscrmbeg, rodstepd
! added end
      
      common /geomrodi/nrodtyp
      common /geomrodia/irodtyp
      common /geomrodfa/rodfullpos,rodstep,rodstepsize,rodfrac, &
	                    rodstepsize0   ! added in ARTOS ver. 0.2 ( for thermal expansion ). 2012_07_03 by SCB
      common /geomrodfa2/rodstep0   ! 2014_12_22 . scb

      integer :: ibcxl,ibcxr,ibcyl,ibcyr,ibczl,ibczr
      common /geombci/ ibcxl,ibcxr,ibcyl,ibcyr,ibczl,ibczr

      real,pointer,dimension(:) :: srstep, srfrac(:,:)  ! added in ARTOS ver. 0.2 ( for SCRAM ). 2012_07_03 by SCB
      common /geomscram/ srstep,srfrac                  ! added in ARTOS ver. 0.2 ( for SCRAM ). 2012_07_03 by SCB

! 2013_05_01 . scb
      logical :: ifrecsp3
	  common /geomlogsp3 / ifrecsp3 
! added end