! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
    subroutine alloctpen

      use allocs
      use param

! Allocate Memory to Dynamic Varible

      include 'global.h'
      include 'cntl.h'
      include 'files.h'
      include 'itrcntl.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'xsec.h'
      include 'defhex.h'
      include 'defsfc.h'
      include 'defpnt.h'
      include 'deffg.h'
      include 'lscoefh.h'
      include 'dummy.h'
      include 'lufac.h'
!
      ncomp=ncomp    !NCOMP
!
      ntph=6*ndivhs*ndivhs
      nasst2=ncomp  !HEXAGON_PR
      nasst3=nat    !ASSY_TYPE
      ncxt=ncomp    !NCOMP
!
! Memory Allocation for Solution Vector in "hsolvec" of geomhfc.h
!
      nhexnode=nassy*nz
      ntrinode=ntph*nhexnode
      nhnodv2=ng*nhexnode
      ntnodv2=ng*ntrinode
      nzhnodv2=ng*2*nhexnode
      nhnodvm=ng*nhexnode
      ntnodvm=ng*ntrinode
      nzhnodvm=ng*2*nhexnode
      ngrid=(2*nxfc+1)*(2*nyfc+1)
      nradinfo=6*nassy
      nzinfo=2*nz
      ncxt=ng*ncomp
      npntv=ng*ncorn*nz
      nsfcv=ng*nsurf*nz
      nzsfcv=ng*nassy*(nz+1)
      ngridy=(2*nyfc+1)

! add
      call dmalloc(igc,ng)	
      call dmalloc(ifbcref,ng)
      call dmalloc(ifbcrefzb,ng)
      call dmalloc(ifbcrefzt,ng)

      call dmalloc0(lfaptr,0,nxy)
      call dmalloc(lfatol,nxy)
! lufac.h
      call dmalloc(del,4,nx)
      call dmalloc(delinv,4,nxy,nz)
      call dmalloc(al,4,nxy,nz)
      call dmalloc(au,4,nx)
      call dmalloc(deliau,4,nxy,nz)
      call dmalloc(ainvl,4,nx)
      call dmalloc(ainvu,4,nx)
      call dmalloc(ainvd,4,nx)
! end
!
      call dmalloc(cnto,ng,ntph,nassy,nz)
      call dmalloc(cntzo,ng,2,nassy,nz)
      call dmalloc(atleak,ng,nassy,nz)
      call dmalloc(srcz,ng,nassy,nz)
!
! Memory Allocation for integer Information in "ihnodinf" of geomhfc.h
!
      call dmalloc(icxt,nassy,nz)
      call dmalloc0(iaass,-nxfc,nxfc,-nyfc,nyfc)
      call dmalloc0(iapoint,-nxfc,nxfc,-nyfc,nyfc)
      call dmalloc0(iasurf,-nxfc,nxfc,-nyfc,nyfc)
      call dmalloc(neignd,6,nassy)
      call dmalloc(neigz,2,nz)
      call dmalloc(neigpt,6,nassy)
      call dmalloc(neigjin,6,nassy)
      call dmalloc(neigsfc,6,nassy)
      call dmalloc(neigsfcz,2,nz)
! Memory Allocation for point information of defpnt.h
      call dmalloc(pflx,ng,ncorn,nz)
      call dmalloc(pbdv,ng,ncorn)
      call dmalloc(pcpbc,ng,ncorn)
      call dmalloc(pcpbsbd,ng,ncorn)
      call dmalloc(cmatpnt,12,ncorn)
      call dmalloc(pflxt,ncorn)
      call dmalloc(codpnt,ncorn)
      call dmalloc(neigppt,12,ncorn)
      call dmalloc(neignpt,ncorn)
      call dmalloc(imatid,5,6,nassy)
! Memory Allocation for coefficients for CMFD solver of defsfc.h
      call dmalloc(dhat,ng,nsurf,nz)
      call dmalloc(dhatz,ng,nassy,nz+1)
      call dmalloc(dfd,ng,nsurf,nz)
      call dmalloc(dfdz,ng,nassy,nz+1)
      call dmalloc(betaphis,ng,nsurf,nz)
      call dmalloc(betaphisz,ng,nassy,nz+1)
!dj+ : add for under-relaxation & reactivity edit
      call dmalloc(dhat0,ng,nsurf,nz)
      call dmalloc(dhatz0,ng,nassy,nz+1)
      call dmalloc(dhatd,ng,nsurf,nz)
      call dmalloc(dhatzd,ng,nassy,nz+1)
      call dmalloc(betaphisd,ng,nsurf,nz)
      call dmalloc(betaphiszd,ng,nassy,nz+1)
!dj-
      call dmalloc(wtdhat,ntph,nassy)
      call dmalloc(dsum,ng,nassy,nz)
! Memory Allocation for neighbor information for each surface of defsfc.h
      call dmalloc(neigsnd,5,nsurf)
      call dmalloc(neigsndz,5,nz+1)
! Memory Allocation for core geometry inforation of geomhfc.h
      call dmalloc0(icorexs,-nyfc,nyfc)
      call dmalloc0(icorexe,-nyfc,nyfc)
      call dmalloc0(icorexfs,-nyfc,nyfc)
      call dmalloc0(icorexfe,-nyfc,nyfc)
      call dmalloc0(icorexps,-nyfc,nyfc)
      call dmalloc0(icorexpe,-nyfc,nyfc)
      call dmalloc0(icorexss,-nyfc,nyfc)
      call dmalloc0(icorexse,-nyfc,nyfc)
! lscoefh.h
      !call dmalloc(hflx,ng,nassy,nz)
      call dmalloc(hflx,ng2,nassy,nz)  ! 2013_05_20 . scb
      call dmalloc(fsrc,nassy,nz)
      call dmalloc(fswiel,nassy,nz)
      call dmalloc(cmat,ng2,8,nassy,nz)
      call dmalloc(dcmat,ng2*ng2,nassy,nz)

      call dmalloc(cmatd,ng,8,nassy,nz)
      !call dmalloc(dcmatd,ng*ng,nassy,nz)     ! 2013_09_11 . scb

      call dmalloc(xlufac,ng,6,nassy,nz)
      call dmalloc0(ipntr,0,6,1,nassy)
      call dmalloc0(ineigcond,0,6,0,6,1,nassy)
      call dmalloc(ilubnd,nassy)
      call dmalloc(iastopnt,nassy,nassy)

! Memory Allocation for Dummy Variable
      call dmalloc(dum1,ng)
      call dmalloc(dum2,ng)
      call dmalloc(dum3,ng)
      call dmalloc(dum4,ng)
      call dmalloc(cftm,ng*ng,7,nassy)
      call dmalloc(dumrv,ng,nassy)
      call dmalloc(dumrs,ng,nassy)
      call dmalloc0(phidum,1,ng,1,nassy,0,nz)
!
      ntrinode=6*nhexnode
      nhnod2gf=ng*nhexnode
      call dmalloc(aflx,ng,ntph,nassy,nz)
      if(ng.ne.2) then
         call dmalloc(xmom,ng,ntph,nassy,nz)
         call dmalloc(ymom,ng,ntph,nassy,nz)
         call dmalloc(hflxf,ng,nassy,nz)
         call dmalloc(phiadjf,ng,nassy,nz)
         call dmalloc(zmom1,ng,nassy,nz)
         call dmalloc(zmom2,ng,nassy,nz)
         call dmalloc(fhflx,ng,nassy,nz)
         call dmalloc(fohflx,ng,nassy,nz)
         call dmalloc(fcnto,ng,6,nassy,nz)
         call dmalloc(fcntzo,ng,2,nassy,nz)
         call dmalloc(focnto,ng,6,nassy,nz)
         call dmalloc(focntzo,ng,2,nassy,nz)
      else
         hflxf=>hflx
      endif

      return
    end subroutine
