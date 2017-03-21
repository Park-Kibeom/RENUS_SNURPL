! added in ARTOS ver. 0.2 ( for Thermal Expansion ). 2012_07_04 by SCB	
    subroutine updgeom
! update the geometry to consider the thermal expansion effect

      use param
      use allocs
      use trinx_cntl
!	use geom, only : hmesh0 => hmesh, volnode0 => volnode, volfuel0 => volfuel

      include 'global.h'
      include 'geom.h'
      include 'xsec.h'
      include 'ffdm.h'
      include 'pinpr.h'
      include 'thcntl.inc'
      include 'thgeom.inc'
      include 'thexpan.inc'

      logical,save :: first=TRUE

! update geometries
      do ia=1,nxa
        hxa(ia)=hxa0(ia)*ar
      enddo

      do ja=1,nya
        hya(ja)=hya0(ja)*ar
      enddo

      do ka=kfbega,kfenda
        hza(ka)=hza0(ka)*azf
      enddo

      if(ifreact) then
        do irod=1,nrodtyp
          rodstepsize(irod)=rodstepsize0(irod)*azc
        enddo
      else
        do irod=1,nrodtyp
          if(isteptr.eq.0) rodstepsize(irod)=rodstepsize0(irod)*azf
        enddo
      endif

! mesh sizes
      do k=1,nz
        ka=ktoka(k)
        do l=1,nxy
          la=ltola(l)
          ia=latoia(la)
          ja=latoja(la)
          hmesh(XDIR,l,k)=hxa(ia)/nmeshx(ia)
          hmesh(YDIR,l,k)=hya(ja)/nmeshy(ja)
          hmesh(ZDIR,l,k)=hza(ka)/nmeshz(ka)

          rhmesh(XDIR,l,k)=1/hmesh(XDIR,l,k)
          rhmesh(YDIR,l,k)=1/hmesh(YDIR,l,k)
          rhmesh(ZDIR,l,k)=1/hmesh(ZDIR,l,k)
        enddo
      enddo
 
      do k=1,nz
        do l=1,nxy
          volnode(l,k)=hmesh(XDIR,l,k)*hmesh(YDIR,l,k)*hmesh(ZDIR,l,k)
        enddo
      enddo


      coreheight=0
      do ka=kfbega,kfenda
         coreheight=coreheight+hza(ka)
      enddo

! volume
      volcore=0
      volfuel=0
      do ka=1,nza
        do la=1,nxya
          ia=latoia(la)
          ja=latoja(la)
          iat=iassytyp(la)
          volnodea(la,ka)=hxa(ia)*hya(ja)*hza(ka)
          volcore=volcore+volnodea(la,ka)
          if(ka.eq.1 .and. iffuela(iat)) volfuel=volfuel+coreheight*hxa(ia)*hya(ja)
        enddo
      enddo
      volfuel0=volfuel

! update geometries for th-calculation
#ifdef FDBK
      if(fdbk) then
        k=kfbeg
        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          if(.not.iffuela(iat)) cycle
          lchan=ltochan(l)
          chanvol(lchan)=chanvol(lchan)+volnode(l,k)
        enddo
      
        do lchan=1,nchan
          chanvol(lchan)=chanvol(lchan)*rhmesh(ZDIR,1,k)
        enddo

        nzthp1=nzth+1                           !number of junctions
        do kth=1,nzth
          hzth(kth)=0
          do k=junb(kth-1)+1,junb(kth)
            ka=ktoka(k)
            hzth(kth)=hzth(kth)+hmesh(ZDIR,1,k)
            ktokth(k)=kth
          enddo
          hzth(kth)=hzth(kth)*0.01             !cm to m
        enddo
      endif
#endif

      return
    end subroutine
