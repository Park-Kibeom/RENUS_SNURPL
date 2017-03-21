    subroutine initgeom
!
      use param
      use allocs
!
      include 'global.h'
      include 'geom.h'
      include 'xsec.h'
      include 'ffdm.h'
      include 'pinpr.h'
! added in ARTOS ver. 0.2 (Stuck rod, Thermal expansion) . 2012_07_06 by SCB   
      include 'trancntl.inc' ! STUCK ROD
      include 'thexpan.inc'
! added end      
!
      nzp1=nz+1
      if(nx.eq.1 .and. hxa(1).eq.0) hxa(1)=1
      if(ny.eq.1 .and. hya(1).eq.0) hya(1)=1
      if(nz.eq.1 .and. hza(1).eq.0) hza(1)=1
      
! radial node index
      la=0
      nxy=0
      do ja=1,nya
         nrowxa(ja)=nxea(ja)-nxsa(ja)+1
         do ia=nxsa(ja),nxea(ja)
            la=la+1
            nodela(ia,ja)=la
            latoia(la)=ia
            latoja(la)=ja
            nmeshxy=nmeshx(ia)*nmeshy(ja)
            latol(la)%nfm=nmeshxy
! 2014_12_18 . scb            
            !call dmalloc(latol(la)%fm,nmeshxy)
            !call dmalloc(latol(la)%ij,nmeshx(ia),nmeshy(ja))
            allocate(latol(la)%fm(nmeshxy))
            allocate(latol(la)%ij(nmeshx(ia),nmeshy(ja)))
            latol(la)%fm=0
            latol(la)%ij=0
! added end            
            nxy=nxy+nmeshxy
         enddo
      enddo
!
      call dmalloc0(ltoi,0,nxy)
      call dmalloc0(ltoj,0,nxy)
      call dmalloc(ltola,nxy)
!
      do la=1,nxya
        latol(la)%nfm=0
      enddo
!
      l=0
      j=0
      do ja=1,nya
        nxst=0
        do ia=1,nxsa(ja)-1
          nxst=nxst+nmeshx(ia)
        enddo
        nxst=nxst+1
        nxet=nxst
        do ia=nxsa(ja),nxea(ja)
          nxet=nxet+nmeshx(ia)
        enddo
        nxet=nxet-1
        nrowxt=nxet-nxst+1
        do ji=1,nmeshy(ja) 
          j=j+1
          nxs(j)=nxst
          nxe(j)=nxet
          nrowx(j)=nrowxt
	    i=0
!	    if(nxsa(ja) .ne. 1) then
            do ia=1,nxsa(ja)-1
	        i=i+nmeshx(ia)
	      enddo
!	    endif
          do ia=nxsa(ja),nxea(ja)
            la=nodela(ia,ja)
            do ii=1,nmeshx(ia)
              i=i+1
              l=l+1
              nodel(i,j)=l
              ltoi(l)=i
              ltoj(l)=j
              latol(la)%nfm=latol(la)%nfm+1
              latol(la)%fm(latol(la)%nfm)=l
              ltola(l)=la
            enddo
          enddo
        enddo
      enddo

!
! determine coordinates of node in a assembly
      do la=1,nxya
        ia=latoia(la)
        ja=latoja(la)
        inum=1
        do imeshy=1,nmeshy(ja)
          do imeshx=1,nmeshx(ia)
            latol(la)%ij(imeshx,imeshy)=latol(la)%fm(inum)
            inum=inum+1
          enddo
        enddo            
      enddo 

!
! determine surface numbers 
      do ia=1,nxa
        nysa(ia)=1
        nyea(ia)=nya
        lan=0
        do ja=1,nya
          la=nodela(ia,ja)
          if(la.ne.0 .and. lan.eq.0) exit
          lan=la
        enddo
        nysa(ia)=ja
        las=0
        do ja=nya,1,-1
          la=nodela(ia,ja)
          if(la.ne.0 .and. las.eq.0) exit
          las=la
        enddo
        nyea(ia)=ja
        nrowya(ia)=nyea(ia)-nysa(ia)+1
      enddo
!
      do i=1,nx
        nys(i)=1
        nye(i)=ny
        ln=0
        do j=1,ny
          l=nodel(i,j)
          if(l.ne.0 .and. ln.eq.0) exit
          ln=l
        enddo
        nys(i)=j
        ls=0
        do j=ny,1,-1
          l=nodel(i,j)
          if(l.ne.0 .and. ls.eq.0) exit
          ls=l
        enddo
        nye(i)=j
        nrowy(i)=nye(i)-nys(i)+1
      enddo
!
      nsurfa=0
      do ja=1,nya
        nsurfa=nsurfa+nrowxa(ja)+1
      enddo
      nsurfxa=nsurfa
      do ia=1,nxa
        nsurfa=nsurfa+nrowya(ia)+1
      enddo
!
      nsurf=0
      do j=1,ny
        nsurf=nsurf+nrowx(j)+1
      enddo
      nsurfx=nsurf
      do i=1,nx
        nsurf=nsurf+nrowy(i)+1
      enddo      
!
      call dmalloc(neibra,nrdir2,nxya)
      call dmalloc(lsfca,nrdir2,nxya)
      call dmalloc(nodema,nsurfa)
      call dmalloc(nodepa,nsurfa)
      call dmalloc(jsfcda,nsurfa)
!
      call dmalloc(neibr,nrdir2,nxy)
      call dmalloc(lsfc,nrdir2,nxy)
!
      do l=1,nxy
        i=ltoi(l)
        j=ltoj(l)
        neibr(1,l)=nodel(i-1,j)
        neibr(2,l)=nodel(i+1,j)
        neibr(3,l)=nodel(i,j-1)
        neibr(4,l)=nodel(i,j+1)
      enddo

!
      call dmalloc(hmesh,ndirmax,nxy,nz)
	call dmalloc(rhmesh,ndirmax,nxy,nz)

! axial neighbor node numbering
      k=1
      neibz(1,k)=0
      kb=k
      do k=2,nz
         neibz(2,kb)=k
         neibz(1,k)=kb
         kb=k
      enddo 
      k=nz
      neibz(2,k)=0
!
      kfineb=1
      ktoka(0)=1
      ktoka(nz+1)=nz
      do ka=1,nza
         fnmesh=nmeshz(ka)
         kfinee=kfineb+nmeshz(ka)-1
         kfmb(ka)=kfineb
         kfme(ka)=kfinee
         ktoka(kfineb:kfinee)=ka
         kfineb=kfinee+1
      enddo
!
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
!
      elev=0
      do ka=1,nza
         elev=elev+hza(ka)*0.5
         heleva(ka)=elev
         elev=elev+hza(ka)*0.5
      enddo
!
      elev=0
      l=1
      do k=1,nz
         elev=elev+hmesh(ZDIR,l,k)*0.5
         helev(k)=elev
         elev=elev+hmesh(ZDIR,l,k)*0.5
      enddo
      
! added in ARTOS ver. 0.2 (Stuck rod) . 2012_07_06 by SCB   
! find stuck rod
      lst=0
      irod=0
      do la=1,nxya
        irodtyp1=irodtyp(la)
        if(irodtyp1.lt.0) then
          lst=lst+1
          lstroda(lst)=la
          lstrodb(lst)=abs(irodtyp1)
        endif
        if(irodtyp(la).ne.0) then
          lcrbptr(abs(irodtyp1))=lcrbptr(abs(irodtyp1))+1
          irod=irod+1
        endif
      enddo
      nstrod=lst
      nrodpos=irod

      row=0
      do irod=1,nrodtyp
        lcrbptr(irod)=sum(lcrbptr(irod-1:irod))
        do la=1,nxya
        irodtyp1=irodtyp(la)
          if(irod.eq.abs(irodtyp1)) then
            row=row+1
            rodtola(row)=la
          endif
        enddo
      enddo
! added end      
	
      return
!-------------------------------------
!------------ entry ------------------
!-------------------------------------
      entry initgeom1

        do iat=1,nassytyp
           if(iffuela(iat)) exit
        enddo
      
        do ka=1,nza
          if(iffuelc(icompnumz(ka,iat))) exit
        enddo
        kfbega=min(ka,nza)
        kfbeg=kfmb(kfbega)
        
        do ka=nza,1,-1 
           if(iffuelc(icompnumz(ka,iat))) exit
        enddo
        kfenda=min(ka,nza)
        kfend=kfme(kfenda)

        coreheight=0
        do ka=kfbega,kfenda
           coreheight=coreheight+hza(ka)
        enddo
!starting and ending radial coordinates
! added in ARTOS ver. 0.2 . 2012_07_06 by SCB   
        jfbeg=ny
        jfend=1
! added end
        nfuela=0
        nxsf=1   ! 2015_03_09 . scb
        nxsfa=1   ! 2015_03_09 . scb
        do j=1,ny
            do i=nxs(j),nxe(j)
                l=nodel(i,j)
                la=ltola(l)
                ja=latoja(la)
                ia=latoia(la)
                iassytyp1=iassytyp(la)
                if(iffuela(iassytyp1)) then
                    jfbeg=min(jfbeg,j)  ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB   
                    nxsf(j)=i
                    nxsfa(ja)=ia
                    exit
                endif
            enddo

            if(nxsf(j).le.0) cycle
            
            !do i=nxe(j),nxs(ja),-1
            do i=nxe(j),nxs(j),-1   ! 2015_06_23 . scb fixed bug
                l=nodel(i,j)
                la=ltola(l)
                ja=latoja(la)
                ia=latoia(la)
                iassytyp1=iassytyp(la)
                if(iffuela(iassytyp1)) then
                    jfend=max(jfend,j)  ! added in ARTOS ver. 0.2 . 2012_07_06 by SCB   
                    nxef(j)=i
                    nxefa(ja)=ia
                    exit
                endif
            enddo
        enddo
! 2012_08_22 . scb
        do la=1,nxya
            iassytyp1=iassytyp(la)
            if(iffuela(iassytyp1)) nfuela=nfuela+1
        enddo
! added end
      
!volume
        call dmalloc(volnode,nxy,nz)
        call dmalloc(volnodea,nxya,nza)
!
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

        do k=1,nz
          do l=1,nxy
            volnode(l,k)=hmesh(XDIR,l,k)*hmesh(YDIR,l,k)*hmesh(ZDIR,l,k)
          enddo
        enddo
!
! determine ndir
        if(ny.eq.1 .and. nz.eq.1) then
          ndir=1
        elseif(ny.gt.1 .and. nz.eq.1) then
          ndir=2
        else
          ndir=3
        endif
!
! r/b map - maprb
!
!    - map for r/b sweep
        call dmalloc(maprb,nxy,2)
        irbls=1
        irble=2
        irbld=1      
        do irbk=1,2
          inow=1  
          do irbl=irbls,irble,irbld 
            irbx=irbl       
            do iy=1,ny      
              do ix=irbx,nxe(iy),2
                if(ix.ge.nxs(iy).and.ix.le.nxe(iy)) then
                  nnode=nodel(ix,iy) 
                  maprb(inow,irbk)=nnode
                  inow=inow+1
                endif
              enddo
              if(irbx.eq.1) then
                irbx=2
              elseif(irbx.eq.2) then
                irbx=1
              endif
            enddo    
          enddo
          irbls=2
          irble=1
          irbld=-1
        enddo

    ! count corners
        nassy=0
        ncorn=0
        do ja=1,nya
          nassyrow=0
          do ia=1,nxa
            if(nodela(ia,ja).ne.0) then
              nassy=nassy+1
              nassyrow=nassyrow+nmeshx(ia) 
            endif 
          enddo
          nassyrow=nassyrow+1
          if(ja.eq.1) then
              ncorn=nassyrow*(nmeshy(ja)+1)
          else
            if(nassyrowd.gt.nassyrow) nassyrowd=nassyrow
            ncorn=ncorn-nassyrowd+nassyrow*(nmeshy(ja)+1)
          endif
          nassyrowd=nassyrow
        enddo
         
      return
        
    end subroutine
