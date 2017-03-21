! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
    subroutine readgeomhex
!
      USE MASTERXSL
      use param
      use allocs
!
      include 'global.h'      
      include 'cards.h'
      include 'files.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'xsec.h'
      include 'ffdm.h'
      include 'thexpan.inc'
      logical ifnumeric
!
      INDEV2 = 900
      indev=io5
      iffile=FALSE
      
      if(probe.eq.DOT)  probe=''   ! 2014_12_17 . scb      
!
  100 continue
      do while (probe.ne.DOT)
         read(indev,'(a512)',end=1000) oneline
         write(io8,'(a)') trim(oneline)
         if(probe.eq.BANG .or. oneline.eq.BLANK .or. ifnumeric(oneline)) cycle
         if(probe.eq.DOT .or. probe.eq.SLASH) exit
         if(probe.ne.BLANK) then
            backspace(indev)
            backspace(io8)
            go to 2000
         endif
         read(oneline,*) cardname
         call toupper(cardname)
         if(cardname.eq.'FILE') then
            indev=io5+100
            call openlf(indev,oneline)
            iffile=TRUE
            go to 100
         endif
         ndataf=nfields(oneline)-1
         select case(cardname) 
            case('GEO_DIM')
            case('FDM')
               ifhexfdm=TRUE
            case('SP3')
               read(oneline,*) cardname,ifhexsp3
            case('RAD_CONF')
               labeg=1
               read(oneline,*) cardname,isymang
               if(ndataf.ge.2) read(oneline,*) cardname,isymang,symopt
               call toupper(symopt)
               symopt=trim(symopt)
               if (isymang.eq.30) then
                  isymtype=2
                  isymmetry=12
               elseif (isymang.eq.60) then
                  isymmetry=6
               elseif (isymang.eq.120) then
                  isymtype=1
                  isymmetry=3
               endif
               call readlay(indev,isymang)
            case('ALBEDO_R')
               read(oneline,*) cardname,temp
               albr(:)=temp        
            case('ALBEDO_ZB')
               read(oneline,*) cardname,temp 
               albzb(:)=temp        
            case('ALBEDO_ZT')
               read(oneline,*) cardname,temp  
               albzt(:)=temp                            
            case('BC_RZ')
              ! 2015_08_05 . scb changed BC index
               read(oneline,*) cardname,ibcr,ibczl,ibczr
               select case(ibcr) 
                  case(0) 
                     albr=0
                  case(1)
                     !albr=big
                     albr=half
                  case(2)
                     !albr=half
                     albr=big
               end select
               select case(ibczl) 
                  case(0) 
                     albzb=0
                  case(1)
                     !albzb=big
                     albzb=half
                  case(2)
                     !albzb=half
                     albzb=big
               end select
               select case(ibczr) 
                  case(0) 
                     albzt=0
                  case(1)
                     !albzt=big
                     albzt=half
                  case(2)
                     !albzt=half
                     albzt=big
               end select 				          
            case('GRID_HEX')
               read(oneline,*) cardname,hf2f
               n3ang=6*4**(ndivhs-1)
               sqrt3=sqrt(3.)
               rt3=sqrt(3.)
               rsqrt3=1/sqrt3
               hside=hf2f*rsqrt3
            case('GRID_Z')
               read(oneline,*) cardname,(hza(k),k=1,nza)
               hza0(:)=hza(:)
            case('NEUTMESH_Z')
               read(oneline,*) cardname,nmeshz(1:nza)
            case('ASSY_TYPE')
               read(oneline,*) cardname,iat
               read(oneline,*) cardname,iatt,(icompnumz(k,iat),k=1,nza)
               do k=1,nza
                  if(iffuelc(icompnumz(k,iat))) iffuela(iat)=TRUE
               enddo
            case('ROD_CONF')              
! 2014_08_20 . scb              
               ndataf=nfields(oneline)-1  
               if(ndataf.eq.1) then       
                 read(oneline,*) cardname,CRTYPE     ! 2014_08_12 . PKB
                 CALL TOUPPER(CRTYPE)
                 IF(CRTYPE.EQ.'CORNER') THEN
                  FLAGCCR = .TRUE.
                  CALL READCCR(INDEV)     
                 ELSE
                   stop
                 END IF
               else                 
                  call readbankh(1,indev,isymang)
               endif
! added end
            case('ROD_TYPE')
               read(oneline,*) cardname,irodtyp1,fullpos,stepsize,step
               rodfullpos(irodtyp1)=fullpos
               rodstep(irodtyp1)=step
               rodstep0(irodtyp1)=step   ! 2014_12_22 . scb
               rodstepsize(irodtyp1)=stepsize
               rodstepsize0(irodtyp1)=stepsize
! 2014_09_01 . SCB
            case('CRADF') 
               read(oneline,*) cardname, flagcradf
               if(flagcradf) then
                 if(.not.flagccr .or. .not.flagmasxsl) stop 'CRADF is valid only for corner control rod model'
                 ndataf=nfields(oneline)-1  
                 if(ndataf.gt.1) read(oneline,*) cardname, flagcradf, fncradf

                 ICRADF=0
                 OPEN(INDEV2,FILE=FNCRADF,STATUS='OLD',IOSTAT=ICRADF)
                 IF(ICRADF.NE.0) STOP 'Input for CRADF does not exist !'
                 
                 CALL READCRADF(INDEV2,ICRADF)      
                 
                 close(INDEV2)   ! 2014_12_16 . scb
               endif
! added end
            case default
               call terminate(trim(cardname)//' Card Not Allowed')
         end select
      enddo
 1000 continue

! reset the input device after done with local file
      if(iffile) then
         close(indev)
         indev=io5
         iffile=FALSE
         ! return to the next card in the input file
         go to 100
      endif

2000  continue

      ndir=3

! initialze assembly weight array
      val=1
      wtass=val

! generate node and other ordering and neighbor info.
      call orderhex
      if(ifhexfdm) call orderhexfdm

      nzp1=nz+1
      if(nz.eq.1 .and. hza(1).eq.0) hza(1)=1

      call dmalloc(ltola,nxy)
      lfa=0
      l=0
      wtasssum=0      
      do la=1,nassy
        l=l+1
        ltola(l)=la
        if(iffuela(iassytyp(la))) then
          lfa=lfa+1
          lfaptr(lfa)=lfa
          lfatol(lfa)=l
          wtasssum=wtasssum+wtass(la) 
        endif
      enddo
      ! end of change
      nfuelfa=lfa
      nfuel=nfuelfa
      nchan=nfuelfa

#ifdef EXTTH
      if( extth ) nchan=nxy
#endif
!
#ifdef FLOT
      do j=1,ny
        do i=nxs(j),nxe(j),nxskip
          lfa=ltolfa(nodel(i,j))

          if((ltolfa(nodel(i-nxskip,j)).eq.0) .and. lfa.ne.0) nxfas(j)=i
          if((ltolfa(nodel(i+nxskip,j)).eq.0) .and. lfa.ne.0) nxfae(j)=i
          nxas(j)=nxs(j)
          nxae(j)=nxe(j)
        enddo
      enddo   
#endif        

! establish bank number to assembly number correspondence
      if(nrodtyp.ne.0) call readbankh(2,indev,isymang)
!
! initial axial flux shape
!
! -- flat flux initial guess using 4:1 fast to thermal flux ratio
!
! for axial plots
#ifdef FLOT
      if(flxlevel.eq.0.) flxlevel=1
      do k=1,nz
        axf(1,k)=0.8*flxlevel
        axf(2,k)=0.2*flxlevel
      enddo
#endif
!
      call dmalloc(hmesh,ndirmax,nxy,nz)
      call dmalloc(rhmesh,ndirmax,nxy,nz)
      call dmalloc0(hz,0,nz+1)
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

      hz=1.
      do k=1,nz
        ka=ktoka(k)
        la=ltola(l)
        hz(k)=hza(ka)/nmeshz(ka)
      enddo


      do iat=1,nassytyp
        if(iffuela(iat)) exit
      enddo

      if(iat.gt.nassytyp) iat=nassytyp  ! 2013_10_08 . scb
      
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

      ! mesh sizes
      do k=1,nz
        ka=ktoka(k)
        do l=1,nxy
          la=ltola(l)
          hmesh(ZDIR,l,k)=hza(ka)/nmeshz(ka)
          rhmesh(ZDIR,l,k)=1/hmesh(ZDIR,l,k)
        enddo
      enddo
! calculate node volumes and axial coordinate of top of each plane
      call dmalloc(volnode,nxy,nz)
      call dmalloc(volnodea,nxya,nza)
      volnodea=0
      hexarea=1.5d0*hf2f*hside
      do k=1,nz
        ka=ktoka(k)
        vol=hexarea*hmesh(ZDIR,1,k)
        do l=1,nxy
          volnode(l,k)=vol*wtass(l)
          volnodea(l,ka)=volnodea(l,ka)+volnode(l,k)
        enddo
      enddo

      hactive=0
      volfuel=0
      do k=kfbeg,kfend
        hactive=hactive+hz(k)
        do la=1,nassy
          if(iffuela(iassytyp(la)))  volfuel=volfuel+volnode(la,k)
        enddo
      enddo
!
      wtasssum=wtasssum*hactive !temp
!

#ifdef HEX
! fill a nonzero for the meshes outside the problem geometry
      hx(0)=one
      hx(nxp1)=one
      hy(0)=one
      hy(nyp1)=one
      hz(0)=one
      hz(nzp1)=one
!
! for axial plots
      zb(0)=0
      do k=1,nz
        zb(k)=znode(k)
        zcent(k)=0.5*(zb(k-1)+zb(k))
      enddo
!
! determine domain structure, dummies for serial runs
      ndomx=1
      ndomy=1
      ndomz=1
      ndomxy=ndomx*ndomy
      ndom=ndomz*ndomxy
      if(ndom.ne.ncpus) stop ' Domain Configuration Inconsistent'
      idom=iam+1
      idomxy=mod(iam,ndomxy)
      idomx=mod(idomxy,ndomx)+1
      idomy=idomxy/ndomx+1
      idomz=iam/ndomxy+1
      idomz=1 !temp for DEC alpha
      idomxy=idomxy+1
#endif
      return
    end subroutine