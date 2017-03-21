    subroutine multiplyff
      
      use allocs
      use param
      
      include 'global.h'
      include 'files.h'
      include 'geom.h'
      include 'xsec.h'
      include 'ffdm.h'
      include 'nodal.h'
      include 'pinpr.h'
      include 'ff.h'
      include 'pow.h'
      
      character*80 form
      real :: avgpinp(ng),powsum(ng),fnormg(ng),avgasyp(ng)

      fuelfac=1 !264.0d0/289.0d0
      
      do la=1,nxya
        powvalr(:,:,la) = 0
      enddo
      
      do la=1,nxya
        ia=latoia(la)
        ja=latoja(la)
        fnodes=latol(la)%nfm
        iat=iassytyp(la)
        if(.not.iffuela(iat)) cycle
        
        npinpermeshx = npin/nmeshx(ja)
        npinpermeshy = npin/nmeshy(ia)

        do k=kfbeg, kfend
          ka=ktoka(k)
          hzk=hmesh(ZDIR,la,ka)/coreheight
          ic=icompnumz(ka,iat)
          iffs=iffset(ic)
          powtot=0

          irodtyp1=irodtyp(la)
          rf = 0
          if(irodtyp1.ne.0) rf=rodfrac(k,irodtyp1)
          rf1=1-rf

! 2015_07_31 . scb recover MOXBENCH          
!          do jp=1,npin
!            li=1
!            if(jp.gt.npinpermeshy) li=3
!            do ip=1,npin
!              if(mod(li,2) .ne. 0 .and. ip.gt.npinpermeshx) li=li+1
!            
!              do m=1,ng
!! control rod
!! powtot : assembly total power
!!#ifdef MOXBENCH
!!                iffrod=iffset(ic+20) ! sucks code.
!!               
!!                ffrod =rf1*ff(ip,jp,m,iffs)+rf*ff(ip,jp,m,iffrod)
!!                powval(ip,jp,m)=phih(ip,jp,la,k,m)*xskpf(m,latol(la)%fm(li),k)*ffrod
!!#else              
!!FIXME don't use xskpf(m,latol(la)%fm(1),k)               
!                powval(ip,jp,m)=phih(ip,jp,la,k,m)*xskpf(m,latol(la)%fm(li),k)*ff(ip,jp,m,iffs)
!!#endif                
!                powtot(m)=powtot(m)+powval(ip,jp,m) 
!              enddo ! of group                               
!            enddo ! of ixp
!          enddo ! of iyp
          

          do jp=1,npin
            li=1
            if(jp.gt.npinpermeshy) li=3
            do ip=1,npin
              if(mod(li,2) .ne. 0 .and. ip.gt.npinpermeshx) li=li+1
            
              do m=1,ng
                if(ixsecver.eq.2) then
                  iffrod=iffset(ic+20) ! sucks code.               
                  ffrod =rf1*ff(ip,jp,m,iffs)+rf*ff(ip,jp,m,iffrod)
                else
                  ffrod=ff(ip,jp,m,iffs)
                endif
                
                powval(ip,jp,m)=phih(ip,jp,la,k,m)*xskpf(m,latol(la)%fm(li),k)*ffrod
        
                powtot(m)=powtot(m)+powval(ip,jp,m) 
              enddo ! of group                               
            enddo ! of ixp
          enddo ! of iyp
! added end          

!         powsum : group-wise total assembly power
          powsum=0
          do li=1,latol(la)%nfm
            l=latol(la)%fm(li)
            do m=1,ng
              powsum(m)=powsum(m)+phif(m,l,k)*xskpf(m,l,k)
            enddo
          enddo
          
!         avgpinp : pin-averaged power
!         avgasyp : absolute node-averaged power
          fnorm=0
          do m=1,ng
            avgpinp(m)=powtot(m)/fntpin
            avgasyp(m)=plevel0*powsum(m)/fnodes
            fnorm=fnorm+avgasyp(m)
            fnormg(m)=avgasyp(m)/avgpinp(m)
          enddo
          rnorm=1/fnorm

!         powval1 : absolute pin-power in a assembly
!         powvalr : absolute radial pin-power
          pmax=0
          do jp=1,npin
            do ip=1,npin
              powval1(ip,jp)=0
              do m=1,ng
                powval1(ip,jp)=powval1(ip,jp)+powval(ip,jp,m)*fnormg(m)
              enddo
              powvalr(ip,jp,la)=powvalr(ip,jp,la)+powval1(ip,jp)*hzk
              if(powval1(ip,jp).gt.pmax) then
                pmax=powval1(ip,jp)
                ipmax=ip
                jpmax=jp
              endif
            enddo
          enddo
          
! store peak pin power
          pppeak(la,k)=pmax
! write pin power distribution for the node          
!!define WRITE_PPR
!#ifdef WRITE_PPR
!          if(ia .gt. nxa/2 .and. ja.gt.nya/2) then
!            write(io15,'("Assembly Coordinate (i,j): ",2i4
!     +,                " , Plane Index & Height: ",i4,f9.3
!     +,                " , Normalization Factor: ",1p,e12.5)')
!     +                 ia,ja,k,hmesh(ZDIR,la,ka),fnorm
!            is=1;js=1;
!            if(mod(nxa,2).ne.0 .and. ia.eq.(nxa/2+1)) is=npin/2+1
!            if(mod(nya,2).ne.0 .and. ja.eq.(nya/2+1)) js=npin/2+1
!            do jp=1,js-1
!              write(io15,'(f7.4,$)') (0.0,ip=1,npin)
!              write(io15,*)            
!            enddo
!            do jp=js,npin
!              if(is.ne.1) write(io15,'(f7.4,$)') (0.0,ip=1,is-1)
!              write(io15,'(f7.4,$)') (powval1(ip,jp)*fuelfac,ip=is,npin)
!              write(io15,*)
!            enddo
!          endif
!#endif          
        enddo !k

        pmax=0
        sum1=0
        do jp=1,npin
           do ip=1,npin
              sum1=sum1+powvalr(ip,jp,la)
              if(powvalr(ip,jp,la).gt.pmax) then
                 pmax=powvalr(ip,jp,la)
                 ipmax=ip
                 jpmax=jp
              endif
           enddo
        enddo
        fnorm=sum1/fntpin
        rnorm=1/fnorm
        pppeak(la,0)=pmax
!FIXME uncomment
        write(io15,'("Assembly Coordinate (i,j): ",2i4,", Plane Index & Height: ",i4,f9.3, " , Normalization Factor: ",1p,e12.5)') &
              ia,ja,0,coreheight,fnorm
        write(form,'("(7x,",i5,"i7)")') npin
        write(io15,form) (ip,ip=1,npin)
        write(form,'("(i7,",i5,"f7.4)")') npin
        do jp=1,npin
           write(io15,form) jp,(powvalr(ip,jp,la),ip=1,npin)
        enddo
        write(io15,'("  Peak at (i,j): ",2i4," , Value: ",f10.4)') ipmax,jpmax,pmax*rnorm
  !
        write(io16,*)'Assembly Coordinate ',ia,ja 
        write(io16,*)
        do jp=1,npin
          do ip=1,npin
            write(io16,'(f7.4,$)') powvalr(ip,jp,la)*rnorm*fuelfac
          enddo
          write(io16,*)
        enddo
        write(io16,*)
      enddo !la
!
      if(nz .ne. 1) then
        do la=1,nxya
          ia=latoia(la)
          ja=latoja(la) 
          iat=iassytyp(la)
          if(.not.iffuela(iat)) cycle
          if(ia .gt. nxa/2 .and. ja.gt.nya/2) then
            write(io15,*)'Horizontally Averaged Power, Assembly Coordinate ',ia,ja 
            is=1;js=1;
            if(mod(nxa,2).ne.0 .and. ia.eq.(nxa/2+1)) is=npin/2+1
            if(mod(nya,2).ne.0 .and. ja.eq.(nya/2+1)) js=npin/2+1
            do jp=1,js-1
              write(io15,'(f7.4,$)') (0.0,ip=1,npin)
              write(io15,*)            
            enddo
            do jp=js,npin
              if(is.ne.1) write(io15,'(f7.4,$)') (0.0,ip=1,is-1)
              write(io15,'(f7.4,$)') (powvalr(ip,jp,la)*fuelfac,ip=is,npin)
              write(io15,*)
            enddo
          endif
        enddo   
      endif
      close(io15)
!!FIXME undefine MATLAB_DEBUG       
!!define MATLAB_DEBUG     
!#ifdef MATLAB_DEBUG
!        iodbg=2001
!        write(filename(1), '("debug/",a,"_POWER.out")') trim(caseid)
!        open(iodbg, file=trim(filename(1)),status='unknown')
!      
!        do ja=nya/2+1,nya
!          do iy=1,npin
!            if(nxsfa(ja) .gt. nxefa(ja)) cycle
!            do ia=(nxefa(ja)+nxsfa(ja))/2+1,nxefa(ja)
!              la=nodela(ia,ja)
!              do ix=1,npin
!                write(iodbg,'(f7.4,$)') powvalr(ix,iy,la)*fuelfac
!              enddo
!            enddo !ia
!            write(iodbg,*)
!          enddo !iy
!        enddo !ja
!        close(iodbg)
!#endif

      return
      
    end subroutine
