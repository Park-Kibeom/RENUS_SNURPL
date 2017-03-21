    subroutine writeff
! generate assembly-wise group constants and discontinuity factors
! this assumes the geometry represents only an assembly.
! if the geometry is a quarter core, the position must be let in 1st quarter.
      use param
      use allocs
      
      include 'global.h'
      include 'files.h'
      include 'geom.h'
      include 'xsec.h'
      include 'ffdm.h'
!      
      real,pointer,dimension(:,:), save :: ff, powvola
      logical, save :: first=TRUE

      real ::  xstra(ng), xsaa(ng), xsnfa(ng), xskpa(ng), &
               xssma(ng,ng), adf(4,ng), rtotphifvol(ng)
      
      if(first) then
        call dmalloc(ff,ng,nxya)
        call dmalloc(powvola,ng,nxya)
        first=FALSE
      endif
      
! generate assembly-wise xsec        
      k=1;ka=ktoka(k);
      xssma=0
      do m=1,ng
        totphifvol=0 !assembly-wise phif*vol
        totvol=0 !assembly-wise vol
        do la=1,nxya
          vola=0
          powvola(m,la)=0
          do li=1,latol(la)%nfm
            l=latol(la)%fm(li)
            vola=vola+volnode(l,k)
            phifvol=volnode(l,k)*phif(m,l,k)            
            totphifvol=totphifvol+phifvol
            
            xstrf1=1/(3*xsdf(m,l,k))
            xstra(m)=xstra(m)+xstrf1*phifvol
            xsaa(m)=xsaa(m)+xsaf(m,l,k)*phifvol
            xsnfa(m)=xsnfa(m)+xsnff(m,l,k)*phifvol
            powvol=xskpf(m,l,k)*phifvol
            xskpa(m)=xskpa(m)+powvol
            powvola(m,la)=powvola(m,la)+powvol
            do ms=xssfs(m,l,k), xssfe(m,l,k) 
              xssma(ms,m)=xssma(ms,m)+xssf(ms,m,l,k)*phif(ms,l,k)*vol
            enddo
          enddo !li
          powvola(m,la)=powvola(m,la)/vola
          totvol=totvol+vola
        enddo !la
        
        rtotphifvol(m)=1/totphifvol
        
        xstra(m)=xstra(m)*rtotphifvol(m)
        xsaa(m)=xsaa(m)*rtotphifvol(m)
        xsnfa(m)=xsnfa(m)*rtotphifvol(m)
        xskpa(m)=xskpa(m)*rtotphifvol(m)

! generate form function
        do la=1,nxya
          ff(m,la)=powvola(m,la)*totvol/(xskpa(m)*totphifvol)
        enddo

!generate ADF
! north
        ravgphif=rtotphifvol(m)*totvol

        irdir=3
        vols=0
        adf(irdir,m)=0
        j=1
        do i=nxs(j), nxe(j)
          l=nodel(i,j)
          la=ltola(l)
          ln=neibr(4,l)
          vol=volnode(l,k)
          vols=vols+vol
          adf(irdir,m)=adf(irdir,m)+(7*phif(m,l,k)-phif(m,ln,k))*RSIX*vol
        enddo  
        adf(irdir,m)=adf(irdir,m)/vols*ravgphif

! west
        irdir=1
        vols=0
        adf(irdir,m)=0
        i=1
        do j=nys(i), nye(i)
          l=nodel(i,j)
          la=ltola(l)
          ln=neibr(2,l)
          vol=volnode(l,k)
          vols=vols+vol
          adf(irdir,m)=adf(irdir,m)+(7*phif(m,l,k)-phif(m,ln,k))*RSIX*vol
        enddo  
        adf(irdir,m)=adf(irdir,m)/vols*ravgphif

        if(isymang.eq.90) then
          adf(2,m)=adf(1,m)
          adf(4,m)=adf(3,m)
        else
! south
          irdir=4
          vols=0
          adf(irdir,m)=0
          j=ny
          do i=nxs(j), nxe(j)
            l=nodel(i,j)
            la=ltola(l)
            ln=neibr(3,l)
            vol=volnode(l,k)
            vols=vols+vol
            adf(irdir,m)=adf(irdir,m)+(7*phif(m,l,k)-phif(m,ln,k))*RSIX*vol
          enddo  
          adf(irdir,m)=adf(irdir,m)/vols*ravgphif
! east
          irdir=2
          vols=0
          adf(irdir,m)=0
          i=nx
          do j=nys(i), nye(i)
            l=nodel(i,j)
            la=ltola(l)
            ln=neibr(1,l)
            vol=volnode(l,k)
            vols=vols+vol
            adf(irdir,m)=adf(irdir,m)+(7*phif(m,l,k)-phif(m,ln,k))*RSIX*vol
          enddo  
          adf(irdir,m)=adf(irdir,m)/vols*ravgphif
        endif ! isymang != 90
      enddo !m
      
      do md=1,ng
        do ms=xssfs(md,l,k), xssfe(md,l,k) 
          xssma(ms,md)=xssma(ms,md)*rtotphifvol(ms)
        enddo       
      enddo
      
! write xsec to file
      ifile=2000
      open(ifile, file=trim(caseid)//'.xsc', status='unknown')

      write(ifile, '(" base")')
      do m=1,ng
        write(ifile, '(4x,i2,1x,1p,5e13.6)') m,xstra(m),xsaa(m),xsnfa(m),xskpa(m),xschif(m,1,1)
      enddo
      do ms=1,ng
        write(ifile,'(4x,i2)') ms
        do md=1,ng
          write(ifile, 600) xssma(ms,md)
        enddo
        write(ifile, *)
      enddo
      write(ifile, *)
      
      write(ifile, '(" adf")')
      do m=1,ng
        write(ifile, '(4x,i2,1x,1p,4e13.6)') m,adf(1,m),adf(2,m),adf(3,m),adf(4,m)
      enddo
      close(ifile)
      
      open(ifile, file=trim(caseid)//'.ff', status='unknown')
      do m=1,ng 
        write(ifile, '("GROUP ",i2)') m
        do ja=1,nya
          do ia=nxsa(ja),nxea(ja)
            la=nodela(ia,ja)
            write(ifile, 600) ff(m,la)
          enddo
          write(ifile, *)
        enddo
      enddo
      close(ifile)

      return
            
600   format(1x,1p,e13.6,$)       

	end subroutine
	