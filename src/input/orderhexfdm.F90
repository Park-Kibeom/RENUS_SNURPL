! added in ARTOS ver. 0.2 ( for Hexagonal ). 2012_07_03 by SCB
    subroutine orderhexfdm

      use param
      use allocs

      include 'global.h'
      include 'itrcntl.h'
      include 'files.h'
      include 'geom.h'
      include 'geomh.h'
      include 'geomhfc.h'
      include 'xsec.h'
      include 'defhex.h'
      include 'defsfc.h'
      include 'defpnt.h'
      include 'hexfdm.h'

      integer :: nend(ndiv), nbeg(ndiv), ntri(ndiv), ndiff(ndiv)
      integer :: iston(ndiv,6)

      call dmalloc(neigtri, 3, nsub, nxy)
      call dmalloc(neigtria,3, nsub, nxy)
!
      do i=1,ndiv
         nbeg(i)=1
         ntri(i)=6*(ndiv-i+1)
         nend(i)=6*(ndiv-i+1)
         ndiff(i)=2*(2*ndiv-i)+1
      enddo
      do i=2,ndiv
         nbeg(i)=nbeg(i)+nend(i-1)+nend(i)
         nend(i)=nend(i-1)+2*nend(i)
      enddo

      m=0
      do idir=1,6
         do i=1,ndiv
            m=m+1
            iston(i,idir)=m
         enddo
      enddo

! numbering
      do l=1,nxy
         m=0
         do iring=1,ndiv
            ibeg=nbeg(iring)
            iend=nend(iring)
            ! outer trianlges
            do idir=1,6
               m2=ntri(iring)-idir
               m3=ntri(iring)-idir+1
               do idiv=1,ndiv-iring+1
                  m=m+1               
                  i1=m-ntri(iring)
                  ia1=l
                  if(iring.eq.1) then
                     if(neignd(idir,l).eq.0) then         
                        i1=0
                        ia1=0
                     else
                        if(isymtype.eq.2) then
                           isfc=neigjin(idir,l)
                           if((isfc.eq.6 .and. idir.eq.6) .or. &
                              (isfc.eq.1 .and. idir.eq.5) .or. &
                              (isfc.eq.5 .and. idir.eq.1) .or. &
                              (isfc.eq.2 .and. idir.eq.3) .or. &
                              (isfc.eq.3 .and. idir.eq.2) ) then                        
                              i1=iston(idiv,isfc)
                              ia1=neignd(idir,l)                        
                           else
                              i1=iston(ndiv-idiv+1,isfc)
                              ia1=neignd(idir,l)
                           endif
                        else
                           isfc=neigjin(idir,l)
                           i1=iston(ndiv-idiv+1,isfc)
                           ia1=neignd(idir,l) 
                        endif                              
                     endif
                  endif
                  i2=m+m2
                  i3=m+m3
                  if(iring.lt.ndiv) then
                     if(idiv.eq.1) then
                        i2=m-1                
                     elseif(idiv.eq.ndiv-iring+1) then
                        i3=m+1  
                     endif
                  else
                     i2=m-1  
                     i3=m+1  
                  endif

                  if(m.eq.ibeg) i2=iend
                  if(m.eq.iend) i3=ibeg
                  neigtri(1,m,l)=i1
                  neigtri(2,m,l)=i2
                  neigtri(3,m,l)=i3

                  neigtria(1,m,l)=ia1
                  neigtria(2,m,l)=l
                  neigtria(3,m,l)=l
               enddo ! idiv
            enddo ! idir

            ! inner triangles
            do idir=1,6
               m2=ntri(iring)-idir
               m3=ntri(iring)-idir+1
               do idiv=1,ndiv-iring
                  m=m+1
                  i1=m+ntri(iring+1)
                  i2=m-m2
                  i3=m-m3

                  neigtri(1,m,l)=i1
                  neigtri(2,m,l)=i2
                  neigtri(3,m,l)=i3

                  neigtria(1,m,l)=l
                  neigtria(2,m,l)=l
                  neigtria(3,m,l)=l
               enddo ! idiv
            enddo ! idir
         enddo ! iring
      enddo ! l

      return
    end subroutine
