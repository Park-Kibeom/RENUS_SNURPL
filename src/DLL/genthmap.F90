! added in ARTOS ver. 0.2 ( for System-Neutronics Coupling ). 2012_07_03 by SCB
    subroutine genthmap(nt,mit,vmgf,iok)
!
      use param
      use chanmap
! 
      include 'global.h'
      include 'geom.h' 
      include 'thlink.inc'
      character*128 oneline
      character*4 astr
      data eps/0.0001/
!
      real :: vmgf(nt)
      integer :: mit(nxy,nz)
!
      allocate(mapnr(nxy),mapnz(nz),mapn(nt))
      allocate(mapthr(nchanmars),mapthz(nlevel),mapth(nthn))
      allocate(rvolthn(nthn),qrel(nthn))
      allocate(dmbar(nlevel),tmbar(nlevel),deltm(nlevel))
!
      open(157,file=filemap,err=1001,status='old')
            
1123  continue    ! 2014_04_09 . scb   
!
! read surface fuel temp weighting factor
      read(157,*)
      read(157,*) tfws
!
! read radial mapping data
      read(157,*)
      mapthr(:)%n=0
      mxchan=0
!
      do l=1,nxy
         read(157,'(a128)') oneline
         nf=nfields2(oneline)           
         if(mod(nf,2)/=1) then
            write(astr,*) l
            oneline='Data Pairs Not Formed for Node'//astr
            write(906,*) oneline
            iok=0
            return
         else
            nf=(nf-1)/2
            mapnr(l)%n=nf
            allocate(mapnr(l)%id(nf))
            allocate(mapnr(l)%frac(nf))
            read(oneline,*) idum,(mapnr(l)%id(i),mapnr(l)%frac(i),i=1,nf)
            fracsum=0
            do i=1,nf
               ichan=mapnr(l)%id(i)
               fracsum=fracsum+mapnr(l)%frac(i)
	           mxchan=max(mxchan,ichan)
               mapthr(ichan)%n=mapthr(ichan)%n+1
            enddo
            if(fracsum.lt.1-eps .or. fracsum.gt.1+eps) then 
               write(astr,'(i4)') l
               oneline='Data Pairs Not Formed for Radial Node'//astr
               write(906,*) oneline
               iok=0
               return
            endif                 
         endif   
      enddo
!       
      if(mxchan/=nchanmars) then
	     iok=0
         oneline="Channel Number Mismatch/Check MAS_MAP"
         write(906,*) oneline
         iok=0
         return
      endif
!
      do ichan=1,nchanmars
  	     allocate(mapthr(ichan)%id(mapthr(ichan)%n))
  	     allocate(mapthr(ichan)%frac(mapthr(ichan)%n))
      enddo
!
      mapthr(:)%n=0
      do l=1,nxy
         do i=1,mapnr(l)%n
            ichan=mapnr(l)%id(i)
            mapthr(ichan)%n=mapthr(ichan)%n+1	      
            mapthr(ichan)%id(mapthr(ichan)%n)=l
            mapthr(ichan)%frac(mapthr(ichan)%n)=mapnr(l)%frac(i)            
         enddo
      enddo
!
! read axial mapping data
      read(157,*)
      mapthz(:)%n=0
      !mlevel=0
      mxlevel=0   ! 2013_09_17 . scb
      do k=1,nz
         read(157,'(a128)') oneline
         nf=nfields2(oneline)           
         if(mod(nf,2)/=1) then
            write(astr,*) l
            oneline='Fractions not Added to Unity for Plane'//astr
            write(906,*) oneline
            iok=0
            return
         else
            nf=(nf-1)/2
            mapnz(k)%n=nf
            allocate(mapnz(k)%id(nf))
            allocate(mapnz(k)%frac(nf))
            read(oneline,*) idum,(mapnz(k)%id(i),mapnz(k)%frac(i),i=1,nf)
            fracsum=0
            do i=1,nf
               ilevel=mapnz(k)%id(i)
               fracsum=fracsum+mapnz(k)%frac(i)
	           mxlevel=max(mxlevel,ilevel)
               mapthz(ilevel)%n=mapthz(ilevel)%n+1
            enddo
            if(fracsum.lt.1-eps .or. fracsum.gt.1+eps) then 
	           write(astr,*) k
               oneline='Fractions not Added to Unity for Plane'//astr
               write(906,*) oneline
               iok=0
               return
            endif                 
         endif   
      enddo
!       
      if(mxlevel/=nlevel) then
         iok=0
         oneline="Plane Number Mismatch/Check Mapping File"
         write(906,*) oneline
         iok=0
         return
      endif
!
      do ilevel=1,nlevel
  	     allocate(mapthz(ilevel)%id(mapthz(ilevel)%n))
  	     allocate(mapthz(ilevel)%frac(mapthz(ilevel)%n))
      enddo
!
      mapthz(:)%n=0
      do k=1,nz
         do i=1,mapnz(k)%n
            ilevel=mapnz(k)%id(i)
            mapthz(ilevel)%n=mapthz(ilevel)%n+1	      
            mapthz(ilevel)%id(mapthz(ilevel)%n)=k
            mapthz(ilevel)%frac(mapthz(ilevel)%n)=mapnz(k)%frac(i)            
	     enddo
      enddo
!
      do l=1,nxy
         do k=1,nz
            m=mit(l,k)
            nthm=mapnr(l)%n*mapnz(k)%n
            allocate(mapn(m)%id(nthm))
            allocate(mapn(m)%frac(nthm))
         enddo
      enddo
!
! assign 3D mapping info
      lth=0
      mapn(:)%n=0
      do ichan=1,nchanmars
         do ilevel=1,nlevel
            lth=lth+1
            mapth(lth)%n=mapthr(ichan)%n*mapthz(ilevel)%n
            allocate(mapth(lth)%id(mapth(lth)%n))
            allocate(mapth(lth)%frac(mapth(lth)%n))
            lthn=0
            vol=0
            lm=0
            do kn=1,mapthz(ilevel)%n
               k=mapthz(ilevel)%id(kn)
               fracz=mapthz(ilevel)%frac(kn)
               do ln=1,mapthr(ichan)%n
                  l=mapthr(ichan)%id(ln)
                  m=mit(l,k)
                  lthn=lthn+1
                  mapth(lth)%id(lthn)=m
                  frac=mapthr(ichan)%frac(ln)*fracz
                  mapth(lth)%frac(lthn)=frac                                 
                  mapn(m)%n=mapn(m)%n+1
                  in=mapn(m)%n
                  mapn(m)%id(in)=lth
                  mapn(m)%frac(in)=frac
                  vol=vol+vmgf(m)*frac
               enddo
            enddo
            rvolthn(lth)=1/vol
         enddo
      enddo						
      iok=1

      if(icouple.eq.2) then
         ifcobra=.TRUE.
      else
         ifcobra=.FALSE.
      endif
!
      close(157)
      return
!
!      do k=1,nz
!         do l=1,nxy
!            m=mit(l,k)
!            write(906,901) k,l,m,mapn(m)%n,(mapn(m)%id(i),mapn(m)%frac(i),i=1,mapn(m)%n)
!         enddo
!      enddo
!      lth=1
!      do ichan=1,nchanmars
!         do ilevel=1,nlevel
!            write(906,902) ichan,ilevel,lth,mapth(lth)%n,1/rvolthn(lth),(mapth(lth)%id(i),mapth(lth)%frac(i),i=1,mapth(lth)%n)
!            lth=lth+1
!         enddo
!      enddo
!  901 format(4i5,200(i5,f5.2),/,20X,200(i5,f5.2))
!  902 format(4i5,f15.3,200(i5,f5.2),/,35X,200(i5,f5.2))
!

! 2014_04_09 . scb
1001 continue
     
     filemap='MAS_MAP'
     open(157,file=filemap,err=1002,status='old')
     
     goto 1123
     
1002 continue
     
     stop 'there is no mapping file named "artos.map" or "MAS_MAP"!!'
! added end     
     
    end subroutine
!=======================================================================================!
    function nfields2(aline)
      character*128  aline
      logical nonblankd,nonblank
!
      nonblankd=.false.
      n=0
      do i=128,1,-1
         if(aline(i:i).eq.' ') then
            nonblank=.false.
         else
            nonblank=.true.
            if(aline(i:i).eq.'!') then
               n=0
               nonblankd=.true.
            endif
         endif
         if(nonblankd.xor.nonblank) n=n+1
         nonblankd=nonblank
      enddo
      nfields2=n/2
      return
    end function