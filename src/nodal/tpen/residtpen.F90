! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    function residtpen()

	    use const
	    use geomhex
	    use tpen
	    use xsec
	    use sfam, only : reigv
	    implicit none

      real,allocatable,save,dimension(:) :: cntz0i,cntz1i,hflxfn,cnti(:,:)
      real,pointer,save,dimension(:,:,:) :: hflxtmp
	    logical,save :: first=TRUE

	    real :: residtpen, totresid, totfsrc, arear, adfl, &
	            reflrhom, vol, sumlkr, sumlkz, sumss,      &
	            sumfs, sumes, sumtxs, sumlk, resid
	    integer :: iz, ih, it, nn, ig, ig2, neigdn, neigup


      if(first) then
        allocate(cnti(ng,ntph),cntz0i(ng),cntz1i(ng),hflxfn(ng))
        cnti=0
        cntz0i=0
        cntz1i=0
        hflxfn=0
        if(ng.ne.2) then
           hflxtmp=>hflxf
        else
           hflxtmp=>hflx
        endif
        first=.FALSE.
      endif

      totresid=0
      totfsrc=0

      do iz=1,nz
        arear=hside*hz(iz)
        do ih=1,nassy
          !if(nassy.gt.1) then   ! 2012_12_17 . scb
            do it=1,6
              nn=neignd(it,ih)
              if(nn.eq.0) then
                do ig=1,ng
!                  adfl=xsadf(ig,ih,iz)
                  adfl=1.
                  reflrhom=(adfl-2*alxr)/(adfl+2*alxr)
                  cnti(ig,it)=cnto(ig,it,ih,iz)*reflrhom
                enddo
              else
                do ig=1,ng
                  cnti(ig,it)=cnto(ig,neigjin(it,ih),nn,iz) 
                enddo
              endif
            enddo
          !endif   ! 2012_12_17 . scb
          !if(nz.gt.1) then   ! 2012_12_17 . scb
            neigdn=neigz(1,iz) 
            neigup=neigz(2,iz)
            if(neigdn.ne.0) then 
              do ig=1,ng
                cntz0i(ig)=cntzo(ig,2,ih,neigdn)
              enddo
            else
              do ig=1,ng
                cntz0i(ig)=reflratzb*cntzo(ig,1,ih,iz)
              enddo
            endif
            if(neigup.ne.0) then 
              do ig=1,ng
                cntz1i(ig)=cntzo(ig,1,ih,neigup)
              enddo
            else
              do ig=1,ng
                cntz1i(ig)=reflratzt*cntzo(ig,2,ih,iz)
              enddo
            endif
          !endif   ! 2012_12_17 . scb
          vol=volnode(ih,iz)
          sumlkr=0
          sumlkz=0
          sumss=0
          sumfs=0
          sumes=0
          sumtxs=0
          do ig=1,ng
            do it=1,6
              sumlkr=sumlkr+cnto(ig,it,ih,iz)-cnti(ig,it)
            enddo
            sumlkz=sumlkz+cntzo(ig,1,ih,iz)-cntz0i(ig) &
                         +cntzo(ig,2,ih,iz)-cntz1i(ig)
            do ig2=xssfs(ig,ih,iz),ig-1
              sumss=sumss+xssf(ig2,ig,ih,iz)*hflxtmp(ig2,ih,iz)
            enddo
            do ig2=ig+1,xssfe(ig,ih,iz)
              sumss=sumss+xssf(ig2,ig,ih,iz)*hflxtmp(ig2,ih,iz)
            enddo
            sumfs=sumfs+xsnff(ig,ih,iz)*hflxtmp(ig,ih,iz)
!            sumes=sumes+sefve(ig,ih,iz)
            sumtxs=sumtxs+xstf(ig,ih,iz)*hflxtmp(ig,ih,iz)
          enddo
          sumlk=(sumlkr*arear+sumlkz*hexarea)*wtass(ih)
          resid=(sumfs*reigv+sumss+sumes-sumtxs)*vol-sumlk
          totresid=totresid+resid*resid
          sumfs=sumfs*vol
          totfsrc=totfsrc+sumfs*sumfs
        enddo
      enddo
!      reigvn=(-resid+sumfs*reigv)/sumfs
      residtpen=sqrt(totresid/totfsrc)

      return
    end
