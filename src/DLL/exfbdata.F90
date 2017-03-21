! added in ARTOS ver. 0.2 ( for System-Neutronics Coupling ). 2012_07_03 by SCB
! exfbdata.f - exchange T/H feedback variables between Neutronics and T/H code
!
    subroutine exfbdata(ifin,nthdim,fbdata,bankpos,iok)
!
!     ifin = true for incoming data processing
!          = false for outgoing data processing
      use param
      use linkdt
      use chanmap
      use wordsize    ! 2014_08_24 . scb
!
      include 'global.h'
      include 'geom.h'  
      include 'xsec.h' 
      include 'pow.h'  
      include 'ffdm.h' 
      include 'thgeom.inc'  
      include 'thfuel.inc'
      include 'thop.inc'  
      include 'thfdbk.inc'
      include 'thlink.inc' 
      include 'thexpan.inc' 
          
      logical ifin
      type(THDIM) nthdim    
      type(FBVAR) ::  fbdata(0:6000)
      !real bankpos(*)
      real(NBF) bankpos(*)    ! 2014_08_24 . scb
!
      logical,save :: first=TRUE, first2=TRUE
      integer :: nt

      real,pointer,dimension(:,:) ::  xsmack, phim
      real,pointer,dimension(:) :: vmgf, tfnew, tmnew, dmnew, rfnod, ppmnod
      integer,pointer ::  mit(:,:), mapm2d(:)
      
      common / exfb_global1 / xsmack, phim
      common / exfb_global2 / vmgf, tfnew, tmnew, dmnew, rfnod, ppmnod
      common / exfb_global3 / mit, mapm2d
      
! 2014_04_08 . scb
      integer,pointer :: ltoidll(:), ltojdll(:), nodeldll(:,:)
      common / geom_hex_dll / ltoidll, ltojdll, nodeldll
      integer,pointer :: nchanscb(:,:)
      common / geom_hex_dll2 / nchanscb
! added end      

! 2014_12_16 . scb
      type(MMIDATA)  art2mmi2      ! mmi display data
      common / artosmmi / art2mmi2
! added end

!
      nt=nxy*nz ! mh  
      if(first) then
        allocate(xsmack(nt,ng), phim(nt,ng), vmgf(nt), mit(nxy,nz), mapm2d(0:nxy), tfnew(nt), tmnew(nt), dmnew(nt), ppmnod(nt), rfnod(nt))
        allocate(nchanscb(0:nzth,0:nchan+1)) ! 2014_04_09 . scb
      endif
!
      if(first) then
! 2014_04_08 . scb
         m=0
         if(hex) then
           allocate(ltoidll(nxy),ltojdll(nxy),nodeldll(-3:nx+4,-1:ny+2))
           ltoidll=0
           ltojdll=0
           nodeldll=0
           
           do iy=-1,ny+2
             do ix=-3,nx+4
               if(nodel(ix,iy).ne.0) then
                 m=m+1
                 ltoidll(m)=ix
                 ltojdll(m)=iy
                 nodeldll(ix,iy)=m
               endif
             enddo
           enddo
         endif
         
         if(m.ne.nxy) stop 'mistmatch in the # of assy! error in exfbdata subroutine'
! added end         
                
         m=0
         do k=1,nz
            do l=1,nxy
                m=m+1
                vmgf(m)=volnode(l,k)
                mit(l,k)=m
            enddo
         enddo
!
! 2014_04_09 . scb commented         
         !do l=0,nxy
         !   mapm2d(l)=ltochan(l)
         !enddo
!
         nchanmars=nthdim%n_chan
         nlevel=nthdim%n_level
         nbank=nthdim%n_bank
! 2014_04_23 . scb     
         nthn=nchanmars*nlevel
         power=fbdata(0)%qrel	 ! total core power in watts
         qtripavg=power/volfuel ! power density in watts/cc
         call genthmap(nt,mit,vmgf,iok)        
!
         tfwc=1-tfws ! centerline fuel temperature weighting factor
      
         first=FALSE
      endif
!
! mapping variables
      if(rect) then
        do ig=1,ng
          m=0
          do k=1,nz
              do l=1,nxy
                  m=m+1
                  xsmack(m,ig)=xskpf(ig,l,k)
                  phim(m,ig)=phif(ig,l,k)
              enddo
          enddo
        enddo   
! 2014_04_08 . scb        
      else
        do ig=1,ng
          m=0
          do k=1,nz
              do l=1,nxy
                  !m=m+1
                  m=nodeldll(ltoi(l),ltoj(l))
                  m=m + nxy*(k-1)
                  
                  xsmack(m,ig)=xskpf(ig,l,k)
                  phim(m,ig)=phif(ig,l,k)
              enddo
          enddo
        enddo   
      endif
! added end      
!
      call exfbdata1(ifin,fbdata,bankpos,xsmack,phim,rfnod,vmgf,mapm2d,mit,  &
                      tfnew,tmnew,dmnew,ppmnod,nt)
!
! mapping variables
      if(ifin) then
! 2014_04_23 . scb
! 2014_08_08 . scb
        if(nbank .eq. nrodtyp) then
          do irod=1,nrodtyp
            rodstep(irod)=bankpos(irod)
          enddo
          if(first2) then
            first2=.false.
            !write(21,'(2a10)')  'Time(sec)', 'Bankpos(1:12)'
          endif
          
          !write(21,'(13f10.3)') trtime, bankpos(1:12)   ! 2014_08_23 . scb
        else
          print *, '# of control rod bank defined in system and neutronics codes is different ! '  ! 2014_08_08 . scb
        endif
! added end
! added end

! 2014_04_09 . scb
        tdopl=0.d0
        tcool=0.d0
        dcool=0.d0
        nchanscb=0
! added end        
        if(rect) then
          m=0
          do k=1,nz
	          do l=1,nxy
			        m=m+1
			        lth=ltochan(l)
			        kth=ktokth(k) 
! 2014_04_09 . scb              
			        !if(lth.ne.nchan+2) then
		         !   tdopl(kth,lth)=sqrt(tfnew(m))
		         !   tcool(kth,lth)=tmnew(m)
		         !   dcool(kth,lth)=dmnew(m)*1000. 
			        !endif	
			        if(lth.le.nchan+1) then
                nchanscb(kth,lth)=nchanscb(kth,lth)+1
		            tdopl(kth,lth)=tdopl(kth,lth)+sqrt(tfnew(m))
		            tcool(kth,lth)=tcool(kth,lth)+tmnew(m)
		            dcool(kth,lth)=dcool(kth,lth)+dmnew(m)*1000. 
                ! non-fuel 영역에서의 T/H값 저장
              elseif(lth.ne.nchan+1) then
                stop 'channel index has problem...'
              endif
! added end              
	          enddo
          enddo
! 2014_04_08 . scb          
        else
          do k=1,nz
	          do l=1,nxy
			        m=nodeldll(ltoi(l),ltoj(l))
              m=m + nxy*(k-1)
			        lth=ltochan(l)
			        kth=ktokth(k)           
              
			        if(lth.le.nchan+1) then
                nchanscb(kth,lth)=nchanscb(kth,lth)+1
		            tdopl(kth,lth)=tdopl(kth,lth)+sqrt(tfnew(m))
		            tcool(kth,lth)=tcool(kth,lth)+tmnew(m)
		            dcool(kth,lth)=dcool(kth,lth)+dmnew(m)*1000. 
              elseif(lth.ne.nchan+1) then
                stop 'channel index has problem...'
              endif	            
	          enddo
          enddo
        endif
        
        do kth=1,nzth
          do lth=1,nchan+1
            tdopl(kth,lth)=tdopl(kth,lth)/nchanscb(kth,lth)
            tcool(kth,lth)=tcool(kth,lth)/nchanscb(kth,lth)
            dcool(kth,lth)=dcool(kth,lth)/nchanscb(kth,lth)
          enddo
        enddo
        
        do kth=1,nzth
          do lth=1,nchan
            if(nchanscb(kth,lth) .gt. 1) stop 'there is some problem'
          enddo
        enddo
! added end

! average temperature
        tfuelavg=0.d0
        tcoolavg=0.d0
        dcoolavg=0.d0   ! 2014_09_05 . scb
        tinavg=0.d0
        toutavg = 0.d0   ! 2014_12_17 . scb
        do l=1,nchan
          do k=kfsth,kfeth
	          temp=tdopl(k,l)*tdopl(k,l)
            tfuelavg=tfuelavg+temp*hzth(k)
	          tcoolavg=tcoolavg+tcool(k,l)*hzth(k)
	          dcoolavg=dcoolavg+dcool(k,l)*hzth(k)   ! 2014_09_05 . scb
          enddo
	        tinavg=tinavg+tcool(1,l)
	        toutavg=toutavg+tcool(nzth,l)
	      enddo
        tfuelavg=tfuelavg/(nchan*0.01*coreheight)
        tcoolavg=tcoolavg/(nchan*0.01*coreheight)
        dcoolavg=dcoolavg/(nchan*0.01*coreheight)   ! 2014_09_05 . scb
        tinavg=tinavg/nchan
        toutavg=toutavg/nchan
            
! 2014_12_16 . scb         
#ifdef DLL_TASS
        art2mmi2%tincore = int(tinavg)
        art2mmi2%toutcore = int(toutavg)
        
        art2mmi2%crbpos(1)=int(bankpos(1))
        art2mmi2%crbpos(2)=int(bankpos(4))
        
        !write(1223, *) tinavg, toutavg, art2mmi2%tincore, art2mmi2%toutcore
#endif         
! added end

        !tcoolavg=tcoolavg+273.15    ! 2014_10_06 . scb
        tinavg=tinavg+273.15
        toutavg=toutavg+273.15
	
        tinavg=0.5*(tinavg+toutavg)
      endif

      return	  
    end subroutine
!=======================================================================================!
!=======================================================================================!
!=======================================================================================!
    subroutine exfbdata1(ifin,fbdata,bankpos,xsmack,phi,rfnod,vmgf,mapm2d,mit,  &
                           tfnew,tmnew,dmnew,ppmnod,nt)
!
      use param
      use linkdt
      use chanmap
      use wordsize    ! 2014_08_24 . scb
!
      include 'global.h'
      include 'xsec.h'
      include 'geom.h'
      include 'pow.h'
      include 'thgeom.inc'
      include 'thfuel.inc'
      include 'thfdbk.inc'
      include 'thlink.inc'
      include 'thop.inc'
      include 'cntl.h'  ! 2014_04_10 . scb
!
      common /decheat/alphabar !hj 24mar99
      logical ifin
      type(FBVAR) fbdata(0:6000)
      !real bankpos(*)
      real(NBF) bankpos(*)    ! 2014_08_24 . scb      
      save twritten
      data twritten/-1.0/
!
      integer nt
      integer nscb ! 2013_09_16 . scb
      real xsmack(nt,ng),phi(nt,ng),rfnod(nt),vmgf(nt), &
           tfnew(nt),tmnew(nt),dmnew(nt),ppmnod(nt), &
           fhc(1000,2), fgc(1000,2), fp(2) 
      integer :: mit(nxy,nz)   
      
! 2014_04_08 . scb
      integer,pointer :: ltoidll(:), ltojdll(:), nodeldll(:,:)
      common / geom_hex_dll / ltoidll, ltojdll, nodeldll
! added end      

      integer,save :: ith2n=0 , in2th=0
!
!***************************************************************
      if(ifin) then
!***************************************************************
         if(ilink.le.1) tpr=sstime
         if(ilink.ge.2) tpr=trtime
         if(.not.(icouple.eq.2 .and. nssstep.gt.0)) then !for mars th
            dmnew=0.
            tmnew=0.
            tfnew=0.
            ppmnod=0.
            lth=1
            !sumdm=0
            !sumtm=0
            !sumtf=0
      			tfmax=0.
            tmmax=0.  ! 2014_10_06 . scb
            do ichan=1,nchanmars
!
! obtain node average values from the exit values            
               if (thcode(1:4) .eq. 'MARS') then
                  lthd=lth
                  dmbar(1)=fbdata(lthd)%dm
                  tmbar(1)=fbdata(lthd)%tm
                  deltm(1)=0 
                  do ilevel=2,nlevel
                     lthn=lthd+1
                     dmbar(ilevel)=0.5*(fbdata(lthn)%dm+fbdata(lthd)%dm) 
                     tmbar(ilevel)=0.5*(fbdata(lthn)%tm+fbdata(lthd)%tm)
                     deltm(ilevel)=fbdata(lthn)%tm-tmbar(ilevel)
                     lthd=lthn
                  enddo   
               elseif (thcode(1:6) .eq. 'RETRAN') then
                  lthn=lth
                  do ilevel=1,nlevel
                     dmbar(ilevel)=fbdata(lthn)%dm 
                     tmbar(ilevel)=fbdata(lthn)%tm 
                     deltm(ilevel)=0
                     lthn=lthn+1
                  enddo   
               endif
!
! assign t/h values            
               do ilevel=1,nlevel
 	                tfeff=tfwc*fbdata(lth)%tfc+tfws*fbdata(lth)%tfs
 	                tfeff=tfeff-deltm(ilevel)
                  !sumdm=sumdm+dmbar(ilevel)
                  !sumtm=sumtm+tmbar(ilevel)
                  !sumtf=sumtf+tfeff
				          tfmax=max(tfmax,fbdata(lth)%tfc)
                  tmmax=max(tmmax,fbdata(lth)%tm)    ! 2014_10_06 . scb
                  
! 2014_04_15 . scb commented            
                  !if(ilevel.eq.1 .or. ilevel.eq.nlevel .or. ichan.eq.nchanmars) then
                  !  tfeff=tmbar(1)+273.15   !in the reflector region                                      
                  !endif
                  
                  nscb = mapth(lth)%n
                  do l=1,nscb
                    m=mapth(lth)%id(l)
                    frac=mapth(lth)%frac(l)
 	                  dmnew(m)=dmnew(m)+dmbar(ilevel)*frac
 	                  tmnew(m)=tmnew(m)+tmbar(ilevel)*frac
 	                  tfnew(m)=tfnew(m)+tfeff*frac
 	                  ppmnod(m)=ppmnod(m)+fbdata(lth)%ppm*frac  !temp
	                enddo
                  lth=lth+1
               enddo
	          enddo
			      tfmax=tfmax-273.15
            !sumdm=sumdm/nthn
            !sumtf=sumtf/nthn
            !if(nssstep.eq.0) then
            if(nssstep.eq.0 .and. .not.rstrt) then  ! 2014_04_10 . scb
               ppmnod=ppm
               go to 210
            endif
         endif
         
! 2014_04_09 . scb commented          
! inlet T/H boundary condition
!         lth=nthn
!         do ichcob=1,nchancob
!            fhc(ichcob,1)=0
!            fgc(ichcob,1)=0
!         enddo
!         sumpres=0
!         sumfrac=0
!         do ichan=1,nchanmars-1
!            lth=lth+1
!            do l=1,mapthr(ichan)%n
!!               ichcob=mapm2d(mapthr(ichan)%id(l))
!               frac=mapthr(ichan)%frac(l)
!               fhc(ichcob,1)=fhc(ichcob,1)+fbdata(lth)%qrel*frac
!               fgc(ichcob,1)=fgc(ichcob,1)+fbdata(lth)%dm*frac
!               sumpres=sumpres+fbdata(lth)%ppm*frac  !pressure averagin
!               sumfrac=sumfrac+frac  !temp
!            enddo
!         enddo
!         fp(1)=sumpres/sumfrac      !average pressure
!         fp(1)=fp(1)*fnormpexit     !normalized pressure
!         fp(2)=fp(1)
!         do ichcob=1,nchancob
!            fhc(ichcob,1)=fhc(ichcob,1)*fnormhin
!            fgc(ichcob,1)=fgc(ichcob,1)*fnormgin
!            fhc(ichcob,2)=fhc(ichcob,1)
!            fgc(ichcob,2)=fgc(ichcob,1)
!         enddo
!
! assign ppm of t/h to Neutronics code
         ppmnod=0
         lth=1
         do ichan=1,nchanmars
            do ilevel=1,nlevel
              ln=mapth(lth)%n
               do l=1,ln
                  m=mapth(lth)%id(l)
                  frac=mapth(lth)%frac(l)
                  ppmnod(m)=ppmnod(m)+fbdata(lth)%ppm*frac  !temp
               enddo
               lth=lth+1
            enddo
         enddo
!            
         ppmmars=0
         volppm=0
         do m=1,nt
            if(ppmnod(m).ne.0) then
               volppm=volppm+vmgf(m)
               ppmmars=ppmmars+vmgf(m)*ppmnod(m)
            endif
         enddo
         if(volppm.eq.0) volppm=1
         ppmmars=ppmmars/volppm
         if(ilink.ge.1) ppm=ppmmars  !temp
         
210      continue
 
! 2014_04_09 . scb commented         
!#define WTRTHBC  
!#ifdef WTRTHBC  
!         if(ilink.ge.2 .and. t.ne.twritten) then                     
!            write(1277,'(f10.3,f15.6)') t,fp(1)
!            write(1278,'(f10.3)') t
!            write(1278,'(5f15.6)') (fhc(ichcob,1),ichcob=1,nchancob)
!            write(1279,'(f10.3)') t
!            write(1279,'(5f15.6)') (fgc(ichcob,1),ichcob=1,nchancob)
!            twritten=t
!         endif
!#endif
         !m=mapth(nlevel/2)%id(1)
         !if (thcode(1:4) .eq. 'MARS') then
         !   write(907,'(" MARS T-H: ",i5,2f10.4,f10.2,5x,0p6f10.3)')  &
         !   nmars,tpr,sumdm,sumtf,fp(2),fhc(1,2),fgc(1,2),fbdata(nlevel/2)%ppm,ppmnod(m),ppm
         !elseif (thcode(1:6) .eq. 'RETRAN') then
         !   write(907,'(" RETRAN T-H: ",i5,2f10.4,f10.2,5x,0p6f10.3)')  &
         !   nmars,tpr,sumdm,sumtf,fp(2),fhc(1,2),fgc(1,2),fbdata(nlevel/2)%ppm,ppmnod(m),ppm
         !endif
         
        ith2n=ith2n+1
        !write(25,*) 'ith2n : ',ith2n
        !do lth=0,nthn         
        !  WRITE(25,'(i5,5f13.5)') lth,fbdata(lth)%dm,fbdata(lth)%tm,fbdata(lth)%tfc,fbdata(lth)%tfs,fbdata(lth)%ppm
        !enddo         
!***************************************************************
      else
!***************************************************************
! normalize flux in eigenvalue calculation
         if(ilink.eq.0) fbdata(0)%ppm=ppm
         fbdata(0)%qrel=power*plevel  ! returns total power in watts

         qtripavg=power/volfuel ! power density in watts/cc
         qrel=0.d0
         idch = 0  ! 2014_04_08 . scb no decay heat 
         !if(idch.eq.0) then ! no decay heat
#define RELP
#ifdef RELP
            do m=1,nt
	            nr=mod(m,nxy)
	            if(nr.eq.0) nr=nxy
	            lchan=ltochan(nr)
              if(hex) lchan=ltochan(nodel(ltoidll(nr),ltojdll(nr)))   ! 2014_04_08 . scb
	            k=(m-nr)/nxy+1
	            kchan=ktokth(k)
	            if(lchan.le.nchan) then
                ln=mapn(m)%n
				        do l=1,ln
					        lth=mapn(m)%id(l)
					        frac=mapn(m)%frac(l)
! 2014_04_15 . scb                  
					        !qrel(lth)=qrel(lth)+relp(kchan,lchan)*frac
					        qrel(lth)=qrel(lth)+relp(kchan,lchan)*frac*vmgf(m)
! added end                  
				        enddo
	            endif
            enddo
            pnorm=1.d0/qtripavg
#else
            do m=1,nt
              ln=mapn(m)%n
               do l=1,ln
                  lth=mapn(m)%id(l)
                  frac=mapn(m)%frac(l)
				          do im=1,ng
                    qrel(lth)=qrel(lth)+(xsmack(m,im)*phi(m,im))*vmgf(m)*frac
				          enddo
               enddo
            enddo
            pnorm=ipart*fnorm/qtripavg
#endif
         !else
! 2014_04_09 . scb commented       
            !stop 'decay heat is not implemented'
!            fnorma=fnorm*alphabar !hj 24mar99
!#ifdef MSLBEDT
!            sumfp=0
!            sumdp=0
!#endif
!            do m=1,nt
!              ln=mapn(m)%n
!              do l=1,ln
!                  lth=mapn(m)%id(l)
!                  frac=mapn(m)%frac(l)
!    	            qrel(lth)=qrel(lth)+((XSMACK(M,1)*PHI(M,1)+XSMACK(M,2)*PHI(M,2))*fnorma+DCHK*RFNOD(M))*VMGF(M)*frac
!#ifdef MSLBEDT
!                  sumfp=sumfp+(XSMACK(M,1)*PHI(M,1)+XSMACK(M,2)*PHI(M,2))*FNORMA*VMGF(M)*FRAC
!                  sumdp=sumdp+DCHK * RFNOD(M) * VMGF(M)* FRAC
!#endif
!	            enddo
!            enddo
!	          pnorm=ipart/qtripavg
	      !endif
!#ifdef MSLBEDT
!        write(711,*) t,sumfp,sumdp,(sumfp+sumdp)/2772.e4
!#endif
!! calculate hot pin power
!        do ichan=1,nchanmars-1
!	        volratio=mapthr(ichan)%frac(1)+mapthr(ichan+1)%frac(1)
!	        if((dabs(volratio-1).lt.0.0001).and.(mapthr(ichan)%frac(1)/mapthr(ichan+1)%frac(1).gt.25)) then
!!#define PIN	           	        
!#ifdef PIN
!          if(pintobox(mchan).ne.1.0) then
!            lth=(ichan-1)*nlevel
!!                 write(906,*) "PTB",nmars,mchan,pintobox(mchan)
!            do ilevel=1,nlevel
!              lth=lth+1
!              lthpin=lth+nlevel
!              pownode=qrel(lth)+qrel(lthpin)
!              qrel(lthpin)=qrel(lthpin)*pintobox(mchan)
!              qrel(lth)=pownode-qrel(lthpin)
!            enddo
!          endif
!#endif	      
!	        endif	      
!	      enddo
! obtain relative nodal power	
!#ifdef DEBUG
!        sump=0.0
!        nfa=0
!        do lth=1,nthn
!          sump=sump+qrel(lth)
!          if(qrel(lth).gt.0.0) nfa=nfa+1
!        enddo
!        pnorm=volfuel/sump
!!	       pnorm=nfa/sump
!
!        sumq=0.0
!        do lth=1,nthn
!          fbdata(lth)%qrel=pnorm*qrel(lth)*rvolthn(lth)
!!		       fbdata(lth)%qrel=pnorm*qrel(lth)
!          sumq=sumq+fbdata(lth)%qrel/rvolthn(lth)
!          WRITE(26,*)lth,fbdata(lth)%qrel,qrel(lth)
!        enddo
!        sumq=sumq/volfuel
!#endif
        lth=0
        do ichan=1,nchanmars
          do ilevel=1,nlevel
            lth=lth+1
            lthn=0
            vol=0.d0
            lm=0
            do kn=1,mapthz(ilevel)%n
              k=mapthz(ilevel)%id(kn)
              fracz=mapthz(ilevel)%frac(kn)
              do ln=1,mapthr(ichan)%n
                l=mapthr(ichan)%id(ln)
                m=mit(l,k)
                lthn=lthn+1
                frac=mapthr(ichan)%frac(ln)*fracz
                in=mapn(m)%n
                vol=vol+volnode(l,k)*frac
              enddo
            enddo
            rvolthn(lth)=1.d0/vol
          enddo
        enddo	

        sump=0.d0
        !nfa=0
        do lth=1,nthn
          sump=sump+qrel(lth)
          !if(qrel(lth).gt.0.d0) nfa=nfa+1
        enddo
        pnorm=volfuel/sump

        !sumq=0.d0
        do lth=1,nthn
          fbdata(lth)%qrel=pnorm*qrel(lth)*rvolthn(lth)
          !sumq=sumq+fbdata(lth)%qrel/rvolthn(lth)
        enddo
 
! 2014_04_11 . scb        
#ifdef DLL_TASS
        if(fbdata(0)%dm.lt.0) then
          do lth=1,nthn
            fbdata(lth)%qrel=fbdata(lth)%qrel*plevel
          enddo
        endif
#endif
! added end
!! 2014_04_16 . scb    for dbg
!        do lth=1,nthn         
!          if(fbdata(lth)%qrel.gt.1.e-10) fbdata(lth)%qrel=1.d0
!        enddo
!! added end        
                 
        in2th=in2th+1
        !write(26,*) 'in2th : ',in2th
        !do lth=0,nthn         
        !  WRITE(26,*)lth,fbdata(lth)%qrel,qrel(lth),1.d0/rvolthn(lth)
        !enddo
          

        !sumq=sumq/volfuel

!***************************************************************
      endif
!***************************************************************
      return
    end subroutine
