! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine upddhathexmg

#define p1_nodal_correction   ! 2014_03_12 . scb    

! Update CMFD coefficient(dhat and beta) from Nodal Solution 

	use const
	use geomhex
	use tpen
	use cmfdhex2g
	use xsec
	use sfam, only : phi, fphi

	implicit none

! variables for nem solver in z-direction
! atleak : transverse leakage for axial direction
! atleaku,atleakl : transverse leakage for upper and lower node
! xsdu,xsdl : xsd(diffusion coefficients) for upper and lower node

    integer :: mvz(2)
    data mvz/2,1/

	integer :: iz, ih, ig, ig2, it, isfc, ind, &
	           is1, is2, id1, id2, iz1, iz2
	real :: sumflx, sumcnt, cnto1, cnto2, flxs, cnts, &
	        flx1, flx2, dfdm, w, fdflxs, adfl, albl,  &
	        reflratfl, hzr, xsd0
	real :: reflratz

	integer :: m, l, k, md, mbeg, mend, m2
	real :: xss2t, xssup

	real :: underrelx

#ifdef UNDER_RELAX
	underrelx = 1.0

      dhatd = dhat
      dhatzd= dhatz
      betaphisd = betaphis
	  betaphiszd= betaphisz
#endif

! update D(diffusion coefficient) for CMFD
      do iz=1,nz 
        do ih=1,nassy
          do ig=1,ng2
            sumflx=0
            do ig2=mgb(ig),mge(ig)
              sumflx=sumflx+hflxf(ig2,ih,iz)
            enddo
            hflx(ig,ih,iz)=sumflx
            phi(ig,ih,iz)=hflx(ig,ih,iz)  !temp
            do ig2=mgb(ig),mge(ig)
              fohflx(ig2,ih,iz)=fhflx(ig2,ih,iz)
              fhflx(ig2,ih,iz)=hflxf(ig2,ih,iz)/sumflx
              fphi(ig2,ih,iz)=fhflx(ig2,ih,iz)
            enddo
          enddo
          do it=1,6
            do ig=1,ng2
              sumcnt=0
              do ig2=mgb(ig),mge(ig)
                sumcnt=sumcnt+cnto(ig2,it,ih,iz)
              enddo
              do ig2=mgb(ig),mge(ig)
                focnto(ig2,it,ih,iz)=fcnto(ig2,it,ih,iz)
                fcnto(ig2,it,ih,iz)=cnto(ig2,it,ih,iz)/sumcnt
              enddo
            enddo
          enddo
        enddo
      enddo
      if(nz.gt.1) then
        do iz=1,nz 
          do ih=1,nassy
            do it=1,2
              do ig=1,ng2
                sumcnt=0
                do ig2=mgb(ig),mge(ig)
                  sumcnt=sumcnt+cntzo(ig2,it,ih,iz)
                enddo
                do ig2=mgb(ig),mge(ig)
                  focntzo(ig2,it,ih,iz)=fcntzo(ig2,it,ih,iz)
                  fcntzo(ig2,it,ih,iz)=cntzo(ig2,it,ih,iz)/sumcnt
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
      
#ifdef p1_nodal_correction   ! 2014_03_12 . scb          
      call colxs(fphi)
      call upddtilhex
#endif      

      do isfc=1,nxsfc
        ind=neigsnd(3,isfc)
        if(ind.eq.12) then
          is1=neigsnd(1,isfc)
          is2=neigsnd(2,isfc)
          id1=neigsnd(4,isfc)
          id2=neigsnd(5,isfc)
          do iz=1,nz
            do ig=1,ng2
              cnto1=0
              cnto2=0
              do ig2=mgb(ig),mge(ig)
                cnto1=cnto1+cnto(ig2,id1,is1,iz)
                cnto2=cnto2+cnto(ig2,id2,is2,iz)
              enddo
              flxs=2.*(cnto1+cnto2)
              cnts=cnto1-cnto2
              flx1=hflx(ig,is1,iz)  
              flx2=hflx(ig,is2,iz)  
              dfdm=dfd(ig,isfc,iz)
              xsd0=xsd(ig,is1,iz)
              w=dfdm/xsd0/2.
              fdflxs=(1.-w)*flx1+w*flx2 
              dhat(ig,isfc,iz)= &
                 -(rt3*hside*cnts+dfdm*(flx2-flx1))/(flx2+flx1)
              betaphis(ig,isfc,iz)=2.*(flxs-fdflxs)/(flx1+flx2)
            enddo
          enddo
        else
          is1=neigsnd(ind,isfc)
          id1=neigsnd(ind+3,isfc)
          do iz=1,nz 
            do ig=1,ng2
              cnto1=0
              cnto2=0
              do ig2=mgb(ig),mge(ig)
                cnto1=cnto1+cnto(ig2,id1,is1,iz)
!                adfl=xsadf(ig2,is1,iz)
                adfl=1. 
                reflratfl=(adfl-2*alxr)/(adfl+2*alxr)
                cnto2=cnto2+reflratfl*cnto(ig2,id1,is1,iz)
              enddo
              flxs=2.*(cnto1+cnto2)
              cnts=cnto1-cnto2
              flx1=hflx(ig,is1,iz) 
              dfdm=dfd(ig,isfc,iz)
              dhat(ig,isfc,iz)=(dfdm*flx1-rt3*hside*cnts)/flx1
              if(ind.eq.2) dhat(ig,isfc,iz)=-dhat(ig,isfc,iz)
              betaphis(ig,isfc,iz)=flxs/flx1
            enddo
          enddo
        endif
      enddo
      if(nz.eq.1) goto 999
! update D in z-direction for CMFD
      do iz=1,nz+1 
        ind=neigsndz(3,iz)
        if(ind.eq.12) then
          iz1=neigsndz(1,iz)
          iz2=neigsndz(2,iz)
          id1=neigsndz(4,iz)
          id2=neigsndz(5,iz)
          hzr=hz(iz2)/hz(iz1)
          do ih=1,nassy 
            do ig=1,ng2
              cnto1=0
              cnto2=0
              do ig2=mgb(ig),mge(ig)
                cnto1=cnto1+cntzo(ig2,id1,ih,iz1)
                cnto2=cnto2+cntzo(ig2,id2,ih,iz2)
              enddo
              flxs=2.*(cnto1+cnto2)
              cnts=cnto1-cnto2
              flx1=hflx(ig,ih,iz1)  
              flx2=hflx(ig,ih,iz2) 
              dfdm=dfdz(ig,ih,iz)
              xsd0=xsd(ig,ih,iz1)
              w=dfdm/xsd0/(1.+hzr)
              fdflxs=(1.-w)*flx1+w*flx2 
              dhatz(ig,ih,iz)=-((hz(iz1)+hz(iz2))*cnts*0.5 &
                        +dfdm*(flx2-flx1))/(flx2+flx1)
              betaphisz(ig,ih,iz)=2.*(flxs-fdflxs)/(flx1+flx2)
            enddo
          enddo
        else
          iz1=neigsndz(ind,iz)
          id1=neigsndz(ind+3,iz)
          if(iz.eq.1) then
            do ig=1,ng
              reflratz=reflratzb
            enddo
          else
            do ig=1,ng
              reflratz=reflratzt
            enddo
          endif
          do ih=1,nassy 
            do ig=1,ng2
              cnto1=0
              cnto2=0
              do ig2=mgb(ig),mge(ig)
                cnto1=cnto1+cntzo(ig2,id1,ih,iz1)
                cnto2=cnto2+reflratz*cntzo(ig2,id1,ih,iz1)
              enddo
              flxs=2.*(cnto1+cnto2)
              cnts=cnto1-cnto2
              flx1=hflx(ig,ih,iz1)
              dfdm=dfdz(ig,ih,iz)
              dhatz(ig,ih,iz)=(dfdm*flx1-hz(iz1)*cnts)/flx1
              if(ind.eq.2) dhatz(ig,ih,iz)=-dhatz(ig,ih,iz)
              betaphisz(ig,ih,iz)=flxs/flx1
            enddo
          enddo
        endif
      enddo

  999 continue

#ifdef COND_NODAL
! update old abs xsec
      do k=1,nz
         do l=1,nxy
            xsslocal=0
            do md=mgb(ng),mge(ng)
               do m=mgb(1),mge(1)
                  xsslocal=xsslocal+xssf(m,md,l,k)*fphi(m,l,k)
               enddo
            enddo
            xssup=0
            do md=mgb(1)+1,mge(1)
               do m=mgb(ng),mge(ng)
                  xssup=xssup+xssf(m-1,md,l,k)*fphi(m,l,k)
               enddo
            enddo
            xssup=xssup*phi(2,l,k)/phi(1,l,k)
            do m=1,ng
               xstd(m,l,k)=0
               do mf=mgb(m),mge(m)
                  xstd(m,l,k)=xstd(m,l,k)+xsaf(mf,l,k)*fphi(mf,l,k)
               enddo
            enddo
            xstd(1,l,k)=xstd(1,l,k)+xsslocal-xssup
         enddo
      enddo
#endif

#ifdef DEBUG
      do iz=1,nz 
        do ih=1,nassy
          do ig=1,ng2
		  do ig2=mgb(ig),mge(ig)
            phif(ig2,ih,iz)=hflxf(ig2,ih,iz)  !temp
	        fphi(ig2,ih,iz)=fhflx(ig2,ih,iz)
	      enddo
	      phi(ig,ih,iz)=hflx(ig,ih,iz)
          enddo
        enddo
      enddo
#endif


#ifdef DEBUG
	print *
!	do ih=1,nxy
	print '(a10,10F9.5)', 'Flux' ,hflx(1,1,1:10) 
	print '(a10,10F9.5)', 'Phi'  ,phi(1,1,1:10)
!	print '(a10,10F9.5)', 'Phif' ,phif(1,1,1:10)  
!	enddo
	print '(/a10,10F9.5)', 'Joutl',cntzo(1,1,1,1:10)
	print '(a10,10F9.5)', 'Joutr',cntzo(1,2,1,1:10)
	print '(11(1pe10.2))',dhatz(1,1,1:11) 
	print '(11(1pe10.2))',dhatz(2,1,1:11) 
#endif

#ifdef UNDER_RELAX
      dhat = (1.-underrelx)*dhatd + underrelx*dhat
      betaphis = (1.-underrelx)*betaphisd + underrelx*betaphis
      if( nz .gt. 1 ) then
          dhatz= (1.-underrelx)*dhatzd+ underrelx*dhatz
          betaphisz= (1.-underrelx)*betaphiszd+ underrelx*betaphisz
      endif
#endif

      return
      end
