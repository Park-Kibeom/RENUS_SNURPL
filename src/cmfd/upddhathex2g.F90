! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    subroutine upddhathex2g      

! Update CMFD coefficient(dhat and beta) from Nodal Solution 

	    use const
	    use geomhex
	    use tpen
	    use cmfdhex2g
	    use xsec
	    use sfam, only : phif, phi

	    implicit none

!dj+ : add under-realaxation of dhat and betaphis
!#ifdef DBL
      double precision underrelx
      data underrelx /1.0/
      
      integer :: k,l  ! 2012_10_23 . scb
!#else
!      real underrelx
!      data underrelx /0.8/
!#endif
!dj-


! variables for nem solver in z-direction
! atleak : transverse leakage for axial direction
! atleaku,atleakl : transverse leakage for upper and lower node
! dcu,dcl : xsd(diffusion coefficients) for upper and lower node

      integer :: mvz(2)
      data mvz/2,1/

	    integer :: isfc, ind, is1, is2, id1, id2, iz, &
	             ig, iz1, iz2, ih
	    real :: cnto1, cnto2, flxs, cnts, flx1, flx2, dfdm, dc1,  &
	          w, fdflxs, adfl, reflratl, hzr
      real :: reflratz

#ifdef UNDEr_RELAX
      dhatd = dhat
      dhatzd= dhatz
      betaphisd = betaphis
	    betaphiszd= betaphisz
#endif
 
      if(ng.ne.2) then
         call upddhathexmg
         goto 999
      endif

! update D(diffusion coefficient) for CMFD
      do isfc=1,nxsfc
        ind=neigsnd(3,isfc)
        if(ind.eq.12) then
          is1=neigsnd(1,isfc)
          is2=neigsnd(2,isfc)
          id1=neigsnd(4,isfc)
          id2=neigsnd(5,isfc)
          do iz=1,nz 
            do ig=1,ng2
              cnto1=cnto(ig,id1,is1,iz)
              cnto2=cnto(ig,id2,is2,iz)
              flxs=2.*(cnto1+cnto2)
              cnts=cnto1-cnto2
              flx1=hflx(ig,is1,iz)  
              flx2=hflx(ig,is2,iz)  
              dfdm=dfd(ig,isfc,iz)
              dc1=xsd(ig,is1,iz)
              w=dfdm/dc1/2.
              fdflxs=(1.-w)*flx1+w*flx2 
              dhat(ig,isfc,iz)=  &
                 -(rt3*hside*cnts+dfdm*(flx2-flx1))/(flx2+flx1)
              betaphis(ig,isfc,iz)=2.*(flxs-fdflxs)/(flx1+flx2)
            enddo
          enddo
        else
          is1=neigsnd(ind,isfc)
          id1=neigsnd(ind+3,isfc)
          do iz=1,nz 
            do ig=1,ng2
              cnto1=cnto(ig,id1,is1,iz)
!              adfl=xsadf(ig,is1,iz)
              adfl=1.
              reflratl=(adfl-2*alxr)/(adfl+2*alxr)
              flxs=2.*cnto1*(1.+reflratl)
              cnts=cnto1*(1.-reflratl)
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
              cnto1=cntzo(ig,id1,ih,iz1)
              cnto2=cntzo(ig,id2,ih,iz2)
              flxs=2.*(cnto1+cnto2)
              cnts=cnto1-cnto2
              flx1=hflx(ig,ih,iz1)  
              flx2=hflx(ig,ih,iz2) 
              dfdm=dfdz(ig,ih,iz)
              dc1=xsd(ig,ih,iz1)
              w=dfdm/dc1/(1.+hzr)
              fdflxs=(1.-w)*flx1+w*flx2 
              dhatz(ig,ih,iz)=-((hz(iz1)+hz(iz2))*cnts/2. &
                        +dfdm*(flx2-flx1))/(flx2+flx1)
              betaphisz(ig,ih,iz)=2.*(flxs-fdflxs)/(flx1+flx2)
            enddo
          enddo
        else
          iz1=neigsndz(ind,iz)
          id1=neigsndz(ind+3,iz)
          if(iz.eq.1) then                                                              
            reflratz=reflratzb
          else
            reflratz=reflratzt
          endif
          do ih=1,nassy 
            do ig=1,ng2
              cnto1=cntzo(ig,id1,ih,iz1)
              flxs=2.*(1.+reflratz)*cnto1
              cnts=(1.-reflratz)*cnto1
              flx1=hflx(ig,ih,iz1)  
              dfdm=dfdz(ig,ih,iz)
              dhatz(ig,ih,iz)=(dfdm*flx1-hz(iz1)*cnts)/flx1
              if(ind.eq.2) dhatz(ig,ih,iz)=-dhatz(ig,ih,iz)
              betaphisz(ig,ih,iz)=flxs/flx1
            enddo
          enddo
        endif
      enddo

! update old abs xsec
#ifdef COND_NODAL         ! commented by 2012_10_23 . scb
      do k=1,nz
         do l=1,nxy
            xstd(1,l,k)=xst(1,l,k)
            xstd(2,l,k)=xst(2,l,k)
!            phi(1,l,k)=hflx(1,l,k)        
!            phi(2,l,k)=hflx(2,l,k)         
         enddo
      enddo
#endif          ! commented by 2012_10_23 . scb

#ifdef DEBUG
	    print *
	    do ih=1,nxy
	    print '(a10,10F9.5)', 'Flux' ,phif(1,ih,1:nz) 
	    print '(a10,10F9.5)', '    ' ,phif(2,ih,1:nz) 
	    enddo
	    print '(/a10,10F9.5)', 'Joutl',cntzo(1,1,1,1:nz)
	    print '(a10,10F9.5)', 'Joutr',cntzo(1,2,1,1:nz)
	    print '(11(1pe10.2))',dhatz(1,1,1:nz+1) 
	    print '(11(1pe10.2))',dhatz(2,1,1:nz+1) 
#endif


  999 continue

#ifdef UNDER_RELAX
      dhat = (1.-underrelx)*dhatd + underrelx*dhat
      betaphis = (1.-underrelx)*betaphisd + underrelx*betaphis
      if( nz .gt. 1 ) then
          dhatz= (1.-underrelx)*dhatzd+ underrelx*dhatz 
          betaphisz= (1.-underrelx)*betaphiszd+ underrelx*betaphisz
      endif
#endif

 1000 return
    end
