! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    subroutine setlshex2g(iftran,iroute)

! Set CMFD Linear System for Hexagonal Geometry

      use const
      use geomhex
      use tpen
      use xsec
      use cmfdhex2g, only : dcmat, cmat, dsum, dfd, dhat, dfdz, dhatz
      use tran,      only : betat, betap
      use tranxsec,  only : xschid

	    use cmfd2g, only : reigvs, af
	    use sfam,   only : psi, phi

	    implicit none

	    integer :: iroute
	    logical :: iftran

      real :: tmpm(6)
	    integer :: iz, ih, it, isfc, ig, ir2, ifd, ir, ifr
	    integer :: k, m, l, i
	    real :: vol, chi, chip, chid, chi1, chi2
      
! 2013_10_11 . scb       for dbg
      !integer,save :: iter=0
      !
      !iter=iter+1
      !
      !if (iter.eq.22) then
      !  continue
      !endif
      
! added end      

      if(iroute.eq.2) then
        do iz=1,nz
          do ih=1,nassy
            vol=volnode(ih,iz)
            chip=reigvs
            chid=0.
            chi1=chip*xschi(1,ih,iz)
            chi2=chip*xschi(2,ih,iz)
            dcmat(1,ih,iz)=(xst(1,ih,iz)-chi1*xsnf(1,ih,iz) &
                           +dsum(1,ih,iz))*vol
            dcmat(2,ih,iz)=-chi1*xsnf(2,ih,iz)*vol-xss(2,1,ih,iz)*vol ! mod
!            dcmat(2,ih,iz)=-chi1*xsnf(2,ih,iz)*vol
            dcmat(3,ih,iz)=-(xss(1,2,ih,iz)+chi2*xsnf(1,ih,iz))*vol
            dcmat(4,ih,iz)=(xst(2,ih,iz)-chi2*xsnf(2,ih,iz) &
                           +dsum(2,ih,iz))*vol
          enddo
        enddo
        return
      endif

      do iz=1,nz
        do ih=1,nassy
          do it=1,6
            isfc=neigsfc(it,ih)
            do ig=1,2
              cmat(ig,it,ih,iz)=r9hs2 &
                *(dfd(ig,isfc,iz)+wtdhat(it,ih)*dhat(ig,isfc,iz))
            enddo
          enddo
          if(nz.eq.1) go to 100 
          isfc=neigsfcz(1,iz)
          do ig=1,ng2
            cmat(ig,7,ih,iz)=rhzbar2(1,iz) &
                      *(dfdz(ig,ih,isfc)-dhatz(ig,ih,iz))
          enddo
          isfc=neigsfcz(2,iz)
          do ig=1,ng2
            cmat(ig,8,ih,iz)=rhzbar2(2,iz) &
                      *(dfdz(ig,ih,isfc)+dhatz(ig,ih,iz+1))
          enddo
  100     continue
        enddo
      enddo

! Dsum Calculation
      do iz=1,nz
        do ih=1,nassy
          do ig=1,2
            dsum(ig,ih,iz)=0 
            do it=1,6
              isfc=neigsfc(it,ih)
              dsum(ig,ih,iz)=dsum(ig,ih,iz) &
                  +dfd(ig,isfc,iz)-wtdhat(it,ih)*dhat(ig,isfc,iz)
            enddo
            dsum(ig,ih,iz)=dsum(ig,ih,iz)*r9hs2
            if(nz.gt.1) then 
              isfc=neigsfcz(1,iz)
              dsum(ig,ih,iz)=dsum(ig,ih,iz) &
                 +rhzbar2(1,iz)*(dfdz(ig,ih,isfc)+dhatz(ig,ih,iz))
              isfc=neigsfcz(2,iz)
              dsum(ig,ih,iz)=dsum(ig,ih,iz) &
                +rhzbar2(2,iz)*(dfdz(ig,ih,isfc)-dhatz(ig,ih,iz+1))
            endif
          enddo               
        enddo
      enddo
! condense cmat on radial direction
      do iz=1,nz
        do ih=1,nassy
          do ir2=1,ineigcond(0,0,ih)
            ifr=ineigcond(ir2,0,ih)
            dsum(1,ih,iz)=dsum(1,ih,iz)-cmat(1,ifr,ih,iz)
            dsum(2,ih,iz)=dsum(2,ih,iz)-cmat(2,ifr,ih,iz)
          enddo
          do ig=1,2
            do ir=1,ipntr(0,ih)
              tmpm(ir)=0
              do ir2=1,ineigcond(0,ir,ih)
                ifr=ineigcond(ir2,ir,ih)
                tmpm(ir)=tmpm(ir)+cmat(ig,ifr,ih,iz)
              enddo
            enddo
            do ir=1,ipntr(0,ih)
              cmat(ig,ir,ih,iz)=tmpm(ir)
            enddo
          enddo
        enddo
      enddo

	if(iftran) then
      do iz=1,nz
        do ih=1,nassy
          chip=(1-betat(ih,iz))*reigvs
          chid=(betap(ih,iz)+betat(ih,iz)-1)*reigvs
          chi1=chip*xschi(1,ih,iz)+chid*xschid(1,ih,iz)
          chi2=chip*xschi(2,ih,iz)+chid*xschid(2,ih,iz)
          dcmat(1,ih,iz)=xst(1,ih,iz) &
                        -chi1*xsnf(1,ih,iz)+dsum(1,ih,iz)
          dcmat(2,ih,iz)=-chi1*xsnf(2,ih,iz)-xss(2,1,ih,iz) ! mod
!          dcmat(2,ih,iz)=-chi1*xsnf(2,ih,iz)
          dcmat(3,ih,iz)=-chi2*xsnf(1,ih,iz)-xss(1,2,ih,iz)
          dcmat(4,ih,iz)=xst(2,ih,iz)-chi2*xsnf(2,ih,iz) &
                        +dsum(2,ih,iz)
        enddo
      enddo
	else
      do iz=1,nz
        do ih=1,nassy
          chip=reigvs
          chid=0.
          chi1=chip*xschi(1,ih,iz)
          chi2=chip*xschi(2,ih,iz)
          dcmat(1,ih,iz)=xst(1,ih,iz) &
                        -chi1*xsnf(1,ih,iz)+dsum(1,ih,iz)
          dcmat(2,ih,iz)=-chi1*xsnf(2,ih,iz)-xss(2,1,ih,iz) ! mod
!          dcmat(2,ih,iz)=-chi1*xsnf(2,ih,iz)
          dcmat(3,ih,iz)=-chi2*xsnf(1,ih,iz)-xss(1,2,ih,iz)
          dcmat(4,ih,iz)=xst(2,ih,iz)-chi2*xsnf(2,ih,iz) &
                        +dsum(2,ih,iz)
        enddo
      enddo
	endif



!	call dmalloc(af,2,nxy,nz)
      do k=1,nz
         do l=1,nxy
            vol=volnode(l,k)  
!            if(ifupdpinv) then
!               pinvd=pinv(m,l,k)
!               pinv(m,l,k)=log(abs(phi(m,l,k)/phid(m,l,k)))/delt
!               rvdelt(m,l,k)=1/(velo(m,l,k)*cetak*delt)
!     +                      +pinv(m,l,k)/velo(m,l,k)
!               expo(m,l,k)=abs(phi(m,l,k)/phid(m,l,k))
!               src(m,l,k)=expo(m,l,k)*(srcfac1(m,l,k)
!     +             -cetakr*pinv(m,l,k)/velo(m,l,k)*vol*phid(m,l,k))
!     +             +srcfac2(m,l,k)
!            endif
            af(1,l,k)=xsnf(1,l,k)*vol
            af(2,l,k)=xsnf(2,l,k)*vol
!            psi(l,k)=af_hex(1,l,k)*phi(1,l,k)+af_hex(2,l,k)*phi(2,l,k)
            do i=1,4
               dcmat(i,l,k)=dcmat(i,l,k)*vol
            enddo
! scb . dcmat for 1 and 4 is different compared to PARCS . But this is covered in drivecmfd2g routine
            do it=1,ipntr(0,l)
               do m=1,ng2
                  cmat(m,it,l,k)=cmat(m,it,l,k)*vol
               enddo
            enddo
            do it=7,8
               do m=1,ng2
                  cmat(m,it,l,k)=cmat(m,it,l,k)*vol
               enddo
            enddo
         enddo
      enddo
      do l=1,nxy
         do m=1,ng2
            cmat(m,7,l,1)=0
            cmat(m,8,l,nz)=0
         enddo
      enddo

      return
      end
