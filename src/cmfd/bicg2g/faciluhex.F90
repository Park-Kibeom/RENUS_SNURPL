! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
      subroutine faciluhex
 
! perform incomplete LU factorization for the 3D coefficient matrices

	   use geomhex
	   use cmfdhex2g
	   implicit none


      real :: rinv(4),rmlt(4)
      integer :: ia1(4),ia2(4),ib1(4),ib2(4)
      data ia1/1,1,3,3/
      data ia2/2,2,4,4/
      data ib1/1,2,1,2/
      data ib2/3,4,3,4/

	   integer :: k, l, ig, ir, ifr, ir2, ipos
	   real :: rdet
     real :: det  ! 2013_10_10 . scb

      do k=1,nz
        do l=1,nxy
          do ir=1,7
            do ig=1,4
              cftm(ig,ir,l)=0
            enddo
          enddo
          do ir=1,ipntr(0,l)
            cftm(1,ir,l)=-cmat(1,ir,l,k)
            cftm(4,ir,l)=-cmat(2,ir,l,k)
          enddo
          do ig=1,4
            cftm(ig,7,l)=dcmat(ig,l,k)
          enddo
        enddo

        do l=2,nxy
          do ir=1,ilubnd(l)
            ifr=ipntr(ir,l)
! 2013_10_10 . scb            
            !det=abs(cftm(1,7,ifr)*cftm(4,7,ifr) &
            !       -cftm(2,7,ifr)*cftm(3,7,ifr))
            !
            !if(det.lt.1.e-30) then
            !  continue
            !endif
! added end            
            rdet=1/(cftm(1,7,ifr)*cftm(4,7,ifr) &
                   -cftm(2,7,ifr)*cftm(3,7,ifr))
            rinv(1)=cftm(4,7,ifr)*rdet
            rinv(2)=-cftm(2,7,ifr)*rdet
            rinv(3)=-cftm(3,7,ifr)*rdet
            rinv(4)=cftm(1,7,ifr)*rdet
            do ig=1,4
              rmlt(ig)=cftm(ia1(ig),ir,l)*rinv(ib1(ig)) &
                      +cftm(ia2(ig),ir,l)*rinv(ib2(ig))
            enddo
            do ig=1,4
              cftm(ig,ir,l)=rmlt(ig)
            enddo
            do ir2=ir+1,ipntr(0,l)
              ipos=iastopnt(ipntr(ir2,l),ifr)
              if(ipos.ne.0) then
                do ig=1,4
                  cftm(ig,ir2,l)=cftm(ig,ir2,l)                        &
                                -rmlt(ia1(ig))*cftm(ib1(ig),ipos,ifr)  &
                                -rmlt(ia2(ig))*cftm(ib2(ig),ipos,ifr)
                enddo
              endif
            enddo
            ipos=iastopnt(l,ifr)
            do ig=1,4
              cftm(ig,7,l)=cftm(ig,7,l)                          &
                          -rmlt(ia1(ig))*cftm(ib1(ig),ipos,ifr)  &
                          -rmlt(ia2(ig))*cftm(ib2(ig),ipos,ifr)
            enddo
          enddo
        enddo
        do l=1,nxy
          do ir=1,6
            xlufac(1,ir,l,k)=cftm(1,ir,l)
            xlufac(2,ir,l,k)=cftm(4,ir,l)
          enddo
! 2013_10_10 . scb     
          !det=abs(cftm(1,7,l)*cftm(4,7,l) &
          !       -cftm(2,7,l)*cftm(3,7,l))
          !  
          !if(det.lt.1.e-30) then
          !  continue
          !endif          
! added end
          rdet=1/(cftm(1,7,l)*cftm(4,7,l) &
                 -cftm(2,7,l)*cftm(3,7,l))
          delinv(1,l,k)=cftm(4,7,l)*rdet
          delinv(2,l,k)=-cftm(2,7,l)*rdet
          delinv(3,l,k)=-cftm(3,7,l)*rdet
          delinv(4,l,k)=cftm(1,7,l)*rdet
        enddo
      enddo

      return
      end
