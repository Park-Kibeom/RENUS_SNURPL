! added in ARTOS ver. 0.2 . 2012_07_24 by SCB      
   subroutine axbhex(b,ab)

! Matrix-vector product for tpen

	   use geomhex
	   use cmfdhex2g
	   implicit none

      real,pointer :: b(:,:,:)
      real,pointer :: ab(:,:,:)

	   integer :: k, kp1, km1, l, is, ln

      do k=1,nz
         kp1=k+1
         km1=k-1
         if(kp1.gt.nz) kp1=nz
         if(km1.lt.1) km1=1
         do l=1,nxy
            ab(1,l,k)=dcmat(1,l,k)*b(1,l,k)+dcmat(2,l,k)*b(2,l,k) &
                     -cmat(1,7,l,k)*b(1,l,km1)-cmat(1,8,l,k)*b(1,l,kp1)
            ab(2,l,k)=dcmat(3,l,k)*b(1,l,k)+dcmat(4,l,k)*b(2,l,k) &
                     -cmat(2,7,l,k)*b(2,l,km1)-cmat(2,8,l,k)*b(2,l,kp1)
            do is=1,ipntr(0,l)
               ln=ipntr(is,l)
               ab(1,l,k)=ab(1,l,k)-cmat(1,is,l,k)*b(1,ln,k)
               ab(2,l,k)=ab(2,l,k)-cmat(2,is,l,k)*b(2,ln,k)
            enddo            
         enddo
      enddo

      return
   end subroutine

