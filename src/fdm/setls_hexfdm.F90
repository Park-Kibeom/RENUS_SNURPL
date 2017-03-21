! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine setls_hexfdm

! Set FDM Linear System for Hexagonal Geometry

	   use const
	   use geomhex
	   use hexfdm
	   use xsec
	   implicit none

      integer :: m,l,k,n,idir,isfc
      integer :: ibeg, iend
      real :: hdiv2, r4hdiv2, offdiag
      real :: offd_ri, offd_ro
      real :: offd_zb, offd_zt   

      do m=1,ng
      do k=1,nz
      do l=1,nassy
         offd_ri=rt3*xsdf(m,l,k)*hz(k)

         offd_zb=dtilzfdm(m,l,neigsfcz(1,k))*triarea
         offd_zt=dtilzfdm(m,l,neigsfcz(2,k))*triarea
         
         ! boundary triangles
         do idir=1,6
            offd_ro=dtilrfdm(m,neigsfc(idir,l),k)*hdiv*hz(k)
            
            ibeg=(idir-1)*ndiv+1
            iend=idir*ndiv
            do n=ibeg,iend
               ccrfdm(1,m,n,l,k)=-offd_ro
               ccrfdm(2,m,n,l,k)=-offd_ri
               ccrfdm(3,m,n,l,k)=-offd_ri
               cczfdm(1,m,n,l,k)=-offd_zb
               cczfdm(2,m,n,l,k)=-offd_zt

               offdiag=offd_ro+2.*offd_ri+offd_zb+offd_zt
          
               diagfdm(m,n,l,k)=xstf(m,l,k)*volnodet(n,l,k) &
                                 +offdiag                                                             
            enddo
         enddo
         ! inner triangles                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
         do n=ndiv*6+1, nsub
            ccrfdm(1,m,n,l,k)=-offd_ri
            ccrfdm(2,m,n,l,k)=-offd_ri
            ccrfdm(3,m,n,l,k)=-offd_ri 
            cczfdm(1,m,n,l,k)=-offd_zb
            cczfdm(2,m,n,l,k)=-offd_zt

            offdiag=3.*offd_ri+offd_zb+offd_zt
            diagfdm(m,n,l,k)=xstf(m,l,k)*volnodet(n,l,k) &
                              +offdiag 
         enddo
      enddo
      enddo
      enddo

      return
      end subroutine
