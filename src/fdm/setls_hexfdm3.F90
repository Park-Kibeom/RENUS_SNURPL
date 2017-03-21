! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine setls_hexfdm3

! Set FDM Linear System for Hexagonal Geometry

	   use const
	   use geomhex
	   use hexfdm3
	   use xsec
	   implicit none

      integer :: m,l,k,n,idir,isfc
      integer :: ibeg, iend
      real :: hdiv2, r4hdiv2
      real,save :: r4o3=4./3., r5o3=5./3., r2o3=2./3.
      real :: sigtrh
      real :: bc(4), rt3rh

      real :: offdiag, offdiag2       
      real :: offd_ri,  offd_ro,  offd_zb,  offd_zt 
      real :: offd_ri2, offd_ro2, offd_zb2, offd_zt2
        
      rt3rh=rt3/hdiv
      do m=1,ng
         do k=1,nz
            do l=1,nassy  
               offd_ri=rt3*xsdf(m,l,k)*hz(k)  
               offd_ri2=rt3*xsdf2(m,l,k)*hz(k)  
  
               offd_zb=dtilzfdm(m,l,neigsfcz(1,k))*triarea
               offd_zt=dtilzfdm(m,l,neigsfcz(2,k))*triarea 
                             
               offd_zb2=dtilzfdm2(m,l,neigsfcz(1,k))*triarea
               offd_zt2=dtilzfdm2(m,l,neigsfcz(2,k))*triarea                 
                 
               sigtrh=r4o3*xstf(m,l,k)+r5o3*xstrf(m,l,k)
               ! boundary triangles
               do idir=1,6
                  offd_ro=dtilrfdm(m,neigsfc(idir,l),k)*hdiv*hz(k)
                  offd_ro2=dtilrfdm2(m,neigsfc(idir,l),k)*hdiv*hz(k)

                  ibeg=(idir-1)*ndiv+1
                  iend=idir*ndiv
                  do n=ibeg,iend  
                     ! 0-th
                     ccrfdm(1,m,n,l,k)=-offd_ro
                     ccrfdm(2,m,n,l,k)=-offd_ri
                     ccrfdm(3,m,n,l,k)=-offd_ri
                     cczfdm(1,m,n,l,k)=-offd_zb
                     cczfdm(2,m,n,l,k)=-offd_zt
                     ! 2-nd
                     ccrfdm2(1,m,n,l,k)=-offd_ro2
                     ccrfdm2(2,m,n,l,k)=-offd_ri2
                     ccrfdm2(3,m,n,l,k)=-offd_ri2
                     cczfdm2(1,m,n,l,k)=-offd_zb2
                     cczfdm2(2,m,n,l,k)=-offd_zt2
                     ! sum of off-dialgonal terms
                     offdiag=offd_ro+2.*offd_ri+offd_zb+offd_zt
                     offdiag2=offd_ro2+2.*offd_ri2+offd_zb2+offd_zt2

                     ! diagonal matrix
                     diagfdm(1,m,n,l,k)=xstf(m,l,k)*volnodet(n,l,k) &
                                          +offdiag
                     diagfdm(4,m,n,l,k)=sigtrh*volnodet(n,l,k) &
                                          +offdiag2
                     diagfdm(2,m,n,l,k)=-2.*xstf(m,l,k)*volnodet(n,l,k)
                     diagfdm(3,m,n,l,k)=-r2o3*xstf(m,l,k)*volnodet(n,l,k)  
                     
                     ! radial bc                           
                     if(neignd(idir,l).eq.0 .and. alxr.eq.0.5) then                
                        call vacuum_bc(xsdf(m,l,k),xsdf2(m,l,k),rt3rh,bc)
                        diagfdm(1,m,n,l,k)=diagfdm(1,m,n,l,k)+bc(1)*hdiv*hz(k)
                        diagfdm(2,m,n,l,k)=diagfdm(2,m,n,l,k)+bc(2)*hdiv*hz(k)
                        diagfdm(3,m,n,l,k)=diagfdm(3,m,n,l,k)+bc(3)*hdiv*hz(k)
                        diagfdm(4,m,n,l,k)=diagfdm(4,m,n,l,k)+bc(4)*hdiv*hz(k)                                                    
                     endif  
                  enddo       
               enddo

               ! inner triangles
               do n=ndiv*6+1, nsub
                  ! 0-th
                  ccrfdm(1,m,n,l,k)=-offd_ri
                  ccrfdm(2,m,n,l,k)=-offd_ri
                  ccrfdm(3,m,n,l,k)=-offd_ri 
                  cczfdm(1,m,n,l,k)=-offd_zb
                  cczfdm(2,m,n,l,k)=-offd_zt
                  ! 2-nd
                  ccrfdm2(1,m,n,l,k)=-offd_ri2
                  ccrfdm2(2,m,n,l,k)=-offd_ri2
                  ccrfdm2(3,m,n,l,k)=-offd_ri2
                  cczfdm2(1,m,n,l,k)=-offd_zb2
                  cczfdm2(2,m,n,l,k)=-offd_zt2
                  ! sum of off-dialgonal terms
                  offdiag=3.*offd_ri+offd_zb+offd_zt
                  offdiag2=3.*offd_ri2+offd_zb2+offd_zt2

                  ! diagonal matrix
                  diagfdm(1,m,n,l,k)=xstf(m,l,k)*volnodet(n,l,k) &
                                       +offdiag
                  diagfdm(4,m,n,l,k)=sigtrh*volnodet(n,l,k) &
                                       +offdiag2
                  diagfdm(2,m,n,l,k)=-2.*xstf(m,l,k)*volnodet(n,l,k)
                  diagfdm(3,m,n,l,k)=-r2o3*xstf(m,l,k)*volnodet(n,l,k)         
               enddo

               ! bottom bc
               if(neigz(1,k).eq.0 .and. alzl.eq.0.5) then 
                  do n=1,nsub  
                     call vacuum_bc(xsdf(m,l,k),xsdf2(m,l,k),1./hz(k),bc)
                     diagfdm(1,m,n,l,k)=diagfdm(1,m,n,l,k)+bc(1)*triarea
                     diagfdm(2,m,n,l,k)=diagfdm(2,m,n,l,k)+bc(2)*triarea
                     diagfdm(3,m,n,l,k)=diagfdm(3,m,n,l,k)+bc(3)*triarea
                     diagfdm(4,m,n,l,k)=diagfdm(4,m,n,l,k)+bc(4)*triarea
                  enddo
               endif
               
               ! top bc
               if(neigz(2,k).eq.0 .and. alzr.eq.0.5) then 
                  do n=1,nsub  
                     call vacuum_bc(xsdf(m,l,k),xsdf2(m,l,k),1./hz(k),bc)
                     diagfdm(1,m,n,l,k)=diagfdm(1,m,n,l,k)+bc(1)*triarea
                     diagfdm(2,m,n,l,k)=diagfdm(2,m,n,l,k)+bc(2)*triarea
                     diagfdm(3,m,n,l,k)=diagfdm(3,m,n,l,k)+bc(3)*triarea
                     diagfdm(4,m,n,l,k)=diagfdm(4,m,n,l,k)+bc(4)*triarea
                  enddo
               endif

            enddo ! l
         enddo ! k
      enddo ! m

      return
      end subroutine


      subroutine vacuum_bc(dc,dc2,rh,bc)
      implicit none

      real :: dc, dc2, rh, bc(4)
      real :: a(4), b(4), c(4), t(4), invb(4), det

      a(1)=0.5
      a(2)=-0.375
      a(3)=-0.125
      a(4)=0.875

      c(1)=2.*rh*dc
      c(2)=0.
      c(3)=0.
      c(4)=2.*rh*dc2

      b=a+c
      det=1./(b(1)*b(4)-b(2)*b(3))

      invb(1)=det*b(4)
      invb(2)=-det*b(2)
      invb(3)=-det*b(3)
      invb(4)=det*b(1)

      t(1)=a(1)*invb(1)+a(2)*invb(3)
      t(2)=a(1)*invb(2)+a(2)*invb(4)
      t(3)=a(3)*invb(1)+a(4)*invb(3)
      t(4)=a(3)*invb(2)+a(4)*invb(4)

      bc(1)=t(1)*c(1)+t(2)*c(3)
      bc(2)=t(1)*c(2)+t(2)*c(4)
      bc(3)=t(3)*c(1)+t(4)*c(3)
      bc(4)=t(3)*c(2)+t(4)*c(4)

      return
      end subroutine
