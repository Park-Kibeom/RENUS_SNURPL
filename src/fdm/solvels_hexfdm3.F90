! added in ARTOS ver. 0.2 . 2012_07_24 by SCB     
   subroutine solvels_hexfdm3(epseig,epsl2)

! Solve Linear System for Hexagonal Geometry
	   use const
	   use geomhex
	   use hexfdm3
      use sfam_cntl,  only : nupmax,maxupscatgr
      use sfam,       only : phif, psi, reigv, eigv
      use xsec,       only : xschif,xssf,xssfs,xssfe,xsnff
      use fdm_solver
	   implicit none

      real :: epseig,epsl2
! local
      integer :: iout, k, l, n, m
      real :: errl2, erreig
      character :: mesg*120 

! initialization
      do k=1,nz
      do l=1,nassy
         do n=1,nsub
            psifdm(n,l,k)=1.
            do m=1,ng
               phifdm(m,n,l,k)=1.
               phifdm2(m,n,l,k)=0.1
            enddo
         enddo
         psi(l,k)=0.
         do m=1,ng
            phif(m,l,k)=0.
         enddo
      enddo
      enddo
      reigv=1.
      call malloc_fdmsolver

      do iout=1,50000
         call sol_gs3(reigv,psifdm,phifdm,phifdm2)
         call sol_power3(reigv,phifdm,phifdm2,errl2,erreig,psifdm)
    
#ifndef CHEBY
         call sol_cheby(phifdm,psifdm,iout)
         write(mesg, '(i10,f15.7,4(1pe15.5),i5)'), iout, eigv, erreig, errl2, domr, erfratio, mcp     
         call message(true,true,mesg)
#else
         print '(i10,f15.7,2(1pe15.5))', iout, eigv, erreig, errl2
#endif
         eigv=1./reigv               
         if(errl2.lt.epsl2 .and. erreig.lt.epseig) exit 
      enddo ! iout

      do k=1,nz
      do l=1,nassy
         do n=1,nsub
            do m=1,ng
               phif(m,l,k)=phif(m,l,k)+(phifdm(m,n,l,k)-2.*phifdm2(m,n,l,k)) 
            enddo
            psi(l,k)=psi(l,k)+psifdm(n,l,k)
         enddo
      enddo
      enddo

      return
   end subroutine
