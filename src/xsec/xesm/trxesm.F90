!      subroutine trxesm
!!
!! Update transient number densities for Iodine/Xenon and Promethium/Samarium.
!      use tran,   only : deltm
!
!      include 'global.h'
!      include 'files.h'
!      include 'xesm.h'
!      include 'geom.h'
!      include 'pow.h'
!      include 'thgeom.inc'
!      include 'ffdm.h'
!      include 'xsec.h'
!!
!      if(ixesmopt.eq.0) return   ! no xe/sm
!      
!      write(*,*) 'Xenon/Samarium Update...'
!      write(io8,*) 'Xenon/Samarium Update...'
!      
!! update node-wise microscopic absorption cross sections for Xenon and Samarium
!      call xsafp
!      
!! update Iodine/Xenon number densities and Promethium/Samarium number densities
!      vol = 0.d0
!      xeave = 0.d0  
!      !flxlevel=1.d0  ! 2013_10_04 . scb for dbg       
!!
!      if(ixesmopt.eq.1) then     !equilibrium xe/sm
!          do k=kfbeg,kfend
!            do l=1,nxy
!              la=ltola(l)
!              iat=iassytyp(la)
!              if(.not.iffuela(iat)) cycle
!              
!              psix=0.d0
!              abxe=0.d0
!              absm=0.d0
!              do m=1,ng
!                phix=flxlevel*phif(m,l,k)
!                psix=psix + xsff(m,l,k)*phix
!                abxe=abxe + xsxeaf(m,l,k)*phix
!                absm=absm + xssmaf(m,l,k)*phix
!              enddo
!              
!              rni(l,k) =(gami(l,k)*psix)/rlambi
!              rnxe(l,k)=(gami(l,k)+gamxe(l,k))*psix/(rlambxe+abxe)
!              rnpm(l,k)=(gampm(l,k)*psix)/rlambpm
!
!              if (absm .le. 1.0e-30) then
!                rnsm(l,k) = 0.0
!              else
!                rnsm(l,k) = (gampm(l,k)*psix)/absm
!              endif
!
!!  total xenon number density
!              vol = vol + volnode(l,k)
!              xeave = xeave + rnxe(l,k)*volnode(l,k)  
!            enddo
!         enddo
!
!!  average xenon number density
!         xeave = xeave/vol
!       elseif(ixesmopt.eq.2 .or. ixesmopt.eq.3) then   ! transient xe/sm
!         deltsec=deltm
!         do k=kfbeg,kfend
!            do l=1,nxy
!              la=ltola(l)
!              iat=iassytyp(la)
!              if(.not.iffuela(iat)) cycle
!              
!              psix=0.d0
!              abxe=0.d0
!              absm=0.d0
!              do m=1,ng
!                phix=flxlevel*phif(m,l,k)
!                psix=psix + xsff(m,l,k)*phix
!                abxe=abxe + xsxeaf(m,l,k)*phix
!                absm=absm + xssmaf(m,l,k)*phix
!              enddo
!              
!              expidt=exp(-rlambi*deltsec)
!              rni(l,k) =rnip(l,k)*expidt+(gami(l,k)*psix)*(1-expidt)/rlambi
!              dexpxe=rlambxe+abxe
!              expxedt=exp(-dexpxe*deltsec)
!              rnxe(l,k)=rnxep(l,k)*expxedt + (gami(l,k)+gamxe(l,k))*psix*(1-expxedt)/dexpxe  &
!                      +(rlambi*rnip(l,k)-gami(l,k)*psix)*(expidt-expxedt)/(dexpxe-rlambi)
!              exppmdt=exp(-rlambpm*deltsec)
!              expsmdt=exp(-absm*deltsec)
!              rnpm(l,k)=rnpmp(l,k)*exppmdt+(gampm(l,k)*psix)*(1-exppmdt)/rlambpm
!              rnsm(l,k)=rnsmp(l,k)*expsmdt+gampm(l,k)*psix*(1-expsmdt)/absm     &
!                      +(rlambpm*rnpmp(l,k)-gampm(l,k)*psix)*(exppmdt-expsmdt)/(absm-rlambpm)
!!  total xenon number density
!              vol = vol + volnode(l,k)
!              xeave = xeave + rnxe(l,k)*volnode(l,k)  
!            enddo
!         enddo
!!  average xenon number density
!         xeave = xeave/vol         
!      endif
!!
!! reset sm
!      if(ismopt.eq.0) then
!         do k=kfbeg,kfend
!            do l=1,nxy
!               rnpm(l,k)=0.0
!               rnsm(l,k)=0.0
!            enddo
!         enddo
!      endif
!
!      return
!      end
!
!