!      subroutine ssxesm
!!
!! Compute steady-state number densities for Iodine/Xenon and Promethium/Samarium.
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
!      include 'thop.inc'
!
!! 2013_10_02 . scb      
!      !logical,save :: first=.true.
!      !
!      !if(first) then
!      !  open(unit=1002,file='xesm',status='unknown')
!      !endif
!! added end      
!      
!      if(ixesmopt.eq.0) return   ! no xe/sm      
!      
!      write(*,*) 'Xenon/Samarium Update...'
!      write(io8,*) 'Xenon/Samarium Update...'
!
!! update node-wise microscopic absorption cross sections for Xenon and Samarium
!      call xsafp
!      
!! compute absolute flux level
!      totpow=0.d0
!      do k=kfbeg,kfend
!        do l=1,nxy
!          la=ltola(l)
!          iat=iassytyp(la)
!          if(.not.iffuela(iat)) cycle
!
!          vol=volnode(l,k)
!          do m=1,ng
!            totpow=totpow+phif(m,l,k)*xskpf(m,l,k)*vol
!          enddo
!        enddo
!      enddo
!      
!      avgpow=totpow/volfuel
!      powfac=avgpow*hac*pfa*pfa*1.e6 !cubic m to cubic cm
!      if(hex) powfac=powfac*0.86602540378444 !txk: hex assy area correction
!      fnorm=(plevel*powfa)/powfac
!! form equilibrium Iodine/Xenon number densities and equilibrium
!! Promethium/Samarium number densities
!!
!      !fnorm=1.d0
!      !phix=1.d0
!      
!      vol = 0.d0
!      xeave = 0.d0  
!!
!      if(ixesmopt.eq.1 .or. ixesmopt.eq.3) then     !equilibrium xe/sm
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
!                phix=fnorm*phif(m,l,k)
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
!         !write(1002,*) xeave ! 2013_10_02 . scb      
!
!       elseif(ixesmopt.eq.2) then   ! transient xe/sm
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
!                phix=fnorm*phif(m,l,k)
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
!!     0:  no Xe/Sm for both steady-state and transient
!!     1:  equilibrium Xe/Sm for both steady-state and transient
!!     2:  transient Xe/Sm for both  steady-state and transient,
!!          if depl=true, the depletion time will be used for transient Xe/Sm at
!!     steadystate,
!!          otherwise, the time from T/H code will be used.
!!     3:  input Xe/Sm density will be used for steady state, these densities will
!!     also used for whole transient.
!!         if depl=true, the densities come from depletion restart file,
!!         otherwise, come from pbtt densities file.
!!     4:   equilibrium Xe/Sm for stready state and transient Xe/Sm for transient,
!!          steady-state densities will be used as initial densities for transient.
!