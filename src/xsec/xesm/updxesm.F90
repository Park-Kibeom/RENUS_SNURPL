      subroutine updxesm
!
! Update previous time step number densities for Iodine/Xenon and Promethium/Samarium.
      include 'global.h'
      include 'xesm.h'
      include 'geom.h'
!
! update Iodine/Xenon number densities and Promethium/Samarium
! number densities
      do k=kfbeg,kfend
         do l=1,nxy
            rnip(l,k) =rni(l,k)
            rnxep(l,k)=rnxe(l,k)
            rnpmp(l,k)=rnpm(l,k)
            rnsmp(l,k)=rnsm(l,k)
         enddo
      enddo
      !! 2014_09_16 . scb for dbg
      !print *, rnip(18,2),rnxep(18,2),rnpmp(18,2),rnsmp(18,2)
      !! added end
      
      return
      
      end
