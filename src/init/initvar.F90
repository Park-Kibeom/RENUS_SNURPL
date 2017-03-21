    subroutine initvar
!
      use param
!
      include 'global.h'
      include 'xsec.h'
      include 'ffdm.h'
      include 'itrcntl.h'
      include 'srchppm.h'
      include 'pow.h'   
!
      mgb(1)=1
      mge(1)=mgb(2)-1
      if(mge(1).eq.0) mge(1)=mgb(1)
      mge(2)=ng
      
      plevel0=plevel    
      plevel00 = plevel0   ! 2014_09_15 . scb  
      
! initialize flux
! 2013_10_02 . scb
      !do k=1,nz
      !  do l=1,nxy
      !    do m=1,ng
      !         !phif(m,l,k)=xschif(m,l,k)
      !      phif(m,l,k)=1.d0   ! 2013_10_02 . scb   
      !    enddo
      !  enddo
      !enddo
      do k=1,nz
        do l=1,nxy
            do m2=1,2
                phif(mgb(m2):mge(m2),l,k)=1.0/(mge(m2)-mgb(m2)+1)
                do m=mgb(m2),mge(m2)
                    phi(m2,l,k)=phi(m2,l,k)+phif(m,l,k) 
                enddo
            enddo
        enddo
      enddo
! added end
!
      return
    end subroutine
