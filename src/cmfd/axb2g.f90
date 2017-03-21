subroutine axb2g(phi,aphi)
! calculate A x phi and store the result to aphi
! this is used for given 2-group
    use const
    use cmfd2g, only : am, ccr, ccz,nx,ny,nz,nxy,indm24
    use geom,   only : neibr,neibz
    implicit none
    
    real,pointer            :: phi(:,:,:)
    real,pointer            :: aphi(:,:,:)
    integer                 :: k, l, m, m1, m2, idir, idirz
    real                    :: aphil


    do k=1,nz
     do l=1,nxy
        do m=1,ng2  ! 1 to 2
          m1=indm24(m,1)
          m2=indm24(m,2)
!am(:,l,k) is 2 x 2 matrix and phi(:,l,k) is 2 x 1.
          aphil=am(m1,l,k)*phi(1,l,k)+am(m2,l,k)*phi(2,l,k)

          do idir=1,nrdir2
            aphil=aphil+ccr(idir,m,l,k)*phi(m,neibr(idir,l),k)
          enddo               
          do idirz=1,2 ! bottom and up
            aphil=aphil+ccz(idirz,m,l,k)*phi(m,l,neibz(idirz,k))

          enddo
           aphi(m,l,k)=aphil
        enddo
     enddo
    enddo

    return 
end subroutine