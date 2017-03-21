subroutine facilumg
    use const
    use bicgmg
    use geom, only : ng, nxs, nxe, nodel
    use cmfdmg, only : diagf, ccrf
    implicit none
    
    integer                 :: m, i, j, k, l, ln, jm1, lnm1, lnp1
    
! loop over groups
    do m=1,ng
    !   loop over planes
        do k=1,nz
            ! first row
            j=1
            do i=nxs(j),nxe(j)
                l=nodel(i,j)
                del1g(i,m)=diagf(l,k,m)
                al1g(i,k,m)=ccrf(1,l,k,m)
                au1g(i,m)=ccrf(2,l,k,m)
            enddo
    !       === loop over rows ===
            do j=2,ny
                jm1=j-1
                ! obtain incomplete lu factor for the 1d matrix of the row
                call facilu1d1g(m,jm1,k)
                ! obtain the inverse of the 1d matrix
                call abi1d1g(m,jm1,k)

                ! d_j+1 = a_j+1 - l_j+1*jnv(d_j)*u_j
                do i=nxs(j),nxe(j)
                    l=nodel(i,j)

                    if (i.ge.nxs(jm1) .and. i.le.nxe(jm1)) then
                        ln=nodel(i,jm1)
                        del1g(i,m)=diagf(l,k,m)-ccrf(3,l,k,m)*ainvd1g(i,m)*ccrf(4,ln,k,m)
                    else
                        del1g(i,m)=diagf(l,k,m)
                    endif

                    if(i.ne.nxs(j) .and. (i-1).ge.nxs(jm1) .and. (i-1).le.nxe(jm1)) then
                        lnm1=nodel(i-1,jm1)
                        al1g(l,k,m)=ccrf(1,l,k,m)-ccrf(3,l,k,m)*ainvl1g(i,m)*ccrf(4,lnm1,k,m)
                    else
                        al1g(l,k,m)=ccrf(1,l,k,m)
                    endif
                    
                    if(i.ne.nxe(j) .and. (i+1).ge.nxs(jm1) .and. (i+1).le.nxe(jm1)) then
                        lnp1=nodel(i+1,jm1)
                        au1g(i,m)=ccrf(2,l,k,m)-ccrf(3,l,k,m)*ainvu1g(i,m)*ccrf(4,lnp1,k,m)
                    else
                        au1g(i,m)=ccrf(2,l,k,m)
                    endif
                enddo
            enddo
    !       obtain incomplete lu factor for the 1d matrix of the last row
            call facilu1d1g(m,ny,k)
        enddo
    enddo
    return      
end subroutine

