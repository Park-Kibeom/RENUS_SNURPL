subroutine upddhat2g(phif, phi)
!
    use const
    use cmfd2g, only : dtilr,dtilz,dhatr,dhatz
    use geom,   only : nx,ny,nz,nxy,nzp1,ng,    &
                       nxs,nxe,nys,nye,         &
                       nodel,lsfc,neibr,neibz
    use xsec,   only : mgb,mge
    use cmfdmg, only : dtilrf,dhatrf,dtilzf,dhatzf
    implicit none
    
    real,pointer            :: phif(:,:,:)
    real,pointer            :: phi(:,:,:)

    integer                 :: i,j,k,l,ll,lr,m,ls,m2,ks
    real                    :: curmg, cur2g
    
    do k=1,nz
        ! x-dir
        do j=1,ny

            !d-hat at the left surface of each mesh.
            do i=nxs(j),nxe(j)
                l=nodel(i,j)
                ls=lsfc(WEST,l)
                ll=neibr(WEST,l)
                do m2=1,ng2
                    curmg = 0
                    do m=mgb(m2),mge(m2)
                        curmg = curmg - dtilrf(m,ls,k)*(phif(m,l,k)-phif(m,ll,k)) &
                                      -dhatrf(m,ls,k)*(phif(m,l,k)+phif(m,ll,k))

                    enddo
                    cur2g=-dtilr(m2,ls,k)*(phi(m2,l,k)-phi(m2,ll,k))
                    dhatr(m2,ls,k) = cur2g-curmg
                    dhatr(m2,ls,k) = dhatr(m2,ls,k) / (phi(m2,l,k)+phi(m2,ll,k))
                enddo 
            enddo

            !d-hat at the last surface
            i=nxe(j)
            l=nodel(i,j)
            ls=lsfc(EAST,l)
            lr = neibr(EAST,l)
            do m2=1,ng2
                if(lr.ne.0) then
                    curmg = 0
                    do m=mgb(m2),mge(m2)
                        curmg = curmg - dtilrf(m,ls,k)*(phif(m,lr,k)-phif(m,l,k)) &
                                      -dhatrf(m,ls,k)*(phif(m,lr,k)+phif(m,l,k))

                    enddo
                    cur2g=-dtilr(m2,ls,k)*(phi(m2,lr,k)-phi(m2,l,k))
                    dhatr(m2,ls,k) = cur2g-curmg
                    dhatr(m2,ls,k) = dhatr(m2,ls,k) / (phi(m2,lr,k)+phi(m2,l,k))
                else
                    curmg = 0
                    do m=mgb(m2),mge(m2)
                        curmg = curmg + dtilrf(m,ls,k)*(phif(m,l,k))    &
                                      -dhatrf(m,ls,k)*phif(m,l,k)
                    enddo
                    cur2g = dtilr(m2,ls,k)*(phi(m2,l,k))
                    dhatr(m2,ls,k) = cur2g - curmg
                    dhatr(m2,ls,k) = dhatr(m2,ls,k) / (phi(m2,l,k))
                endif
            enddo
        enddo

        ! y-dir
        do i=1,nx
            do j=nys(i),nye(i)
                l=nodel(i,j)
                ls=lsfc(NORTH,l)
                ll=neibr(NORTH,l)
                do m2=1,ng2
                    curmg = 0
                    do m=mgb(m2),mge(m2)
                    curmg = curmg - dtilrf(m,ls,k)*(phif(m,l,k)-phif(m,ll,k))   &
                                  - dhatrf(m,ls,k)*(phif(m,l,k)+phif(m,ll,k))
                    enddo
                    cur2g=-dtilr(m2,ls,k)*(phi(m2,l,k)-phi(m2,ll,k))
                    dhatr(m2,ls,k) = cur2g-curmg
                    dhatr(m2,ls,k) = dhatr(m2,ls,k) / (phi(m2,l,k)+phi(m2,ll,k))
                enddo 
            enddo
            
            j=nye(i)
            l=nodel(i,j)
            ls=lsfc(SOUTH,l)
            lr = neibr(SOUTH,l)
            do m2=1,ng2
                if(lr.ne.0) then
                    curmg = 0
                    do m=mgb(m2),mge(m2)
                        curmg = curmg - dtilrf(m,ls,k)*(phif(m,lr,k)-phif(m,l,k)) &
                                      -dhatrf(m,ls,k)*(phif(m,lr,k)+phif(m,l,k))

                    enddo
                    cur2g=-dtilr(m2,ls,k)*(phi(m2,lr,k)-phi(m2,l,k))
                    dhatr(m2,ls,k) = cur2g-curmg
                    dhatr(m2,ls,k) = dhatr(m2,ls,k) / (phi(m2,lr,k)+phi(m2,l,k))
                else
                    curmg = 0
                    do m=mgb(m2),mge(m2)
                        curmg = curmg + dtilrf(m,ls,k)*(phif(m,l,k))    &
                                      -dhatrf(m,ls,k)*phif(m,l,k)
                    enddo
                    cur2g = dtilr(m2,ls,k)*(phi(m2,l,k))
                    dhatr(m2,ls,k) = cur2g - curmg
                    dhatr(m2,ls,k) = dhatr(m2,ls,k) / (phi(m2,l,k))
                endif
            enddo      
        enddo
    enddo

    ! z-dir
    do l=1,nxy
        do k=1,nz
            ks=k
            do m2=1,ng2
                curmg = 0
                do m=mgb(m2),mge(m2)
                curmg = curmg - dtilzf(m,l,ks)*(phif(m,l,k)-phif(m,l,neibz(1,ks)))  &
                              - dhatzf(m,l,ks)*(phif(m,l,k)+phif(m,l,neibz(1,ks)))                
                enddo
                cur2g=-dtilz(m2,l,ks)*(phi(m2,l,k)-phi(m2,l,neibz(1,ks)))
                dhatz(m2,l,ks)= cur2g-curmg
                dhatz(m2,l,ks) = dhatz(m2,l,ks) / (phi(m2,l,k)+phi(m2,l,neibz(1,ks)))
            enddo 
        enddo

        k=nz
        ks=nzp1
        do m2=1,ng2
            curmg = 0
            do m=mgb(m2),mge(m2)
            curmg = curmg + dtilzf(m,l,ks)*(phif(m,l,k))    &
                          - dhatzf(m,l,ks)*phif(m,l,k)
            enddo
            cur2g = dtilz(m2,l,ks)*(phi(m2,l,k))
            dhatz(m2,l,ks) = cur2g - curmg
            dhatz(m2,l,ks) = dhatz(m2,l,ks) / (phi(m2,l,k))
        enddo  
    enddo
    
    return
end subroutine
