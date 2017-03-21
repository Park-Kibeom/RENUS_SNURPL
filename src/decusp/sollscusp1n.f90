subroutine sollscusp1n(ccz, diag, src, phi1n)
    use const
    use matop
    use mat2x2
    use decusping1n
    use geom, only : ng
    implicit none
    
    real              :: ccz(ng,BOTTOM:TOP,ntfine)
    real              :: diag(ng,ng,ntfine), phi1n(ng,ntfine), src(ng,ntfine)
    
    integer           :: k, j, m, ms
    real              :: ld(ng,ng), ldu(ng,ng), ldb(ng)
    
    if(ng.eq.ng2) then
        do k=1,ntfine-1
            call invmat2x2(diag(:,:,k),diag(:,:,k))   ! inverse of matrix
            call diagxmat(ccz(:,BOTTOM,k+1),diag(:,:,k),ld(:,:),ng)
            call matxdiag(ld(:,:),ccz(:,TOP,k),ldu,ng)
            diag(:,:,k+1) = diag(:,:,k+1) - ldu(:,:)
            call matxvec2x2(ld(:,:),src(:,k),ldb(:))
            src(:,k+1) = src(:,k+1)-ldb(:)
        enddo
        call invmat2x2(diag(:,:,ntfine),diag(:,:,ntfine))   ! inverse of matrix
        
        call matxvec2x2(diag(:,:,ntfine),src(:,ntfine),phi1n(:,ntfine))
        do k=ntfine-1, 1, -1
            ldb(:)=src(:,k)-ccz(:,TOP,k)*phi1n(:,k+1)
            call matxvec2x2(diag(:,:,k),ldb(:),phi1n(:,k))
        enddo        
    else
        do k=1,ntfine-1
            call invmat(diag(:,:,k),diag(:,:,k),ng)   ! inverse of matrix
            call diagxmat(ccz(:,BOTTOM,k+1),diag(:,:,k),ld(:,:),ng)
            call matxdiag(ld(:,:),ccz(:,TOP,k),ldu,ng)
            diag(:,:,k+1) = diag(:,:,k+1) - ldu(:,:)
            call matxvec(ld(:,:),src(:,k),ldb(:),ng)
            src(:,k+1) = src(:,k+1)-ldb(:)
        enddo
        call invmat(diag(:,:,ntfine),diag(:,:,ntfine),ng)   ! inverse of matrix
        
        call matxvec(diag(:,:,ntfine),src(:,ntfine),phi1n(:,ntfine),ng)
        do k=ntfine-1, 1, -1
            ldb(:)=src(:,k)-ccz(:,TOP,k)*phi1n(:,k+1)
            call matxvec(diag(:,:,k),ldb(:),phi1n(:,k),ng)
        enddo
    endif
    
!define DEBUG ! write linear system into a file.
#ifdef DEBUG
            k=1
            do m=1,ng
                do ms=1,ng
                    write(1000, '(1p,e15.6,$)') diag(ms,m,k)
                enddo
                do ms=1,m-1
                    write(1000, '(1p,e15.6,$)') 0.
                enddo
                write(1000, '(1p,e15.6,$)') ccz(m,TOP,k)
                do ms=m+1,ng
                    write(1000, '(1p,e15.6,$)') 0.
                enddo
                do j=(k+1)*ng+1,ntfine*ng
                    write(1000, '(1p,e15.6,$)') 0.
                enddo                
                write(1000, '()')
            enddo
            
            do k=2,ntfine-1
                do m=1,ng
                    do j=1,(k-2)*ng
                        write(1000, '(1p,e15.6,$)') 0.
                    enddo
                    do ms=1,m-1
                        write(1000, '(1p,e15.6,$)') 0.
                    enddo
                    write(1000, '(1p,e15.6,$)') ccz(m,BOTTOM,k)
                    do ms=m+1,ng
                        write(1000, '(1p,e15.6,$)') 0.
                    enddo
                    do ms=1,ng
                        write(1000, '(1p,e15.6,$)') diag(ms,m,k)
                    enddo
                    do ms=1,m-1
                        write(1000, '(1p,e15.6,$)') 0.
                    enddo
                    write(1000, '(1p,e15.6,$)') ccz(m,TOP,k)
                    do ms=m+1,ng
                        write(1000, '(1p,e15.6,$)') 0.
                    enddo
                    do j=(k+1)*ng+1,ntfine*ng
                        write(1000, '(1p,e15.6,$)') 0.
                    enddo
                    
                    write(1000, '()')
                enddo
            enddo
            k=ntfine
            do m=1,ng
                do j=1,(k-2)*ng
                    write(1000, '(1p,e15.6,$)') 0.
                enddo
                do ms=1,m-1
                    write(1000, '(1p,e15.6,$)') 0.
                enddo
                write(1000, '(1p,e15.6,$)') ccz(m,BOTTOM,k)
                do ms=m+1,ng
                    write(1000, '(1p,e15.6,$)') 0.
                enddo
                do ms=1,ng
                    write(1000, '(1p,e15.6,$)') diag(ms,m,k)
                enddo
                write(1000, '()')
            enddo
            write(1000, '(1p,e15.6)') src(:,:)
#endif          
end subroutine