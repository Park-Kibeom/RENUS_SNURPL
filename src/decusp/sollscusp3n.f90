subroutine sollscusp3n
!   solve 3 node problem with blockwise gauss elimination
    use const
    use decusping3n
    use mat2x2,         only  : invmat2x2
    implicit none
    
    
    integer                   :: k,kp1,kl,klp1,onezero,nmesh2
    real                      :: mult(2,2)
    real                      :: coll(2,2,2:ntfine),rowl(4,2:ntfine+1),x1(ng2),rhs(ng2)
    
    nmesh2=nmesh*2
    
    diag1d(:,1)=(/1,0,0,1/)
    diag1d(:,ntfine)=diag1d(:,1)
    
    coll(:,:,2)=0; coll(1,1,2)=bot1d(1,2); coll(2,2,2)=bot1d(2,2)
    coll(:,:,ntfine)=0
    rowl(:,2)=(/1,0,0,1/)
    
    
    do k=2, ntfine-2
        kp1=k+1
        call invmat2x2(diag1d(:,k), diag1d(:,k))
        
!      eliminate lower diagonal
        mult(1,1)=bot1d(1,kp1)*diag1d(1,k)
        mult(2,1)=bot1d(1,kp1)*diag1d(2,k)
        mult(1,2)=bot1d(2,kp1)*diag1d(3,k)
        mult(2,2)=bot1d(2,kp1)*diag1d(4,k)        
         
        diag1d(1,kp1)=diag1d(1,kp1)-mult(1,1)*top1d(1,k)
        diag1d(2,kp1)=diag1d(2,kp1)-mult(2,1)*top1d(2,k)
        diag1d(3,kp1)=diag1d(3,kp1)-mult(1,2)*top1d(1,k)
        diag1d(4,kp1)=diag1d(4,kp1)-mult(2,2)*top1d(2,k)
         
        coll(1,1,kp1)=-mult(1,1)*coll(1,1,k)-mult(2,1)*coll(1,2,k)
        coll(2,1,kp1)=-mult(1,1)*coll(2,1,k)-mult(2,1)*coll(2,2,k)
        coll(1,2,kp1)=-mult(1,2)*coll(1,1,k)-mult(2,2)*coll(1,2,k)
        coll(2,2,kp1)=-mult(1,2)*coll(2,1,k)-mult(2,2)*coll(2,2,k)
         
        src1d(1,kp1)=src1d(1,kp1)-mult(1,1)*src1d(1,k)-mult(2,1)*src1d(2,k)
        src1d(2,kp1)=src1d(2,kp1)-mult(1,2)*src1d(1,k)-mult(2,2)*src1d(2,k)

!      eliminate element on the last row
         mult(1,1)=rowl(1,k)*diag1d(1,k)+rowl(2,k)*diag1d(3,k)
         mult(2,1)=rowl(1,k)*diag1d(2,k)+rowl(2,k)*diag1d(4,k)
         mult(1,2)=rowl(3,k)*diag1d(1,k)+rowl(4,k)*diag1d(3,k)
         mult(2,2)=rowl(3,k)*diag1d(2,k)+rowl(4,k)*diag1d(4,k)

         onezero=0
         if(k.lt.nmesh) onezero=1
         
         rowl(1,kp1)=onezero-mult(1,1)*top1d(1,k)
         rowl(2,kp1)= -mult(2,1)*top1d(2,k)
         rowl(3,kp1)= -mult(1,2)*top1d(1,k)
         rowl(4,kp1)=onezero-mult(2,2)*top1d(2,k)
         
         diag1d(1,1)=diag1d(1,1)-mult(1,1)*coll(1,1,k)-mult(2,1)*coll(1,2,k)
         diag1d(2,1)=diag1d(2,1)-mult(1,1)*coll(2,1,k)-mult(2,1)*coll(2,2,k)
         diag1d(3,1)=diag1d(3,1)-mult(1,2)*coll(1,1,k)-mult(2,2)*coll(1,2,k)
         diag1d(4,1)=diag1d(4,1)-mult(1,2)*coll(2,1,k)-mult(2,2)*coll(2,2,k)
         src1d(1,1)=src1d(1,1)-mult(1,1)*src1d(1,k)-mult(2,1)*src1d(2,k)
         src1d(2,1)=src1d(2,1)-mult(1,2)*src1d(1,k)-mult(2,2)*src1d(2,k)

!      eliminate element on the last-1 row
         if(k.gt.nmesh2 .and. k.lt.ntfine-1) then
            kl=k-nmesh2+1
            klp1=kl+1
            
            mult(1,1)=rowl(1,kl)*diag1d(1,k)+rowl(2,kl)*diag1d(3,k)
            mult(2,1)=rowl(1,kl)*diag1d(2,k)+rowl(2,kl)*diag1d(4,k)
            mult(1,2)=rowl(3,kl)*diag1d(1,k)+rowl(4,kl)*diag1d(3,k)
            mult(2,2)=rowl(3,kl)*diag1d(2,k)+rowl(4,kl)*diag1d(4,k)
            
            rowl(1,klp1)=1-mult(1,1)*top1d(1,k)
            rowl(2,klp1)= -mult(2,1)*top1d(2,k)
            rowl(3,klp1)= -mult(1,2)*top1d(1,k)
            rowl(4,klp1)=1-mult(2,2)*top1d(2,k)
            
            coll(1,1,ntfine)=coll(1,1,ntfine)-mult(1,1)*coll(1,1,k)-mult(2,1)*coll(1,2,k)
            coll(2,1,ntfine)=coll(2,1,ntfine)-mult(1,1)*coll(2,1,k)-mult(2,1)*coll(2,2,k)
            coll(1,2,ntfine)=coll(1,2,ntfine)-mult(1,2)*coll(1,1,k)-mult(2,2)*coll(1,2,k)
            coll(2,2,ntfine)=coll(2,2,ntfine)-mult(1,2)*coll(2,1,k)-mult(2,2)*coll(2,2,k)
            src1d(1,ntfine)=src1d(1,ntfine)-mult(1,1)*src1d(1,k)-mult(2,1)*src1d(2,k)
            src1d(2,ntfine)=src1d(2,ntfine)-mult(1,2)*src1d(1,k)-mult(2,2)*src1d(2,k)
         endif
    enddo
    
!   last-1 row
    k=ntfine-1
    kp1=k+1
    kl=k-nmesh2+1
!
!      store inverse
    call invmat2x2(diag1d(:,k),diag1d(:,k))
!
!      eliminate lower diagonal
    mult(1,1)=rowl(1,kl)*diag1d(1,k)+rowl(2,kl)*diag1d(3,k)
    mult(2,1)=rowl(1,kl)*diag1d(2,k)+rowl(2,kl)*diag1d(4,k)
    mult(1,2)=rowl(3,kl)*diag1d(1,k)+rowl(4,kl)*diag1d(3,k)
    mult(2,2)=rowl(3,kl)*diag1d(2,k)+rowl(4,kl)*diag1d(4,k)

    diag1d(1,kp1)=diag1d(1,kp1)-mult(1,1)*top1d(1,k)
    diag1d(2,kp1)=diag1d(2,kp1)-mult(2,1)*top1d(2,k)
    diag1d(3,kp1)=diag1d(3,kp1)-mult(1,2)*top1d(1,k)
    diag1d(4,kp1)=diag1d(4,kp1)-mult(2,2)*top1d(2,k)

    coll(1,1,kp1)=coll(1,1,kp1)-mult(1,1)*coll(1,1,k)-mult(2,1)*coll(1,2,k)
    coll(2,1,kp1)=coll(2,1,kp1)-mult(1,1)*coll(2,1,k)-mult(2,1)*coll(2,2,k)
    coll(1,2,kp1)=coll(1,2,kp1)-mult(1,2)*coll(1,1,k)-mult(2,2)*coll(1,2,k)
    coll(2,2,kp1)=coll(2,2,kp1)-mult(1,2)*coll(2,1,k)-mult(2,2)*coll(2,2,k)
    
    src1d(1,kp1)=src1d(1,kp1)-mult(1,1)*src1d(1,k)-mult(2,1)*src1d(2,k)
    src1d(2,kp1)=src1d(2,kp1)-mult(1,2)*src1d(1,k)-mult(2,2)*src1d(2,k)
    
    mult(1,1)=rowl(1,k)*diag1d(1,k)+rowl(2,k)*diag1d(3,k)
    mult(2,1)=rowl(1,k)*diag1d(2,k)+rowl(2,k)*diag1d(4,k)
    mult(1,2)=rowl(3,k)*diag1d(1,k)+rowl(4,k)*diag1d(3,k)
    mult(2,2)=rowl(3,k)*diag1d(2,k)+rowl(4,k)*diag1d(4,k)

    rowl(1,kp1)= -mult(1,1)*top1d(1,k)
    rowl(2,kp1)= -mult(2,1)*top1d(2,k)
    rowl(3,kp1)= -mult(1,2)*top1d(1,k)
    rowl(4,kp1)= -mult(2,2)*top1d(2,k)

    diag1d(1,1)=diag1d(1,1)-mult(1,1)*coll(1,1,k)-mult(2,1)*coll(1,2,k)
    diag1d(2,1)=diag1d(2,1)-mult(1,1)*coll(2,1,k)-mult(2,1)*coll(2,2,k)
    diag1d(3,1)=diag1d(3,1)-mult(1,2)*coll(1,1,k)-mult(2,2)*coll(1,2,k)
    diag1d(4,1)=diag1d(4,1)-mult(1,2)*coll(2,1,k)-mult(2,2)*coll(2,2,k)

    src1d(1,1)=src1d(1,1)-mult(1,1)*src1d(1,k)-mult(2,1)*src1d(2,k)
    src1d(2,1)=src1d(2,1)-mult(1,2)*src1d(1,k)-mult(2,2)*src1d(2,k)
!
! last row
    k=ntfine
    
!   store inverse
    call invmat2x2(diag1d(:,k), diag1d(:,k))
!
!      eliminate lower diagonal
    mult(1,1)=rowl(1,k)*diag1d(1,k)+rowl(2,k)*diag1d(3,k)
    mult(2,1)=rowl(1,k)*diag1d(2,k)+rowl(2,k)*diag1d(4,k)
    mult(1,2)=rowl(3,k)*diag1d(1,k)+rowl(4,k)*diag1d(3,k)
    mult(2,2)=rowl(3,k)*diag1d(2,k)+rowl(4,k)*diag1d(4,k)
    
    diag1d(1,1)=diag1d(1,1)-mult(1,1)*coll(1,1,k)-mult(2,1)*coll(1,2,k)
    diag1d(2,1)=diag1d(2,1)-mult(1,1)*coll(2,1,k)-mult(2,1)*coll(2,2,k)
    diag1d(3,1)=diag1d(3,1)-mult(1,2)*coll(1,1,k)-mult(2,2)*coll(1,2,k)
    diag1d(4,1)=diag1d(4,1)-mult(1,2)*coll(2,1,k)-mult(2,2)*coll(2,2,k)
    
    src1d(1,1)=src1d(1,1)-mult(1,1)*src1d(1,k)-mult(2,1)*src1d(2,k)
    src1d(2,1)=src1d(2,1)-mult(1,2)*src1d(1,k)-mult(2,2)*src1d(2,k)

!
! backward substitution
    k=1
!   store inverse
    call invmat2x2(diag1d(:,k),diag1d(:,k))

!   solution at the first node
    x1(1)=diag1d(1,k)*src1d(1,k)+diag1d(2,k)*src1d(2,k)
    x1(2)=diag1d(3,k)*src1d(1,k)+diag1d(4,k)*src1d(2,k)
    phi1d(1,k)=x1(1)
    phi1d(2,k)=x1(2)
!
!   solution at the last node
    k=ntfine
    rhs(1)=src1d(1,k)-coll(1,1,k)*x1(1)-coll(2,1,k)*x1(2)
    rhs(2)=src1d(2,k)-coll(1,2,k)*x1(1)-coll(2,2,k)*x1(2)
    phi1d(1,k)=diag1d(1,k)*rhs(1)+diag1d(2,k)*rhs(2)
    phi1d(2,k)=diag1d(3,k)*rhs(1)+diag1d(4,k)*rhs(2)
!
    kp1=k
    do k=ntfine-1,2,-1
       rhs(1)=src1d(1,k)-coll(1,1,k)*x1(1)-coll(2,1,k)*x1(2)                &
   &         -top1d(1,k)*phi1d(1,kp1)
       rhs(2)=src1d(2,k)-coll(1,2,k)*x1(1)-coll(2,2,k)*x1(2)                &
   &         -top1d(2,k)*phi1d(2,kp1)
       phi1d(1,k)=diag1d(1,k)*rhs(1)+diag1d(2,k)*rhs(2)
       phi1d(2,k)=diag1d(3,k)*rhs(1)+diag1d(4,k)*rhs(2)
       kp1=k
    enddo

end subroutine