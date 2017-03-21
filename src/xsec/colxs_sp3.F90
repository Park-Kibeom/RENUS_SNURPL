! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
!subroutine colxs_sp3(fphi2)
subroutine colxs_sp3(fphi3)
    use const
    use xsec
    use geom,   only        :   ng,nxy,nz
    use sfam,   only        :   fphiscb => fphi
    implicit none
    
    !real,pointer            ::  fphi2(:,:,:)
    real,pointer            ::  fphi3(:,:,:)
    
    integer                 :: l,k,m,m2
    
    if(ng.eq.ng2) return
   
! collapse into two group xsec
    do k=1,nz
        do l=1,nxy
            do m2=1,ng2
                xsd2(m2,l,k)=0
                xstr(m2,l,k)=0
                do m=mgb(m2),mge(m2)
                    xsd2(m2,l,k)=xsd2(m2,l,k)+fphiscb(m,l,k)*xsdf2(m,l,k)
                    xstr(m2,l,k)=xstr(m2,l,k)+fphiscb(m,l,k)*xstrf(m,l,k)
                enddo
            enddo ! of m2
        enddo !l
    enddo !k

    return
end subroutine
