subroutine setlscusp1n(hfine, xsd, xst, xss,xsnf,xschi, trlfine,  &
                   diag, ccz)
    use const
    use geom,   only : ng
    use decusping1n
    implicit none
    
    real                :: hfine(ntfine)
    real                :: xsd(ng,ntfine), xst(ng,ntfine),    &
                           xsnf(ng,ntfine), xss(ng,ng,ntfine),&
                           xschi(ng), trlfine(ng,ntfine)
    real                :: ccz(ng,BOTTOM:TOP,ntfine), diag(ng,ng,ntfine)

    integer             :: k, ks, m, mm, ms
    real                :: dtil1n(ng,ntfine+1)
    real                :: betal(ng), betar(ng)
    real                :: offdiag

! d-tilde
    k=1
    ks=1
    do m=1,ng
        betal(m)=xsd(m,k)/hfine(k)
        dtil1n(m,ks)=betal(m)*0.5
    enddo            

    do k=2,ntfine
        ks=k
        do m=1,ng
            betar(m)=xsd(m,k)/hfine(k)
            dtil1n(m,ks)=2*betal(m)*betar(m)/(betal(m)+betar(m))
            betal(m)=betar(m)
        enddo 
    enddo

    k=ntfine
    ks=ntfine+1
    do m=1,ng
        dtil1n(m,ks)=betar(m)*0.5
    enddo 

!   contructing M-matrix
    do k=1, ntfine
        do m=1,ng
!             fission matrix        
              do mm=1,ng
                  diag(mm,m,k)=-xschi(m)*xsnf(mm,k)*hfine(k)
              enddo
!             scattering matrix
              do ms=1,ng
                  diag(ms,m,k)=diag(ms,m,k)-xss(ms,m,k)*hfine(k)
              enddo
!             removal matrix              
              diag(m,m,k)=diag(m,m,k)+xst(m,k)*hfine(k)

!             coupling term to currents
              offdiag = dtil1n(m,k)
              ccz(m,BOTTOM,k)=-offdiag
              diag(m,m,k)=diag(m,m,k)+offdiag

              offdiag = dtil1n(m,k+1)
              ccz(m,TOP,k)=-offdiag
              diag(m,m,k)=diag(m,m,k)+offdiag
        enddo
    enddo

    
    return
end subroutine