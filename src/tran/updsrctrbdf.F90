! 2012_09_26 . scb
!
! This subroutine determines parameters about transient source when bdf method is used
!
subroutine updsrctrbdf(phi,phif,psi,bdforder)
!   determine initial parameters
    use const
    use geom,       only : ng,nxy,nz,volnode,neibr,neibz
    use xsec,       only : xschif,xssfs,xssfe,xssf,mgb,mge
    use cmfdmg,     only : diagf,ccrf,cczf
    use cmfd2g,     only : am,ccr,ccz
    use tran
    use tranxsec,   only : lmbdk,xschifd,xschid
    use bdf
    use geomhex,    only : ifhexsp3
    use tpen_sp3,   only : hflx2   ! 2012_12_07 . scb
    use geom,       only : ifrecsp3  ! 2013_08_12 . scb
    implicit none
    
    real,   pointer         :: phi(:,:,:)
    real,   pointer         :: phif(:,:,:)
    real,   pointer         :: psi(:,:)

    integer :: bdforder       

    integer                 :: m,l,k,mp,ms,m2,idir,idirz,ii
    real                    :: rldt,rldtgp1,cappab,cappap1,         &
                               capbrldt,capbrldt2,srctrd,           &
                               sdtilsum,wn,srctrdtot,               &
                               gamma,rgamma,rgammap1,               &
                               wnp1,scatsrc,aphi,Mphif
    
    real :: omega(PREV:CURR, nprec), lamprec
    
!    real,pointer :: bdfcoef(:),phibdf(:,:,:,:),phifbdf(:,:,:,:)
!    
!    bdfcoef => bdfarray%bdfcoef
!    phibdf => bdfarray%phibdf
!    phifbdf => bdfarray%phifbdf

    gamma=deltm/deltmd
    rgamma=ONE/gamma
    rgammap1=ONE/(gamma+ONE)

!!!!!  scb added
    if(.not.automain) then
        psitd=psit
        psit=psi
    endif
!!!!!  added end
    
    do k=1,nz
    do l=1,nxy

!   calculate omega that is used for calculating source by precursor terms
!   betaomegat - sum of beta*omega that is used while delayed source term
!   sdtilsum - xschifd * sum of (lambda*cappa*prec + omega(prev)*psitd + omega(curr)*psit
!   wnp1 - sum of ( beta * omega_n+1)
        sdtil(:,l,k)=0
        wnp1=0
        
        do mp=1,nprec
!       local time dependent constants
            rldt=0
            if(lmbdk(mp,l,k).ne.0) rldt=ONE/(lmbdk(mp,l,k)*deltmd) 
            ! 1/(lambdk*deltm*(gamma+1))
            rldtgp1=rldt*rgammap1   

            !cappa bar=1-cappa
            !cappa is made in updcnsttr.f90
            cappab=ONE-cappa(mp,l,k)     
            cappap1=cappa(mp,l,k)+ONE
            capbrldt=cappab*rldt
            capbrldt2=cappab*(1.-2*rldt)

            omega(PREV,mp)=rldtgp1*(2*capbrldt-gamma*cappap1)
            omega(CURR,mp)=rldt*(cappap1+capbrldt2*rgamma)-cappa(mp,l,k)
            ! betaomegan=beta/lam*omega(n+1)
            betaomegan(mp,l,k)=betalam(mp,l,k)*(1-rldtgp1*(2+capbrldt2*rgamma))
            ! betalam : betak / lamdk 
            sdtil(mp,l,k)=betalam(mp,l,k)*(omega(PREV,mp)*psitd(l,k)+omega(CURR,mp)*psit(l,k))
            wnp1=wnp1+betaomegan(mp,l,k)*lmbdk(mp,l,k)
        enddo

!   srctrd = the previous transient source. xschif*(1-betat)*psif + xschifd * sum of lamda*prec + scattering src
        srctrd=0
        sdtilsum=0
        do mp=1,nprec
            lamprec=lmbdk(mp,l,k)*prec(mp,l,k)
            srctrd=srctrd+lamprec
            sdtil(mp,l,k)=sdtil(mp,l,k)+cappa(mp,l,k)*prec(mp,l,k)
            sdtilsum=sdtilsum+sdtil(mp,l,k)*lmbdk(mp,l,k)
        enddo

!   wn - previous wnp1. betap includes wn        
        wn=betap(l,k)+betat(l,k)-ONE
        
        if(ng.ne.ng2) then
            ! 2013_08_12 . scb
            !if(ifhexsp3) then
            !  stop "MG BDF transient is not implemented"
            !endif
            ! added end
            do m=1,ng
!                srctrdtot=xschifd(m,l,k)*srctrd+xschif(m,l,k)*(ONE-betat(l,k))*psi(l,k)
                srctrdtot=xschifd(m,l,k)*srctrd+xschifp(m,l,k)*(ONE-betat(l,k))*psi(l,k) 
                scatsrc=0.d0
                do ms=xssfs(m,l,k),xssfe(m,l,k)
                    scatsrc=scatsrc+phif(ms,l,k)*xssf(ms,m,l,k)
                enddo

                srctrdtot=srctrdtot+scatsrc*volnode(l,k)
                
                srctrf(m,l,k)=0.d0
                do ii=1,bdforder
                  srctrf(m,l,k)=srctrf(m,l,k)-bdfcoef(ii)*phifbdf(m,l,k,ii)
                enddo
                srctrf(m,l,k)=srctrf(m,l,k)*volnode(l,k)*rveloftm(m,l,k)/bdfarray%bdfcoef(0)+  &
                              xschifd(m,l,k)*sdtilsum      
!               chibetapf - used for calculating prompt fission source in CMFD.
!                chibetapf(m,l,k)=xschif(m,l,k)*(1-betat(l,k))+xschifd(m,l,k)*wnp1
                chibetapf(m,l,k)=xschifp(m,l,k)*(1-betat(l,k))+xschifd(m,l,k)*wnp1
            enddo
            
            do m2=1,ng2
                srctr(m2,l,k)=0.d0
                do m=mgb(m2),mge(m2)
                    srctr(m2,l,k)=srctr(m2,l,k)+srctrf(m,l,k)
                enddo
            enddo
            
            phifd(:,l,k)=phif(:,l,k)
            phid(:,l,k)=phi(:,l,k)
            
! 2013_08_12 . scb            
            if(ifhexsp3 .or. ifrecsp3) then
              do m=1,ng
                srctrf2(m,l,k)=0.d0
                do ii=1,bdforder
                  srctrf2(m,l,k)=srctrf2(m,l,k)-bdfcoef(ii)*phifbdf2(m,l,k,ii)
                enddo             
                !srctrf2(m,l,k)=5.d0/3.d0*rvelof(m,l,k)*srctrf2(m,l,k)
                srctrf2(m,l,k)=rvelof(m,l,k)*srctrf2(m,l,k)
              enddo
              
              do m2=1,ng2
                  srctr2(m2,l,k)=0.d0
                  do m=mgb(m2),mge(m2)
                      srctr2(m2,l,k)=srctr2(m2,l,k)+srctrf2(m,l,k)
                  enddo
              enddo
            endif  
! added end            
        else
            do m2=1,ng2
              srctr(m2,l,k)=0.d0
              do ii=1,bdforder
                srctr(m2,l,k)=srctr(m2,l,k)-bdfcoef(ii)*phibdf(m2,l,k,ii)
              enddo
              srctr(m2,l,k)=srctr(m2,l,k)*volnode(l,k)*rvelotm(m2,l,k)/bdfarray%bdfcoef(0)+  &
                            xschid(m2,l,k)*sdtilsum                        
            enddo
            ! srctr = s_dtil + 1/v*sum(bdfcoef*phi), rveloftm=a0/v in BDF routines
            
            if(ifhexsp3 .or. ifrecsp3) then  ! 2013_08_12 . scb
              do m2=1,ng2
                srctr2(m2,l,k)=0.d0
                do ii=1,bdforder
                  srctr2(m2,l,k)=srctr2(m2,l,k)-bdfcoef(ii)*phibdf2(m2,l,k,ii)
                enddo             
                !srctr2(m2,l,k)=srctr2(m2,l,k)-bdfcoef(0)*hflx2(2,m2,l,k)
                !srctr2(m2,l,k)=5.d0/3.d0*rvelo(m2,l,k)*srctr2(m2,l,k)   ! 2013_08_12 . scb
                srctr2(m2,l,k)=rvelo(m2,l,k)*srctr2(m2,l,k)   ! 2013_08_12 . scb
              enddo
            endif            
        endif

!       betap includes wnp1       
        betap(l,k)=one-betat(l,k)+wnp1
        
    enddo
    enddo   
    
    return 
end subroutine
! added end