! added in ARTOS ver. 0.2 . 2012_07_19 by SCB.
subroutine updsrctr_hex(phi,phif,psi)
!   determine initial parameters
    use const
    use geom,       only : ng,nxy,nz,volnode,neibr,neibz
    use xsec,       only : xschif,xssfs,xssfe,xssf,mgb,mge
    use cmfdmg,     only : diagf,ccrf,cczf
    use cmfd2g,     only : am,ccr,ccz
    use tran
    use tranxsec,   only : lmbdk,xschifd,xschid
	  use cmfdhex2g,  only : dcmat, cmat
	  use geomhex,    only : ipntr

    implicit none
    
    real,   pointer         :: phi(:,:,:)
    real,   pointer         :: phif(:,:,:)
    real,   pointer         :: psi(:,:)

    integer                 :: m,l,k,mp,ms,m2,idir,idirz,ln
    real                    :: rldt,rldtgp1,cappab,cappap1,         &
                               capbrldt,capbrldt2,srctrd,           &
                               sdtilsum,wn,srctrdtot,               &
                               gamma,rgamma,rgammap1,               &
                               wnp1,scatsrc,aphi,Mphif
    
    real :: omega(PREV:CURR, nprec), lamprec


    gamma=deltm/deltmd
    rgamma=ONE/gamma
    rgammap1=ONE/(gamma+ONE)

    do k=1,nz
    do l=1,nxy

        psitd(l,k)=psit(l,k)
        psit(l,k)=psi(l,k)        

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
            rldtgp1=rldt*rgammap1   

            cappab=ONE-cappa(mp,l,k)     
            cappap1=cappa(mp,l,k)+ONE
            capbrldt=cappab*rldt
            capbrldt2=cappab*(1.-2*rldt)

            omega(PREV,mp)=rldtgp1*(2*capbrldt-gamma*cappap1)
            omega(CURR,mp)=rldt*(cappap1+capbrldt2*rgamma)-cappa(mp,l,k)
            betaomegan(mp,l,k)=betalam(mp,l,k)*(1-rldtgp1*(2+capbrldt2*rgamma))
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
            do m=1,ng
!                srctrdtot=xschifd(m,l,k)*srctrd+xschif(m,l,k)*(ONE-betat(l,k))*psi(l,k)    ! added in ARTOS ver. 0.2 . 2012_08_06 by SCB    
				        srctrdtot=xschifd(m,l,k)*srctrd+xschifp(m,l,k)*(ONE-betat(l,k))*psi(l,k)    ! added in ARTOS ver. 0.2 . 2012_08_06 by SCB     ! mhtr

                scatsrc=0
                do ms=xssfs(m,l,k),xssfe(m,l,k)
                    scatsrc=scatsrc+phif(ms,l,k)*xssf(ms,m,l,k)
                enddo

                srctrdtot=srctrdtot+scatsrc*volnode(l,k)

                Mphif=diagf(l,k,m)*phif(m,l,k)                       
                do idir=1,nrdir2
                    Mphif=Mphif+ccrf(idir,l,k,m)*phif(m,neibr(idir,l),k)
                enddo
                do idirz=1,2
                    Mphif=Mphif+cczf(idirz,l,k,m)*phif(m,l,neibz(idirz,k))
                enddo     
                srctrf(m,l,k)=  xschifd(m,l,k)*sdtilsum+                        &
                                expof(m,l,k)*(                                         &
                                    volnode(l,k)*rveloftm(m,l,k)*phif(m,l,k)   &
                                   +thetabar*(srctrdtot-Mphif)                  &
                                   -thetabar*pinvf(m,l,k)*volnode(l,k)*phif(m,l,k) &
                                )

!               chibetapf - used for calculating prompt fission source in CMFD.
!                chibetapf(m,l,k)=xschif(m,l,k)*(1-betat(l,k))+xschifd(m,l,k)*wnp1    ! added in ARTOS ver. 0.2 . 2012_08_06 by SCB    
  				      chibetapf(m,l,k)=xschifp(m,l,k)*(1-betat(l,k))+xschifd(m,l,k)*wnp1    ! added in ARTOS ver. 0.2 . 2012_08_06 by SCB     ! mhtr
            enddo
            
            do m2=1,ng2
                srctr(m2,l,k)=0
                do m=mgb(m2),mge(m2)
                    srctr(m2,l,k)=srctr(m2,l,k)+srctrf(m,l,k)
                enddo
            enddo
            
            phid(:,l,k)=phi(:,l,k)
            phifd(:,l,k)=phif(:,l,k)
                        
            if(exptheta) then
!                phif(:,l,k)=phif(:,l,k)*expof(:,l,k)    ! 2012_11_09 . scb
				      do m2=1,ng2
					      phi(m2,l,k)=0
					      do m=mgb(m2),mge(m2)
						      phi(m2,l,k)=phi(m2,l,k)+phif(m,l,k)
					      enddo
				      enddo
            endif
        else
            do m2=1,ng2
                aphi=dcmat(indm24(m2,FAST),l,k)*phi(FAST,l,k)+                     &
                     dcmat(indm24(m2,THERMAL),l,k)*phi(THERMAL,l,k)
                     
				        do idir=1,ipntr(0,l)
					        aphi=aphi-cmat(m2,idir,l,k)*phi(m2,ipntr(idir,l),k)
				        enddo
				        do idirz=1,2
					        aphi=aphi-cmat(m2,idirz+6,l,k)*phi(m2,l,neibz(idirz,k))
				        enddo

                aphi=aphi+xschid(m2,l,k)*wn*psi(l,k)
                srctr(m2,l,k)=  xschid(m2,l,k)*sdtilsum+                        &
                                expo(m2,l,k)*(                                         &
                                    volnode(l,k)*rvelotm(m2,l,k)*phi(m2,l,k)+   &
                                    thetabar*(xschid(m2,l,k)*srctrd-aphi)-       &
                                    thetabar*pinv(m2,l,k)*volnode(l,k)*phi(m2,l,k)  &
                                )
            enddo
            if(exptheta) then
                phid(:,l,k)=phi(:,l,k)
!                phi(:,l,k)=phi(:,l,k)*expo(:,l,k)       ! 2012_11_09 . scb
            endif        
        endif

!       betap includes wnp1
        betap(l,k)=one-betat(l,k)+wnp1
        
    enddo
    enddo   
    
    return 
end subroutine
