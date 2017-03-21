subroutine resetsenm2n
    use const
    use senm2n
    use geom,   only : ng,nxy,nz,hmesh
    use xsec,   only : xstf,xsdf
    implicit none
    
    integer                 :: idir,l,k,m
    real                    :: rkp4c
    
    do idir=1, ndirmax
    do k=1,nz
        do l=1,nxy
            do m=1,ng
                kp2(m,l,k,idir)=xstf(m,l,k)*hmesh(idir,l,k)**2/(4*xsdf(m,l,k))
                kp(m,l,k,idir)=sqrt(kp2(m,l,k,idir))
                if(not(kp2(m,l,k,idir) == 0)) rkp2(m,l,k,idir)=1/kp2(m,l,k,idir) !!!!!!!!! error for hexagonal geometry
                sinhkp(m,l,k,idir)=sinh(kp(m,l,k,idir))
                coshkp(m,l,k,idir)=cosh(kp(m,l,k,idir))
                if(not(sinhkp(m,l,k,idir) == 0)) cschkp(m,l,k,idir)=1/sinhkp(m,l,k,idir) !!!!!!!!! error for hexagonal geometry
                rkp4c = rkp2(m,l,k,idir) * rkp2(m,l,k,idir)
! obtain constants, used for updating projected flux polynomial coeff
                phicnst(0,m,l,k,idir)=sinhkp(m,l,k,idir)*sqrt(rkp2(m,l,k,idir))
                phicnst(1,m,l,k,idir)=3*rkp2(m,l,k,idir)*(                                          &
                                        coshkp(m,l,k,idir)*kp(m,l,k,idir)-                          &
                                        sinhkp(m,l,k,idir)                                          &
                                      )
                phicnst(2,m,l,k,idir)=5*rkp4c*(                                                     &
                                        -3*coshkp(m,l,k,idir)*kp2(m,l,k,idir)+                      &
                                        sinhkp(m,l,k,idir)*kp(m,l,k,idir)*(                         &
                                            3+kp2(m,l,k,idir)                                       &
                                        )                                                           &
                                      )
                phicnst(3,m,l,k,idir)=7*rkp4c*(                                                     &
                                        coshkp(m,l,k,idir)*kp(m,l,k,idir)*(                         &
                                            15+kp2(m,l,k,idir)                                      &
                                        )-3*sinhkp(m,l,k,idir)*(                                    &
                                            5+2*kp2(m,l,k,idir)                                     &
                                        )                                                           &
                                      )
                phicnst(4,m,l,k,idir)=9*rkp4c*(                                                     &
                                        -5*coshkp(m,l,k,idir)*(                                     &
                                            21+2*kp2(m,l,k,idir)                                    &
                                        )+sinhkp(m,l,k,idir)*(                                      &
                                            105*sqrt(rkp2(m,l,k,idir))+45*kp(m,l,k,idir)+           &
                                            kp2(m,l,k,idir)*kp(m,l,k,idir)                          &
                                        )                                                           &
                                      )
            enddo ! m
        enddo 
    enddo 
    enddo ! idir
    
    return
	end subroutine