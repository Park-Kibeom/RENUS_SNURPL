subroutine calAby1n(idir,l,k,m,bnd,phifs,psol,hsolb,hsola)
    use const
    use senm2n
    use geom,   only : albedo
    use xsec,   only : xbeta
    use xsec,   only : xsadf   ! 2015_08_04 . scb
    implicit none
    
    integer, intent(in)     :: bnd              ! LEFT, RIGHT
    real,    intent(in)     :: phifs
    real,    intent(in)     :: psol(0:4),hsolb  ! PARTICULAR, B
    real,    intent(out)    :: hsola            ! A

    integer                 :: l,k,m,idir
    real                    :: cschkps,sinhkps,coshkps,kps,betas
    real                    :: adf   ! 2015_08_04 . scb

    cschkps=cschkp(m,l,k,idir)
    sinhkps=sinhkp(m,l,k,idir)
    coshkps=coshkp(m,l,k,idir)
    kps=kp(m,l,k,idir)
    betas=xbeta(m,l,k,idir)
    
! 2015_08_04 . scb    
    adf=1
    !if(idir.le.YDIR) then
    !    adf=xsadf(idir*2-2+bnd,m,l,k)   ! Left B : ADF(1 or 3), Right B : ADF(2 or 4)
    !endif      
! added end
! In 150804, ADF was also applied to the system boundary. 
! However, it can be unphiysical if reflector assembly is loaded at the core pheriphery.
! Thus remove ADF again. 

    !!!!! GO TO !!!!!!
    select case(bnd)
!   LEFT BOUNDARY
    case(LEFT)
        hsola=2*betas*(kps*sinhkps*hsolb-psol(1)+3*psol(2)-6*psol(3)+10*psol(4))+       &
              !(coshkps*hsolb+psol(0)-psol(1)+psol(2)-psol(3)+psol(4))*albedo(LEFT,idir)
              (coshkps*hsolb+psol(0)-psol(1)+psol(2)-psol(3)+psol(4))*albedo(LEFT,idir)*adf
        !hsola=hsola/(albedo(LEFT,idir)*sinhkps+2*betas*kps*coshkps)
        hsola=hsola/(albedo(LEFT,idir)*adf*sinhkps+2*betas*kps*coshkps)
    goto 100            
!   RIGHT BOUNDARY
    case(RIGHT)
        hsola=-2*betas*(kps*sinhkps*hsolb+psol(1)+3*psol(2)+6*psol(3)+10*psol(4))-      &
              !(coshkps*hsolb+psol(0)+psol(1)+psol(2)+psol(3)+psol(4))*albedo(RIGHT,idir)
              (coshkps*hsolb+psol(0)+psol(1)+psol(2)+psol(3)+psol(4))*albedo(RIGHT,idir)*adf
        !hsola=hsola/(albedo(RIGHT,idir)*sinhkps + 2*betas*kps*coshkps)
        hsola=hsola/(albedo(RIGHT,idir)*adf*sinhkps + 2*betas*kps*coshkps)
    goto 100            
    end select
    
100 return
end subroutine
