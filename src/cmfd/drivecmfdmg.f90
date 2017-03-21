subroutine drivecmfdmg(iftran,chkconv,ncmfd,ibeg,nintot,eigv,reigv,chi,phif,epsl2,erreig,errl2)

    use const
    use cmfdmg
    use sfam,       only : psi
    use sfam_cntl,  only : nupmax,maxupscatgr,ninmax,epsbicgmg,epscmfdmg
    
    use geom,       only : ng,nxy,nz,kfbeg,kfend,volnode,nxsf,nxef,nodel, &
                           neibr,neibz  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB. ! mhtr
    use xsec,       only : xschif,xssf,xssfs,xssfe,xsnff
    use tran,       only : rveloftm,srctrf
    use bicgmg,     only : initbicg1g,solbicg1g 
    use sfam,       only : psid    ! 2013_05_14 . scb
    use geom,       only : ifrecsp3    ! 2013_07_29 . scb
    use geomhex,    only : ifhexsp3    ! 2014_02_24 . scb
    implicit none

    logical                 :: iftran,chkconv
    integer                 :: ncmfd
    integer                 :: ibeg,nintot
    real                    :: eigv,reigv
    real,pointer            :: chi(:,:,:)
    real,pointer            :: phif(:,:,:)
    real                    :: epsl2,erreig,errl2
    
    integer                 :: iout,iup,m,l,k,mbeg,ms,iin,i,j
    real                    :: err2d,err2,err,ss,r20,r2,  &
                               !psipsid,psipsi,psid,domr,eigvd,  &
                               psipsid,psipsi,domr,eigvd,  &   ! 2013_05_14 . scb
                               resid,resid0,relresid
    character               :: mesg*120
    logical                 :: converged

! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.        
    integer                 :: negative, idir, iphi
    real                    :: sumphi
! added end

! 2013_07_30 . scb
    real                    :: eigvfirst
    common / convcheck / eigvfirst
! added end    
	
    err2=big

!    write(mesg,'(a5,a10,2a12,3a10,a14)') 'iout','eigv','err_eig','l2_norm','domr','resid','relresid','flux sum'
    write(mesg,'(a5,a10,2a12,3a10,a10,a10)') 'MG_i','eigv','err_eig','l2_norm','domr','resid','relresid','flux_sum','nega_MG'  ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    call message(true,true,mesg)

! 2013_07_29 . scb    
    if(ifrecsp3 .or. ifhexsp3) then
      ncmfd=5
      ninmax=5
    endif
! added end

    do iout=1,ncmfd
    !do iout=1,10    ! 2013_07_22 . scb
        do iup=1,nupmax
            mbeg=maxupscatgr
            if(iup.eq.1) mbeg=1

            do m=mbeg,ng
                do k=1,nz
                do l=1,nxy
                    srcf(m,l,k)=reigv*chi(m,l,k)*psi(l,k)
                    if(iftran) srcf(m,l,k)=srcf(m,l,k)+srctrf(m,l,k)
                    ss=0
                    do ms=xssfs(m,l,k),xssfe(m,l,k)
                        ss=ss+phif(ms,l,k)*xssf(ms,m,l,k)
                    enddo
                    srcf(m,l,k)=srcf(m,l,k)+ss*volnode(l,k)
                enddo
                enddo

                call initbicg1g(m,phif,srcf,r20)
                
                do iin=1, ninmax
                !do iin=1, 10    ! 2013_07_22 . scb
                    call solbicg1g(m,r20,phif,r2)
                    if(r2.le.epsbicgmg) exit
                enddo
                nintot=nintot+min(iin,ninmax)
            enddo ! of m              
        enddo ! of iup
        
        err2d=err2;err2=0;
        psipsid=0;psipsi=0
        do k=kfbeg,kfend
            do j=1,ny
            do i=nxsf(j),nxef(j)
                l=nodel(i,j)
                
                psid(l,k)=psi(l,k)
                psi(l,k)=0
                do m=1,ng
                    psi(l,k)=psi(l,k)+xsnff(m,l,k)*phif(m,l,k)
                enddo
                psi(l,k)=psi(l,k)*volnode(l,k)
                
                psipsi=psipsi+psi(l,k)*psi(l,k)               
                psipsid=psipsid+psi(l,k)*psid(l,k)
                !
                err=psid(l,k)-psi(l,k)
                err2=err2+err*err
            enddo ! l
            enddo
        enddo ! k

        ! estimate error reduction factor (dominance ratio if no extrapolation)
        if(err2d .ne. 0) domr=sqrt(err2/err2d)
        
        eigvd=eigv
! update eigenvalue

        if(.not.iftran) then
            eigv=eigv*psipsi/psipsid
            reigv=1/eigv        
        endif
        
        if(iout.eq.1)  eigvfirst=eigv  ! 2013_07_30 . scb          

! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
        negative=0
        do k=1,nz
        do l=1,nxy
        do m=1,ng
        if(phif(m,l,k) .lt. 0) then
            negative=negative+1
        endif
        enddo
        enddo
        enddo

        if(ng.eq.ng2 .and. negative .ne. 0 .and. negative .ne. nxy*nz*ng) then
            write(mesg,'(a,i6,"/",i6)') 'NEGATIVE FLUX : ', negative, nxy*nz*ng
            call message(true,true,mesg)
        endif
! added end
        
        ! update residual error
        call residualmg(phif,psi,resid)
        if(iout.eq.1) resid0 = resid
        relresid = resid/resid0 
        errl2=sqrt(err2/psipsid)
        erreig=abs(eigv-eigvd)
!        write(mesg,100) iout+ibeg,eigv,erreig,errl2,domr,resid,relresid, sum(phif)
        write(mesg,100) iout+ibeg,eigv,erreig,errl2,domr,resid,relresid,sum(phif),negative ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
        call message(true,true,mesg)

        !if(errl2.lt.epsl2) exit
        if(.not. ifrecsp3 .and. errl2.lt.epsl2) exit   ! 2013_07_30 . scb
    enddo !of iout 
    iout=min(iout,ncmfd)
    ibeg=ibeg+iout
    
    return
!100 format(i5,f10.7,1p,2e12.5,1p,3e10.3,1p,e14.5)
100 format(i5,f10.7,1p,2e12.5,1p,3e10.3,1p,e10.3,i10) ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
end subroutine