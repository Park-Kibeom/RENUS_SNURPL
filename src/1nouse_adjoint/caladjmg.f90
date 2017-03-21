subroutine caladjmg(nout)
    use const
    use allocs
    use geom,   only : ng, nxy, nx, ny, nz,   &
                       nxs, nxe,              &
                       nys, nye,              &
                       nxsf , nxef,           &
                       kfbeg, kfend,          &
                       nodel, volnode

    use xsec,   only : xsnff     ,             &
                       xschif    ,             &
                       xssf

    use sfam,   only : phif
    
    use sfam_cntl, only : ninmax,             &
                          epsbicgmg
    
                       
    use cmfdmg, only : ccrf,                  &
                       cczf,                  &
                       srcf

    use bicgmg, only : initbicg1g ,           &
                       solbicg1g
    use adjoint
    use chebyacc, only : runcheby
                           
    implicit none

    integer       ::  nout
    
    integer       ::  i,j,k,le,l,m,kp1,ls,iout,iin,me

    real          ::  temp    ,         &
                      r20     ,         &
                      r2      ,         &
                      psiad(ng),        & 
                      fs      ,         &
                      ss      ,         &
                      psipsid ,         &
                      psipsi  ,         &
                      vol     ,         &
                      err2d   ,         &
                      err2    ,         &
                      domr

    character*80  ::  mesg


    err2 = 1.

    do k=1,nz
    do l=1,nxy
        fs=sum(xschif(:,l,k)*phiaf(:,l,k))*volnode(l,k)
        psiaf(:,l,k)=xsnff(:,l,k)*fs
    enddo
    enddo
        
    do iout=1,nout
        do m=ng,1,-1
            do k=1,nz
            do l=1,nxy
                srcf(m,l,k)=reigv*psiaf(m,l,k)
                ss = 0
                do me=1,ng
                    ss=ss+phiaf(me,l,k)*xssf(m,me,l,k)
                enddo
                srcf(m,l,k)=srcf(m,l,k)+ss*volnode(l,k)
            enddo
            enddo      

            call initbicg1g(m,phiaf,srcf,r20)      

            do iin=1,ninmax
                call solbicg1g(m, r20, phiaf, r2)
                if(r2.lt.epsbicgmg) exit
            enddo
        enddo
        
        psipsid = 0
        psipsi  = 0
        err2d   = err2
        err2    = 0
        do k=kfbeg,kfend
        do j=1,ny
        do i=nxsf(j),nxef(j)
            l=nodel(i,j)
            
            psiad(:)=psiaf(:,l,k)

            fs=sum(xschif(:,l,k)*phiaf(:,l,k))*volnode(l,k)
            psiaf(:,l,k)=xsnff(:,l,k)*fs
            
            err2=err2+sum((psiad(:)-psiaf(:,l,k))**2)
            
            psipsi  = 0
            psipsid = 0
            do m=1,ng
                psipsi=psipsi+psiaf(m,l,k)*psiaf(m,l,k)               
                psipsid=psipsid+psiaf(m,l,k)*psiad(m)
            enddo
        enddo
        enddo
        enddo

        domr=sqrt(err2/err2d)

        eigvd=eigv        
        eigv=eigv*psipsi/psipsid
        reigv=1/eigv      

        call runcheby(iout, domr, phiaf, psiaf)
                  
        erreig = abs(eigvd-eigv)
        if(erreig.le.epseig) exit    

        write(mesg,'(i4,f12.7,1p,e12.5)') iout, eigv, erreig
        call message(TRUE,TRUE,mesg)
    enddo
    
    print *, sum(phiaf(:,:,10))/sum(phiaf)
                           
end subroutine