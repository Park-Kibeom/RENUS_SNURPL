    subroutine updmoxbchxs(iffdbk, calbchxs)
      use param
      use sfam,     only : fphi
      use xsec,     only : colxs
      use xsec,     only : xstfsp3, xsrsp3   ! 2015_08_03 . scb
      include 'global.h'
      include 'files.h'
      include 'geom.h'
      include 'xsec.h'
      include 'srchppm.h'
      include 'moxbch.inc'
      include 'thgeom.inc'
      include 'thfdbk.inc'
      include 'frzfdbk.inc'
      include 'thop.inc'
      logical,save :: first=TRUE
      logical, intent(in) :: iffdbk
      
! update xsec
      do k=1,nz
        ka=ktoka(k)
        if(first) then
          do la=1,nxya
            iat=iassytyp(la)
            icomp=icompnumz(ka,iat)
!            nfmsqrt=sqrt(latol(la)%nfm+0.0)   ! 15_11_16 . scb
            do li=1,latol(la)%nfm
              l=latol(la)%fm(li)     
              
              xsmax(l,k)=2   ! 2015_08_03 . scb
              do m=1,ng         
                xschif(m,l,k)=sigchi(m,icomp)
                xsadf(:,m,l,k)=sigadf(:,m,icomp)
!  1      2
!   ------
!  |      |
!  |      |   : cdf index
!   ------
!  3      4
                xscdf(:,m,l,k)=sigmdf(:,m,icomp)              

! 15_11_16 . scb
                !if(li.eq.1) xscdf(1,m,l,k)=sigcdf(1,m,icomp)
                !if(li.eq.nfmsqrt) xscdf(2,m,l,k)=sigcdf(1,m,icomp)
                !if(li.eq.(latol(la)%nfm-nfmsqrt+1)) xscdf(3,m,l,k)=sigcdf(1,m,icomp)
                !if(li.eq.latol(la)%nfm) xscdf(4,m,l,k)=sigcdf(1,m,icomp)
! delete end
                
                xssfs(m,l,k)=sigsms(m,icomp)
                xssfe(m,l,k)=sigsme(m,icomp)
!                xscdf(:,m,l,k) = xsadf(:,m,l,k)
              enddo ! m
            enddo ! li
          enddo ! la
        endif ! first
                
        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          ia=latoia(la)
          ja=latoja(la)
          icomp=icompnumz(ka,iat)
          irodtyp1=irodtyp(la)
          
          if(iffdbk) then
            kth=ktokth(k)
            lchan=ltochan(l)        
            dm=dcool(kth,lchan)  
            tf=tdopl(kth,lchan)
          else
            dm=unifdm
            tf=sqrt(uniftf)
          endif
          
          do m=1,ng
            xstrf1=calbchxs(icomp,dm,ppm,tf, sigbchtr(:,m,icomp))
            xsaf(m,l,k)=calbchxs(icomp,dm,ppm,tf, sigbcha(:,m,icomp))
            xsnff(m,l,k)=calbchxs(icomp,dm,ppm,tf, sigbchnf(:,m,icomp))
            xskpf(m,l,k)=calbchxs(icomp,dm,ppm,tf, sigbchkf(:,m,icomp))
            xsadf(:,m,l,k)=calbchxs(icomp,dm,ppm,tf, sigbchadf(:,m,icomp))
            
            if(first) then
              xschif(m,l,k)=sigchi(m,icomp)
              xssfs(m,l,k)=1
              xssfe(m,l,k)=ng
            endif
            do ms=xssfs(m,l,k),xssfe(m,l,k)
              xssf(ms,m,l,k)=calbchxs(icomp,dm,ppm,tf,sigbchs(:,ms,m,icomp))
            enddo

! control rod            
            if(irodtyp1 .ne. 0 .and. iffuelc(icomp)) then
              if(rodfrac(k,irodtyp1) .ne. 0) then
                !icrod=icomp+20 ! sucks code.
                icrod=icomp+28 !  2015_08_03 . scb modified
                rf=rodfrac(k,irodtyp1)
                rf1=1-rf
                
                xstrfrod=calbchxs(icrod,dm,ppm,tf, sigbchtr(:,m,icrod))
                xstrf1=rf1*xstrf1+rf*xstrfrod

                xsafrod=calbchxs(icrod,dm,ppm,tf, sigbcha(:,m,icrod))
                xsaf(m,l,k)=rf1*xsaf(m,l,k)+rf*xsafrod
                
                xsnffrod=calbchxs(icrod,dm,ppm,tf, sigbchnf(:,m,icrod))
                xsnff(m,l,k)=rf1*xsnff(m,l,k)+rf*xsnffrod
                
                xskpfrod=calbchxs(icrod,dm,ppm,tf, sigbchkf(:,m,icrod))
                xskpf(m,l,k)=rf1*xskpf(m,l,k)+rf*xskpfrod
                
                xsadfrod=calbchxs(icrod,dm,ppm,tf, sigbchadf(:,m,icrod))
                xsadf(:,m,l,k)=rf1*xsadf(:,m,l,k)+rf*xsadfrod
                
                do ms=xssfs(m,l,k),xssfe(m,l,k)
                  xssmrod=calbchxs(icrod,dm,ppm,tf,sigbchs(:,ms,m,icrod))
                  xssf(ms,m,l,k)=rf1*xssf(ms,m,l,k)+rf*xssmrod
                enddo
              endif !rodfrac(k,irodtyp1) .ne. 0
            endif ! irodtyp1 .ne. 0 

            xsdf(m,l,k)=1.d0/(3.d0*xstrf1)            
            xstrf(m,l,k)=xstrf1   ! 2015_08_03 . scb            
            
	          xbeta(m,l,k,XDIR)=xsdf(m,l,k)*rhmesh(XDIR,l,k)
      	    xbeta(m,l,k,YDIR)=xsdf(m,l,k)*rhmesh(YDIR,l,k)
	          xbeta(m,l,k,ZDIR)=xsdf(m,l,k)*rhmesh(ZDIR,l,k)
! 2015_08_03 . scb            
            if(ifsp3) then
              if(.not. usesigt_p3 .or. xstfsp3(m,l,k) .lt. 1.e-15) xstfsp3(m,l,k)=(xstrf(m,l,k)-xsaf(m,l,k))/xsrsp3(m,1) + xsaf(m,l,k)  ! 2015_08_03 . scb
              xsdf2(m,l,k)=3.d0/(7.d0*xstfsp3(m,l,k))  ! 2015_08_03 . scb
	            xbeta2(m,l,k,XDIR)=xsdf2(m,l,k)*rhmesh(XDIR,l,k)
              xbeta2(m,l,k,YDIR)=xsdf2(m,l,k)*rhmesh(YDIR,l,k)
	            xbeta2(m,l,k,ZDIR)=xsdf2(m,l,k)*rhmesh(ZDIR,l,k)
            endif            
! added end            
          enddo !m
          
!define GENXSEC
#ifdef GENXSEC    
          write(filename,'("comp",i2.2,".xsc")'), icomp
          open(1000,file=filename,status='unknown')
          write(1000, *) ' XSEC_ORDER sigtr siga signf sigkf sigchi'
          write(1000, *) ' base'
          do m=1,ng
            write(1000, '(1x,1p,i1,1x,5e12.5)') m, 1/(3*xsdf(m,l,k)), xsaf(m,l,k),xsnff(m,l,k),xskpf(m,l,k),xschif(m,l,k)
          enddo
          
          do m=1,ng
            write(1000, '(1x,i1,$)') m
            do md=1,ng
              write(1000, '(1x,1p,e12.5,$)') xssf(m,md,l,k)
            enddo
            write(1000,*)
          enddo
          
          write(1000, *) ' adf'
          do m=1,ng
            write(1000, '(1x,i1,$)') m
            write(1000, '(1x,1p,4e12.5)') xsadf(:,m,l,k)
          enddo
          close(1000)
#endif
          do m=1,ng
            xssf(m,m,l,k)=0 !remove self-scattering
          enddo
        enddo
      enddo
      
! 15_11_16 . scb
      if(rect) then
        do k=1,nz
          do la=1,nxya
            iat=iassytyp(la)
            icomp=icompnumz(ka,iat)
            nmx=nmeshx(latoia(la));   nmx1=nmx-1
            nmy=nmeshy(latoja(la));   nmy1=nmy-1
            do ly=1,nmy;  do lx=1,nmx1
              li1=nmx*(ly-1)+lx;     li2=nmx*(ly-1)+lx+1
              l1=latol(la)%fm(li1);  l2=latol(la)%fm(li2)
              do m=1,ng
                xsadf(2,m,l1,k)=1. ;  xsadf(1,m,l2,k)=1.
              enddo
            enddo; enddo
                    
            do ly=1,nmy1;  do lx=1,nmx
              li1=nmx*(ly-1)+lx;     li2=nmx*ly+lx
              l1=latol(la)%fm(li1);  l2=latol(la)%fm(li2)
              do m=1,ng
                xsadf(4,m,l1,k)=1. ;  xsadf(3,m,l2,k)=1.
              enddo
            enddo; enddo
          enddo
        enddo
      endif      
! add end

      do k=1,nz
        do l=1,nxy
          do m=1,ng
            summ=0
            do md=1,ng
              if(m .ge. xssfs(md,l,k) .and. m.le.xssfe(md,l,k)) then
                summ=summ+xssf(m,md,l,k)
              endif
            enddo
            xstf(m,l,k)=xsaf(m,l,k)+summ
          enddo
        enddo
      enddo      
      
      call colxs(fphi)
      first=false

      return
    end subroutine