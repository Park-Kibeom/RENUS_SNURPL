! added in ARTOS ver. 0.2 . 2012_07_24 by SCB   
   subroutine runsteady_hex_nodal(flag2g,                &
                           noutbegmg,  ninbegmg,  &
                           noutbeg2g,  ninbeg2g,  &
                           epsl2, epseig,     errl2)

      use const
      use sfam
      use timer
      use itrinfo
      use sfam_cntl,  only : nmaxcmfd2g                             
      use geom,       only : ng,nxy,nz,kfbeg,kfend,   &
                              volnode
      use xsec
      use tpen,       only : drivetpen, hflx, hflxf
      use tpen_sp3,   only : drivetpen_sp3, hflx2, hflxf2, hflxf2d
      use geomhex,    only : ifhexfdm, ifhexsp3
      use cmfd2g,     only : drivecmfd2g
      use cmfdhex2g_sp3,  only : setlshex2g_sp3, ilufachex_sp3,  &
                                 drivecmfd2g_sp3
       
      implicit none

      logical                 :: flag2g
      integer                 :: noutbegmg,  ninbegmg
      integer                 :: noutbeg2g,  ninbeg2g
      real                    :: epsl2,erreig,errl2,epseig
      integer                 :: iout,iup,m,l,k,mbeg,ms,iin,i,j
      real                    :: err2d,err2,err,ss,r20,r2,  &
                               !psipsid,psipsi,psid,domr,eigvd,  &
                               psipsid,psipsi,domr,eigvd,  &   ! 2013_05_14 . scb
                               resid,resid0,relresid, phim,domr2
      character               :: mesg*120  
      real                    :: tstart, tend
      real                    :: phiphi0, phiphi2, errphi0, errphi2
      integer                 :: negative

      real,pointer            :: psid1(:,:), psid2(:,:)
      real                    :: errd, errdd, avgdomr
      
      integer :: nout
      
      logical,save :: first=.true.  ! 2013_10_08 . scb
      
! 2013_10_08 . scb      
      if(first) then
        first=.false.
        write(mesg,'(a5,a10,1p,2a12,a15)') 'iout','eigv','erreig','errl2','negative flux'
        call message(true,true,mesg)        
      endif
! added end

      allocate(psid1(nxy,nz), psid2(nxy,nz))
        
      call cpu_time(tstart)
       
      if(ifhexsp3) then
        nout=1
      else
        nout=500
      endif
      
      avgdomr=0.d0              
      do iout=1,nout                            	      
         if(ifhexsp3) then
            do k=1,nz
               do l=1,nxy
                  do m=1,ng
                     hflxf2d(1,m,l,k)=hflxf2(1,m,l,k) 
                     hflxf2d(2,m,l,k)=hflxf2(2,m,l,k)
                  enddo
               enddo
            enddo 
            call drivetpen_sp3(FALSE)
         else
            call drivetpen(FALSE)            

           err2d=err2;err2=0.d0;
           psipsid=0.d0;psipsi=0.d0
         
           errd=0.d0
           errdd=0.d0
           avgdomr=0.d0
           do k=kfbeg,kfend
              do l=1,nxy
                 psid2(l,k)=psid1(l,k)
                 psid1(l,k)=psi(l,k)
         
                 psid(l,k)=psi(l,k)
                 psi(l,k)=0.d0
                 do m=1,ng
                    phim=hflxf(m,l,k)                              
                    psi(l,k)=psi(l,k)+xsnff(m,l,k)*phim
                 enddo
                 psi(l,k)=psi(l,k)*volnode(l,k)
         
                 psipsi=psipsi+psi(l,k)*psi(l,k)               
                 psipsid=psipsid+psi(l,k)*psid(l,k)
            
                 err=psid(l,k)-psi(l,k)
                 err2=err2+err*err
         
                 errd=errd+(psid1(l,k)-psi(l,k))*(psid1(l,k)-psi(l,k))
                 errdd=errdd+(psid2(l,k)-psid1(l,k))*(psid2(l,k)-psid1(l,k))
              enddo ! l
           enddo ! k
         
           ! estimate error reduction factor 
           if(err2d .ne. 0) domr=sqrt(err2/err2d)
           domr2=sqrt(errd/errdd)
           avgdomr=avgdomr+domr
         
           eigvd=eigv
           ! update eigenvalue
           eigv=eigv*psipsi/psipsid
           if(iout.le.5) eigv=1.d0
           reigv=1.d0/eigv 
         
           errl2=sqrt(err2/psipsid)
           erreig=abs(eigv-eigvd)/eigv

! 2013_10_08 . scb           
            do k=1,nz
               do l=1,nxy
                  do m=1,ng
                     if(ifhexsp3) then
                        phim=hflxf2(1,m,l,k)-2.*hflxf2(2,m,l,k)
                     else
                        phim=hflxf(m,l,k)                     
                     endif
                     phif(m,l,k)=phim
                  enddo
               enddo
            enddo

            negative=0
            do k=1,nz
               do l=1,nxy
                  do m=1,ng
                     if(phif(m,l,k) .lt. 0) then
                        negative=negative+1
                        !phif(m,l,k)=1.e-30  ! 2013_10_08 . scb. negative fixup
                     endif
                  enddo
               enddo
            enddo         
! added end            
          
           write(mesg,100) iout,eigv,erreig,errl2,negative   ! 2013_10_08 . scb
           call message(true,true,mesg)
         
           if(erreig.lt.epseig .and. errl2.lt.epsl2 .and. iout.gt.10) exit ! mhtr            
         endif         
      enddo 

      do k=1,nz
         do l=1,nxy
            do m=1,ng
               if(ifhexsp3) then
                  phim=hflxf2(1,m,l,k)-2.*hflxf2(2,m,l,k)
               else
                  phim=hflxf(m,l,k)                     
               endif
               phif(m,l,k)=phim
            enddo
         enddo
      enddo

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
      print *, negative

      call cpu_time(tend)
      print '(/a8,f8.4,a4/)', ' Time : ', tend-tstart, 'sec'
           
100   format(i5,f10.7,1p,2e12.5,i15)   ! 2013_10_08 . scb
      
		return 
   end subroutine
