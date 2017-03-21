    function calbchxs1st(ic,dm,ppm,tf,sigxs)
      use param
      include 'global.h'
      include 'moxbch.inc'
      real              :: calbchxs1st
      real, intent(in)  :: sigxs(ntotpnt)

      idmmax=1
      idmmin=1
      ddm=0
      if(npnt(IDM,ic).ne.0) then
        do ipnt=npnt(IDM,ic),1,-1
          if(dm .ge. (pnts(ipnt,IDM,ic)*(1-1e-6))) then
            idmmin=ipnt
            exit
          endif
        enddo

        if(idmmin.eq.npnt(IDM,ic)) then
          idmmin=idmmin-1       
          idmmax=npnt(IDM,ic)
        else
          idmmax=idmmin+1
        endif
        
        dmmin=pnts(idmmin,IDM,ic)
        dmmax=pnts(idmmax,IDM,ic)
        ddm=(dm-dmmin)/(dmmax-dmmin)
        if(ddm.le.1e-6) then
          ddm=0
          idmmax=idmmin
        elseif(ddm.le.(1+1e6) .and. ddm.ge.(1-1e-6)) then
          ddm=1
          idmmin=idmmax
        endif
      endif
      
      ippmmax=1
      ippmmin=1
      dppm=0
      if(npnt(IPPM,ic).ne.0) then
        do ipnt=npnt(IPPM,ic),1,-1
          if(ppm .ge. (pnts(ipnt,IPPM,ic)*(1-1e-6))) then
            ippmmin=ipnt
            exit 
          endif
        enddo
        
        if(ippmmin.eq.npnt(IPPM,ic)) then
          ippmmin=ippmmin-1       
          ippmmax=npnt(IPPM,ic)
        else
          ippmmax=ippmmin+1
        endif
        
        ppmmin=pnts(ippmmin,IPPM,ic)
        ppmmax=pnts(ippmmax,IPPM,ic)
        dppm=(ppm-ppmmin)/(ppmmax-ppmmin)
        if(dppm.le.1e-6) then
          dppm=0
          ippmmax=ippmmin
        elseif(dppm.le.(1+1e6) .and. dppm.ge.(1-1e-6)) then
          dppm=1
          ippmmin=ippmmax
        endif        
      endif

      itfmax=1
      itfmin=1
      dtf=0
      if(npnt(ITF,ic).ne.0) then
        do ipnt=npnt(ITF,ic),1,-1
          if(tf .ge. (pnts(ipnt,ITF,ic)*(1-1e-6))) then
            itfmin=ipnt
            exit
          endif
        enddo

        if(itfmin.eq.npnt(ITF,ic)) then
          itfmin=itfmin-1      
          itfmax=npnt(ITF,ic)
        else
          itfmax=itfmin+1
        endif
        
        tfmin=pnts(itfmin,ITF,ic)
        tfmax=pnts(itfmax,ITF,ic)
        dtf=(tf-tfmin)/(tfmax-tfmin)
      endif
            
      ilow=idxofxs(ic, idmmin,ippmmin,itfmin)
      iup=idxofxs(ic, idmmin,ippmmin,itfmax)
      v1=sigxs(ilow)*(1-dtf)+sigxs(iup)*dtf

      ilow=idxofxs(ic, idmmin,ippmmax,itfmin)
      iup=idxofxs(ic, idmmin,ippmmax,itfmax)      
      v2=sigxs(ilow)*(1-dtf)+sigxs(iup)*dtf
      
      ilow=idxofxs(ic, idmmax,ippmmin,itfmin)
      iup=idxofxs(ic, idmmax,ippmmin,itfmax)
      v3=sigxs(ilow)*(1-dtf)+sigxs(iup)*dtf

      ilow=idxofxs(ic, idmmax,ippmmax,itfmin)
      iup=idxofxs(ic, idmmax,ippmmax,itfmax)
      v4=sigxs(ilow)*(1-dtf)+sigxs(iup)*dtf
      
      b5=v1*(1-dppm)+v2*dppm
      b6=v3*(1-dppm)+v4*dppm
      
      calbchxs1st=b5*(1-ddm)+b6*ddm
      
    end function      
!=======================================================================================!
    function calbchxs2nd(ic,dm,ppm,tf,sigxs)
      use param
      include 'global.h'
      include 'moxbch.inc'
      real, intent(in)  :: sigxs(ntotpnt)
      
      real              :: del(2),x(3),y(3),cff(3)
      real              :: f9(3,3), f3(3)
      real              :: calbchxs2nd
      real              :: lagrange

      if(npnt(IDM,ic).ne.0) then
        x(:)=pnts(1:3,ITF,ic)
        do iidm=1,3
          do iippm=1,3
            y(1)=sigxs(idxofxs(ic,iidm,iippm,1))
            y(2)=sigxs(idxofxs(ic,iidm,iippm,2))
            y(3)=sigxs(idxofxs(ic,iidm,iippm,3))
            f9(iippm,iidm) = lagrange(x,y,tf)
          enddo
        enddo
        x(:)=pnts(1:3,IPPM,ic)
        do iidm=1,3
            y(:)=f9(:,iidm)
            f3(iidm) = lagrange(x,y,ppm)
        enddo
        calbchxs2nd = lagrange(pnts(1:3,IDM,ic),f3,dm)
      else ! using ppm only
        x(:)=pnts(1:3,IPPM,ic)
        y(1)=sigxs(idxofxs(ic,1,1,1))
        y(2)=sigxs(idxofxs(ic,1,2,1))
        y(3)=sigxs(idxofxs(ic,1,3,1))

        calbchxs2nd = lagrange(x,y,ppm)
      endif      

      return
    end function
!=======================================================================================!
    function lagrange(x,y,xpos)

      real          :: lagrange
      real          :: x(3), y(3)
      real          :: xpos
      real          :: del(2), cff(3)
      
      del(:)=x(2:3)-x(1:2)

      del1 = 1/(del(1)*(del(1)+del(2)))
      del2 = 1/(del(1)*del(2))
      del3 = 1/((del(1)+del(2))*del(2))

      cff(1) = del1*y(1)-del2*y(2)+del3*y(3)
      cff(2) = -(x(2)+x(3))*del1*y(1)+(x(1)+x(3))*del2*y(2)-(x(1)+x(2))*del3*y(3)
      cff(3) = (x(2)*x(3))*del1*y(1)-(x(1)*x(3))*del2*y(2)+(x(1)*x(2))*del3*y(3)
      
      lagrange = cff(1)*xpos*xpos+cff(2)*xpos+cff(3)
      
    end function
!=======================================================================================!
    function idxofxs(ic,idm,ippm,itf)
      use param
      include 'global.h'    

      integer idxofxs
      
      if(idm.eq.0) then
        idxofxs=ippm
      else
        idxofxs=(itf-1)*9
        idxofxs=idxofxs+(ippm-1)*3
        idxofxs=idxofxs+idm
      endif
    end function      
