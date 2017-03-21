    subroutine updppm(eigv,targeteigv)
      use param
!      
      include 'global.h'
      include 'itrcntl.h'
      include 'xsec.h'
      include 'srchppm.h'
      
      logical :: first=TRUE
      real    :: eigvd,ppmd,ppmn,slope
      save       eigvd,ppmd,ppmn,first
      
      if(first) then
         ppmd=ppm
         ppm=ppm+100

         eigvd=eigv
         eigv=eigv-0.01

         first=FALSE
      else
         if( abs(eigvd-eigv) .lt. 1.0e-7 ) then   !r7m0
            slope = -10000.0
         else
            slope = (ppmd-ppm)/(eigvd-eigv)
            if( slope .lt. -20000.0 ) slope = -20000.0
            if( slope .gt. -500.0 )   slope = -500.0
         endif
         ppmd=ppm
         ppm=(targeteigv-eigv)*slope+ppm
         eigvd=eigv
         eigv=targeteigv
      endif
           
      write(mesg, '(a25,2f16.3)') 'UPDATED BORON(PPM / Diff):', ppm,ppm-ppmd
      call message(true,true,mesg)
    end subroutine