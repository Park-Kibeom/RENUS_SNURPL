    subroutine updmwt(eigv,targeteigv,count)
      ! Search critical power level 
      use param
!      
      include 'global.h'
      include 'itrcntl.h'
      include 'srchppm.h'
      include 'pow.h'
      
      integer :: count
      real    :: eigvd,ppmd,ppmn,slope
      save       eigvd,ppmd,ppmn

      select case (count)
      case (0)
         plevel = plevel
      case (1)
         ppmd=plevel
         plevel=plevel+1.E-1

         eigvd=eigv
         eigv=eigv-0.01
      case default
         slope = (ppmd-plevel)/(eigvd-eigv)
         ppmd=plevel
         plevel=(targeteigv-eigv)*slope+plevel
         eigvd=eigv
         eigv=targeteigv
      end select

    end subroutine