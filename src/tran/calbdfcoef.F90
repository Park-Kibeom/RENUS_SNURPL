! 2012_09_07 . scb
    subroutine calbdfcoef(order)
    
      use bdf
      
      integer :: order
      
      if ( order == 1 ) then
          bdfcoef(0) = 1./deltmarray(1)
          bdfcoef(1) = -1./deltmarray(1)
      elseif ( order == 2 ) then
          bdfcoef(0) = 1./deltmarray(1) + 1./(deltmarray(1)+deltmarray(2))
          bdfcoef(1) = -(deltmarray(1)+deltmarray(2))/(deltmarray(1)*deltmarray(2))
          bdfcoef(2) = deltmarray(1)/((deltmarray(1)+deltmarray(2))*deltmarray(2))
      elseif ( order == 3 ) then
          bdfcoef(0) = 1./deltmarray(1) + 1./(deltmarray(1)+deltmarray(2))   &
                     +1./(deltmarray(1)+deltmarray(2)+deltmarray(3))
          bdfcoef(1) = -(deltmarray(1)+deltmarray(2))*(deltmarray(1)+deltmarray(2)+deltmarray(3))   &
                     /(deltmarray(1)*deltmarray(2)*(deltmarray(2)+deltmarray(3)))
          bdfcoef(2) = deltmarray(1)*(deltmarray(1)+deltmarray(2)+deltmarray(3))   &
                     /((deltmarray(1)+deltmarray(2))*deltmarray(2)*deltmarray(3))
          bdfcoef(3) = -deltmarray(1)*(deltmarray(1)+deltmarray(2))   &
                     /((deltmarray(1)+deltmarray(2)+deltmarray(3))*(deltmarray(2)+deltmarray(3))*deltmarray(3))
      elseif ( order == 4 ) then
          bdfcoef(0) = 1./deltmarray(1) + 1./(deltmarray(1)+deltmarray(2))   &
                     +1./(deltmarray(1)+deltmarray(2)+deltmarray(3))   &
                     +1./(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4))
          bdfcoef(1) = -(deltmarray(1)+deltmarray(2))*(deltmarray(1)+deltmarray(2)+deltmarray(3))    &
                     *(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4))    &
                     /(deltmarray(1)*deltmarray(2)*(deltmarray(2)+deltmarray(3))*(deltmarray(2)+deltmarray(3)+deltmarray(4)))
          bdfcoef(2) = deltmarray(1)*(deltmarray(1)+deltmarray(2)+deltmarray(3))    &
                     *(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4))    &
                     /((deltmarray(1)+deltmarray(2))*deltmarray(2)*deltmarray(3)*(deltmarray(3)+deltmarray(4)))
          bdfcoef(3) = -deltmarray(1)*(deltmarray(1)+deltmarray(2))*(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4))    &
                     /((deltmarray(1)+deltmarray(2)+deltmarray(3))*(deltmarray(2)+deltmarray(3))*deltmarray(3)*deltmarray(4))
          bdfcoef(4) = deltmarray(1)*(deltmarray(1)+deltmarray(2))*(deltmarray(1)+deltmarray(2)+deltmarray(3))    &
                     /((deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4))*(deltmarray(2)+deltmarray(3)+deltmarray(4))    &
                     *(deltmarray(3)+deltmarray(4))*deltmarray(4))
      elseif ( order == 5 ) then
          bdfcoef(0) = 1./deltmarray(1) + 1./(deltmarray(1)+deltmarray(2))    &
                     +1./(deltmarray(1)+deltmarray(2)+deltmarray(3))     &
                     +1./(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4))    &
                     +1./(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4)+deltmarray(5))
          bdfcoef(1) = -(deltmarray(1)+deltmarray(2))*(deltmarray(1)+deltmarray(2)+deltmarray(3))    &
                     *(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4))    &
                     *(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4)+deltmarray(5))    &
                     /(deltmarray(1)*deltmarray(2)*(deltmarray(2)+deltmarray(3))    &
                     *(deltmarray(2)+deltmarray(3)+deltmarray(4))    &
                     *(deltmarray(2)+deltmarray(3)+deltmarray(4)+deltmarray(5)))
          bdfcoef(2) = deltmarray(1)*(deltmarray(1)+deltmarray(2)+deltmarray(3))    &
                     *(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4))    &
                     *(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4)+deltmarray(5))    &
                     /((deltmarray(1)+deltmarray(2))*deltmarray(2)*deltmarray(3)*(deltmarray(3)+deltmarray(4))    &
                     *(deltmarray(3)+deltmarray(4)+deltmarray(5))) 
          bdfcoef(3) = -deltmarray(1)*(deltmarray(1)+deltmarray(2))*(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4))    &
                     *(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4)+deltmarray(5))    &
                     /((deltmarray(1)+deltmarray(2)+deltmarray(3))*(deltmarray(2)+deltmarray(3))    &
                     *deltmarray(3)*deltmarray(4)*(deltmarray(4)+deltmarray(5)))
          bdfcoef(4) = deltmarray(1)*(deltmarray(1)+deltmarray(2))*(deltmarray(1)+deltmarray(2)+deltmarray(3))    &
                     *(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4)+deltmarray(5))    &
                     /((deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4))*(deltmarray(2)+deltmarray(3)+deltmarray(4))    &
                     *(deltmarray(3)+deltmarray(4))*deltmarray(4)*deltmarray(5))
          bdfcoef(5) = -deltmarray(1)*(deltmarray(1)+deltmarray(2))*(deltmarray(1)+deltmarray(2)+deltmarray(3))    &
                     *(deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4))    &
                     /((deltmarray(1)+deltmarray(2)+deltmarray(3)+deltmarray(4)+deltmarray(5))    &
                     *(deltmarray(2)+deltmarray(3)+deltmarray(4)+deltmarray(5))    &
                     *(deltmarray(3)+deltmarray(4)+deltmarray(5))*(deltmarray(4)+deltmarray(5))*deltmarray(5))
      endif
  
    end subroutine