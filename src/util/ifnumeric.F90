    function ifnumeric(aline)
!
      use param
!      
      include 'cards.h'
      character aline*(*)
      logical :: ifnumeric
!
      ifnumeric=FALSE
      oneline=aline
      do i=1,mxncol
        iascii=ichar(sc(i))
        if(sc(i).ne.BLANK .and. iascii.ne.9 ) then  !determine if the first character is numeric
          if((iascii-48)*(iascii-57) .le. 0) ifnumeric=TRUE
          return
        endif
      enddo
!
      return
    end function