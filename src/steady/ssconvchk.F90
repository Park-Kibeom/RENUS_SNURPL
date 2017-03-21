    logical function ssconvchk(eigv,errl2)
      use param
      
      include 'global.h'
      include 'itrcntl.h'
      include 'srchppm.h'
      
      logical,save :: doublecheck=FALSE
      
      ssconvchk=FALSE

      !if( srchppm .and. srchmwt) then
      if( srchppm .or. srchmwt) then    ! 2016.9.28 pkb
        if(abs(eigv-targetk).lt.epseig .and. errl2.lt.epsl2) then
          ssconvchk=TRUE
        endif
      else
        if(errl2.lt.epsl2) then
          ssconvchk=TRUE
        endif      
      endif
      
      if(.not.doublecheck .and. ssconvchk) then
        doublecheck=TRUE
        ssconvchk=FALSE
      elseif(.not.ssconvchk) then
        doublecheck=FALSE
      endif
      
    end function
