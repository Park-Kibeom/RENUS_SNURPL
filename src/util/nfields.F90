    function nfields(aline)
!
      use param
!
      include 'cards.h'
      include 'cntl.h'
!
      character*512 :: aline
      logical :: nonblankd,nonblank,multidata
      integer,save :: istep=0
!
      istep=istep+1
            
      nonblankd=.false.
      multidata=.false.
      oneline=aline
      ncol=len_trim(aline)
      n=0
      do i=ncol,1,-1
         if(sc(i).eq.BLANK .or. ichar(sc(i)).eq.9) then !ichar(tab)=9
            nonblank=FALSE
         else
            if(sc(i).eq.AST .and. .not.tracinp) then
               multidata=TRUE
               multidcol=i
            endif
            nonblank=TRUE
            if(sc(i).eq.BANG .or. (sc(i).EQ.AST .and. tracinp)) then
               n=-1
               nonblankd=TRUE
            endif
         endif
         if((.not.nonblankd.and.nonblank) .or. (nonblankd.and..not.nonblank)) then
            n=n+1
            if(multidata.and.(nonblankd.and. .not.nonblank)) then
               read(oneline(i+1:multidcol-1),*) nmult
               n=n+(nmult-1)*2
               multidata=FALSE
            endif
         endif   
         nonblankd=nonblank
      enddo
      !
      !write(903,*) istep
      !backspace(903)      
      
      if(mod(n,2).ne.0) then
        nfields=n/2+1
      else
        nfields=n/2
      endif
!
      return
      
    end function