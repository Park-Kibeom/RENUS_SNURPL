! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB
    subroutine trimw6( string1, string2 )
      use param
      use BasicParam

      character, dimension(ISTN) :: string1, string2
      integer(NBI) :: i,j

      character(ISTN), parameter :: BLANK = '      '
      character,       parameter :: SPACE = ' '

      string1 = BLANK
      j=0
      do i=1,ISTN
        if(string2(i) .ne. SPACE) then
          j=j+1
          string1(j) = string2(i)
        endif
      enddo

    end subroutine trimw6
