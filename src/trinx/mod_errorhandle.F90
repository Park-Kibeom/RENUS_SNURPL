! added in ARTOS ver. 0.2 ( for TRINKX ). 2012_07_05 by SCB  
    module ErrorHandle
      use param
      use BasicParam

      contains

!=======================================================================================!
      subroutine write_error( string )

        character(*), intent(in) :: string

        write( nerror, * ) string, ' error'
        stop

      end subroutine write_error

!=======================================================================================!
      subroutine write_error2( string, key_word )

        character(*), intent(in) :: string
        character(*), intent(in) :: key_word

        write( nerror, * ) string, ' error : ', key_word
        stop

      end subroutine write_error2

!=======================================================================================!
      subroutine write_error3( string, key_word )

        character(*), intent(in) :: string
        integer(NBI), intent(in) :: key_word

        write( nerror, * ) string, ' error : ', key_word
        stop

      end subroutine write_error3

!=======================================================================================!
      subroutine write_error4( string, key_word )

        character(*), intent(in) :: string
        real(NBF), intent(in)    :: key_word

        write( nerror, * ) string, ' error : ', key_word
        stop

      end subroutine write_error4

!=======================================================================================!
      subroutine check_word( string1, string2, check )
        character(*), intent(in)  :: string1, string2
        logical,      intent(out) :: check


        if ( trim(string1) == trim(string2) ) then
        check = .TRUE.
        else
        check = .FALSE.
        endif

      end subroutine check_word

    end module ErrorHandle