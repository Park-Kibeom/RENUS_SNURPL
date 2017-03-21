    subroutine Line_upper(aa)
!
      use param
!
! convert lowercase string to uppercase
      parameter (INDXA=97,IDNXZ=122)
      character aa*(*)
!
      lenaa=len_trim(aa)
      i=1
      do while (i.le.lenaa)
         ia=ichar(aa(i:i))
         if(ia.ge.INDXA) aa(i:i)=char(ia-32)
         i=i+1
         if(i.gt.lenaa) return
      enddo
    
    end subroutine