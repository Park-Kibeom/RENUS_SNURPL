    subroutine getfn(aline,i,str)
!
      use param
!
! get i-th string on a line of strings
!
      character*(*) aline,str
!
      ncol=len(aline)
      icharjm1=ichar(aline(1:1))  !tab
      if(aline(1:1).eq.' ' .or. icharjm1.eq.9) then
         jf=0
      else
         jf=1
      endif
      jm1=1
      do j=2,ncol
         icharj=ichar(aline(j:j))
         if((aline(jm1:jm1).eq.' ' .or. icharjm1.eq.9) .and. (aline(j:j).ne.' ' .and. icharj.ne.9)) then
            jf=jf+1
            jbeg=j
         endif
         if(jf.eq.i) then
            if(aline(j:j).eq.' ' .or. icharj.eq.9) then
               str=aline(jbeg:j-1)
               return
            elseif(j.eq.ncol) then
               str=aline(jbeg:j)
               return
            endif
         endif
         jm1=j
         icharjm1=icharj
      enddo
!
      write(mesg,'("File name not present at",i3,"-th field.")') i 
      call terminate(mesg)
!
      return
!
    end subroutine