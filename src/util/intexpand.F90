    subroutine intexpand(n,ival,ifield,itarget,ntarget)
!
      use param
!
! expand list of integers containing minus signs and assigns ival
! to target array
      dimension ifield(*),itarget(*)
      iprev=ifield(1)
      do i=1,n
         if(ifield(i).gt.0) then
            if(ifield(i).le.ntarget) itarget(ifield(i))=ival
            iprev=ifield(i)
         elseif(ifield(i).lt.0) then
            do j=iprev+1,-ifield(i)
               if(j.le.ntarget) itarget(j)=ival
            enddo
         else
            mesg='Value 0 Not Allowed'
            call terminate(mesg)
         endif
      enddo
!
      return
      
    end subroutine