  subroutine Decide_delta(ibranch,branch,nbranch,delta, dXS1,dXS)    ! rearrange branch(s) including BASECOND and calculate gradient
    USE PARAM
    USE MASTERXSL
    use readN2A
  
  
    implicit none
    
    integer :: i, j, k, s,g, id
    integer, intent(in) :: ibranch, nbranch
    real, intent(in) :: branch(nbranch), delta, dXS1(nbranch) 
    real, intent(out) :: dXS
    real :: temp(nbranch+1), dtemp(nbranch), tmp

    temp=0._8; dtemp=0._8; tmp=0._8
    
    if (nbranch<=1) then
      id=1
      tmp=dxs1(1)
    else if (nbranch>1) then
      ! rearrangement including base condition(=0) => getting temp(:)
      do s=1, nbranch
        temp(s)=branch(s)      
        if (branch(s)>=0._8) then
          j=s
          temp(j)=0._8
          k=j          ! k=base condition
          do j=s+1, nbranch+1
            temp(j)=branch(j-1)
          end do
          exit
        end if
        if(branch(nbranch)<0._8) then
          do j=1, nbranch
            temp(j)=branch(j)
          end do
          temp(nbranch+1)=0._8
          k=nbranch+1
          exit
        end if
      end do

      ! getting interval of branch
      do i=1, nbranch
        dtemp(i) = temp(i+1)-temp(i)
      end do
      
      
      ! decide number of ID and recalculate gradient of dppm
      if (delta<0._8) then
        do s=2, nbranch
          if(delta<temp(s)) then
            id=s-1
            exit 
          else if(delta>=temp(nbranch)) then
            id=nbranch
            exit 
          end if
        end do
        
        if ((k-id)/=1) then
          do i=k-1, id+1, -1   ! k=base condition
            tmp = tmp - dXS1(i)*(temp(i+1)-temp(i))
          end do
        end if
        tmp = tmp - dXS1(id)*(temp(id+1)-delta) 

      else if (delta>=0._8) then
        do s=2, nbranch
          if(delta<=temp(s)) then
            id=s
            exit
          else if(delta>temp(nbranch)) then
            id=nbranch
            exit
          end if
        end do 
        
        if ((id-k)/=1) then
          do i=k, id-1   ! k=base condition
            tmp = tmp + dXS1(i)*(temp(i+1)-temp(i))
          end do
          tmp = tmp + dXS1(id)*(delta-temp(id))
        else if ((id-k)==1) then
          tmp = tmp + dXS1(id-1)*(delta-temp(id-1))
        end if       
      end if
    end if    
    
    
    
    
  
    if (delta/=0._8) then
      dXS=tmp/delta 
    else if (delta==0._8) then
      dXS=0._8
    end if
    
  end subroutine