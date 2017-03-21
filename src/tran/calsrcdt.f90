! 2014_09_23 . scb  
  subroutine calsrcdt(m,l,k,bdforder,spsrcdt)
  
    use xsec,       only : xstrf
    use tranxsec,   only : rvelof
    use bdf,        only : srcdt, srcdtbdf, bdfcoef
    use geom,       only : ng,nz,nxy
    
    integer :: bdforder
    integer :: m,l,k
    real    :: spsrcdt(3)
    integer :: i,idir
    
    real :: sol
    
    do idir=1,3
      sol=0.d0    
      do i=0,bdforder
        sol=sol + bdfcoef(i)*srcdtbdf(m,l,k,idir,i)
      enddo
      
      spsrcdt(idir)=sol*rvelof(m,l,k)/xstrf(m,l,k)
    enddo             
  
  end subroutine
! added end  