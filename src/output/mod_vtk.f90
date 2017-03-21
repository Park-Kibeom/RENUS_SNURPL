! 2014_01_08 . scb for visualization . .   
  module vtk
    
    logical :: flagvtk=.false.
    integer :: ivtkfreq=0,ivtkstep=0
    integer :: np, nc, nxyz
    
    real,pointer :: px(:), py(:), pz(:), pr(:,:)
    
    integer,pointer :: iovtk(:)
    integer,pointer,dimension(:,:,:) :: indp, indc
    
    type ctop
      sequence
      integer :: ix, iy, iz
      integer,Allocatable :: ippp(:)
    end type
    
    type(ctop),pointer :: ncinfo(:)
    
  end module
! added end  