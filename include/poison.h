type poisondata
  sequence
  integer:: dim
  real(8) :: lxe
  real(8) :: lio
  real(8), dimension(:), pointer :: yio
  real(8), dimension(:), pointer :: yxe
!  yio_,yxe_,arate_,frate_,&
!        dim_, data, io0_
end type poisondata
