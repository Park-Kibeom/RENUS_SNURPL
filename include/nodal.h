! file name = nodal.h
      integer(4) :: funcinit, funcinit2, funcdrive
      
      character*8 :: nodkern
      logical :: is1n,nodalonly
      integer :: nmaxswp, nlupd,nsrcexp
      real    :: nlswp,epsnodal,underrelx


      common /nodal_funcp/funcinit,funcinit2,funcdrive
      common /nodal_settingc/ nodkern
      common /nodal_settingl/ is1n,nodalonly
      common /nodal_settingi/ nmaxswp, nlupd,nsrcexp
      common /nodal_settingf/ nlswp,epsnodal,underrelx  
  
      
      real, pointer, dimension(:,:,:,:) :: curil,  & !(ndir,l,k,m)
                                           curir,  & !(ndir,l,k,m)
                                           curol,  & !(ndir,l,k,m)
                                           curor     !(ndir,l,k,m)
     
      real, pointer, dimension(:,:,:) :: wr,  & !(ls,k,m)
                                         wz     !(l,ks,m)

      common /nodal_cur/ curil,curir,curol,curor,trlcff0,trlcff1,trlcff2

      common /nodal_weighting/ wr, wz

      integer, parameter :: SELF=1, NEIB=2
