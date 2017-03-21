subroutine Sol1DSENM(ng, eigv, iouter)
  use senm1d
  use senmop
  use directsenm, only : setsenmmat, UpdPsin, UpdSrcq, UpdPsol, SetRhs, UpdFlx, UpdCurrent, ChkPsiShp
  use sp3senm, only : jnet0, jnet2, lkg0, lkg2, lkg0d, lkg2d, avglkg
  use itrcntl, only :noutsenm, epsshp
  use sp3senm,  only : nmax  ! 2013_07_18 . scb
  implicit none
  
  integer, intent(in) :: ng, iouter
  real(8), intent(in) :: eigv
  !integer :: ndummy
  
  integer :: iter, m, matsize
  real(8) :: errsum
  
  integer :: ig, ih, iz   ! 2015_01_14 . scb  

! 2013_07_18 . scb for debugging
  !integer :: ig,in,i,j
  !real,pointer,dimension(:,:,:,:) :: jnetdbg, jneterr
  !real :: jnetrms, jnetmax
  !logical,save :: first=.true.
  !integer,save :: istep=0
  !
  !allocate(jnetdbg(2,2,nmax,ng))
  !allocate(jneterr(2,2,nmax,ng))
  !if(first) then
  !  first=.false.
  !  open(718,file='E:\myProject\SP3SENNONLY\run\IAEA\current',status='old')
  !endif
! added end  
    
  call SetSimTrans(ng, qt, ksq, s)
  call SetSENMmat(ng, ns, ne, bc, ksq, s)
  call UpdPsin(ng, ns, ne, phin, psin, psidn) ! 2013/01/28 added
  matsize=2*(ne-ns+1)
  !do iter=1,noutsenm
  !do iter=1,5
  do iter=1,10
    do m=1,ng
      !call UpdSrcq(m, ng, ns, ne, eigv, phin, psin, tlkg, qt, mqh)
      !call UpdSrcq(m, ng, ns, ne, eigv, phin, psin, tlkg, tspc, qt, mqh)  ! 2013_08_12 . scb
      call UpdSrcq(m, ng, ns, ne, eigv, phin, psin, tlkg, tspc, tdmd, qt, mqh)  ! 2014_10_06 . scb
      call UpdPsol(m, ng, ns, ne, ksq, mqh, mpsol)
      call SetRhs(m, ng, ns, ne, bc, ksq, s, mpsol, rhs)
      call SolSenmMatLU(m, matsize, rhs, coeff)
      call UpdFlx(m, ns, ne, ksq, s, mpsol, coeff, phin, mA, mB)
    end do
    
    call UpdPsin(ng, ns, ne, phin, psin, psidn)
    call ChkPsiShp(ng, ns, ne, psin, psidn, errsum)
    
    if (errsum.ne.0.and.errsum.lt.1.e-4) then
   !  print*, iter, errsum
      exit   
    end if
  end do
  
  call UpdCurrent(ng, ns, ne, bc, ksq, s, mA, mB, mpsol, jnet1d)
! 2013_07_18 . scb for debugging  
  !jnetrms=0.d0
  !do ig=1,ng
  !  do in=1,nmax
  !    do j=1,2
  !      do i=1,2
  !        istep=istep+1
  !        read(718,'(e40.30)') jnetdbg(i,j,in,ig)
  !        jneterr(i,j,in,ig)=abs(jnet1d(i,j,in,ig)-jnetdbg(i,j,in,ig))/abs(jnet1d(i,j,in,ig))
  !        if(jneterr(i,j,in,ig).gt.jnetmax)  jnetmax=jneterr(i,j,in,ig)
  !      enddo
  !    enddo
  !  enddo
  !enddo
! added end  
  
  continue
  
  !if (iouter.lt.20) call UpdTlkg(jnet0, jnet2, lkg0, lkg2, lkg0d, lkg2d, avglkg)
  
!  call UpdTlkg(jnet0, jnet2, lkg0, lkg2, lkg0d, lkg2d, avglkg)
end subroutine