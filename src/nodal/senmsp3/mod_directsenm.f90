module directsenm
  use senmop
  use matop_sp3
  contains
  subroutine SetSENMmat(ng,ns,ne,bc,ksq,s)
    implicit none
    integer :: ng,ns,ne,bc(2)
    real(8),pointer,dimension(:,:,:) :: ksq, s
    
    real(8),dimension(4,ng) :: AL,AR,BL,BR
    integer :: m,i,icol,irow
    
    irow=1
    call BcCondition(AL, BL, ns, bc(1), 1, ng, ksq, s)
    icol=2*ns-1;   call UpdtMatElmt(ng, AL, irow, icol)   ! fill in the large linear system
    icol=2*ns;   call UpdtMatElmt(ng, BL, irow, icol)
    
    do i=ns,ne-1
      irow=irow+1
      call InterfaceMmtCondition(ng,AL,AR,BL,BR,i,ksq,s)
      icol=2*i-1;     call UpdtMatElmt(ng,AL,irow,icol)
      icol=2*i;       call UpdtMatElmt(ng,BL,irow,icol)
      icol=2*i+1;     call UpdtMatElmt(ng,AR,irow,icol)
      icol=2*i+2;     call UpdtMatElmt(ng,BR,irow,icol)
      
      irow=irow+1
      call InterfaceCurrentCondition(ng,AL,AR,BL,BR,i,ksq,s)
      icol=2*i-1;     call UpdtMatElmt(ng,AL,irow,icol)
      icol=2*i;       call UpdtMatElmt(ng,BL,irow,icol)
      icol=2*i+1;     call UpdtMatElmt(ng,AR,irow,icol)
      icol=2*i+2;     call UpdtMatElmt(ng,BR,irow,icol)
    end do
    irow=irow+1
    call BcCondition(AR,BR,ne,bc(2),2,ng,ksq, s)
    icol=2*ne-1;   call UpdtMatElmt(ng,AR,irow,icol)
    icol=2*ne;     call UpdtMatElmt(ng,BR,irow,icol)
    
    ! LU factorization
    do m=1,ng
      call Senmmat_LU(m, 2*(ne-ns+1))
    end do
    
  end subroutine
  
  subroutine InterfaceCurrentCondition(ng, A1, A2, B1, B2, i, ksq, s)
    use senm1d, only : xsdn, xsd2n, h
    use const, only : LEFT, RIGHT
    implicit none
    integer :: ng, i
    real(8),dimension(4,ng) :: A1, A2, B1, B2
    real(8),pointer :: ksq(:,:,:),s(:,:,:)
    
    integer :: m, is
    real(8) :: k(2,2)
    real(8),dimension(4,2) :: km, sm, shk, chk, mm
    
    do m=1,ng
      do is=1,2
        k(1,is)=sqrt(ksq(1,i+is-1,m))
        k(2,is)=sqrt(ksq(2,i+is-1,m))
        km(:,is)=kmat(k(1,is),k(2,is))
        sm(:,is)=s(:,i+is-1,m)
        shk(:,is)=sinhmat(k(1,is),k(2,is),1._8)
        chk(:,is)=coshmat(k(1,is),k(2,is),1._8)
        mm(:,is)=mmat(xsdn(i+is-1,m),xsd2n(i+is-1,m),h(i+is-1))
      end do
      A1(:,m)=0;   A2(:,m)=0
      B1(:,m)=0;   B2(:,m)=0
      A1(:,m)=multi_matop(mm(:,LEFT),sm(:,LEFT), km(:,LEFT), shk(:,LEFT));
      A2(:,m)=multi_matop(mm(:,RIGHT),sm(:,RIGHT), km(:,RIGHT), shk(:,RIGHT));
      B1(:,m)=multi_matop(mm(:,LEFT),sm(:,LEFT), km(:,LEFT), chk(:,LEFT));
      B2(:,m)=-multi_matop(mm(:,RIGHT),sm(:,RIGHT), km(:,RIGHT), chk(:,RIGHT));
    end do
  end subroutine
  
  subroutine InterfaceMmtCondition(ng, A1, A2, B1, B2, i, ksq, s)
    use const, only : LEFT, RIGHT
    use senm1d, only : xsdn, xsd2n, h
    use senm1d, only : xsadfn    ! 2015_08_05 . scb
    implicit none
    integer :: ng, i
    real(8),dimension(4,ng) :: A1, A2, B1, B2
    real(8),pointer :: ksq(:,:,:), s(:,:,:)
    
    real(8) :: k(2,2)
    real(8),dimension(4,2) :: km, sm, shk, chk, mm
    integer :: m, is
    integer :: l   ! 2015_08_05 . scb
    
    do m=1,ng
      do is=1,2
        l=i+is-1
        
        k(1,is)=sqrt(ksq(1,l,m))
        k(2,is)=sqrt(ksq(2,l,m))
        !k(1,is)=sqrt(ksq(1,i+is-1,m))
        !k(2,is)=sqrt(ksq(2,i+is-1,m))
        km(:,is)=kmat(k(1,is),k(2,is))
! 2015_08_05 . scb sp3 adf         
        !sm(:,is)=s(:,i+is-1,m)
        sm(1:2,is)=xsadfn(is,l,m)*s(1:2,l,m)
        sm(3:4,is)=1.d0
! added end        
        shk(:,is)=sinhmat(k(1,is),k(2,is),1._8)
        chk(:,is)=coshmat(k(1,is),k(2,is),1._8)
        !mm(:,is)=mmat(xsdn(i+is-1,m),xsd2n(i+is-1,m),h(i+is-1))    ! 2015_08_05 . scb commented because it is not used
      end do
      A1(:,m)=0;   A2(:,m)=0
      B1(:,m)=0;   B2(:,m)=0
      A1(:,m)=-submatop(sm(:,LEFT),chk(:,LEFT));  A2(:,m)=submatop(sm(:,RIGHT),chk(:,RIGHT))
      B1(:,m)=-submatop(sm(:,LEFT),shk(:,LEFT));  B2(:,m)=-submatop(sm(:,RIGHT),shk(:,RIGHT))      
    end do
    
  end subroutine
  
  subroutine UpdtMatElmt(ng, inElmt, irow, icol)
  ! inElmt : Element matrix of Large linear system
  ! irow : row number
  ! icol : column number
    use senm1d, only : nElmt,ElmtLoc,DiagIdx,nLowerElmt,LowerElmtIdx,Elmt
    implicit none
    integer :: ng,irow,icol
    real(8) :: inElmt(4,ng)
    
    integer :: m,j    
    
    do m=1,ng
      nElmt(irow,m)=nElmt(irow,m)+1   ! the number of elements in irow
      j=nElmt(irow,m)
      ElmtLoc(j,irow,m)=icol
      Elmt(:,j,irow,m)=inElmt(:,m)    ! linear system !
      if (irow.eq.icol) DiagIdx(irow,m)=j
      if (irow.gt.icol) then
        nLowerElmt(icol,m)=nLowerElmt(icol,m)+1
        j=nLowerElmt(icol,m)
        LowerElmtIdx(j,icol,m)=irow
      end if
    end do
    
  end subroutine
  
  subroutine BcCondition(A,B,i,bcc,idx,ng,ksq,s)
    use senm1d, only : xsdn, xsd2n, h
    use senmop, only : kmat,sinhmat,coshmat,mmat
    use senm1d, only : xsadfn   ! 2015_08_05 . scb for sp3senm adf
    implicit none
    integer,intent(in) :: i, bcc, idx, ng
    real(8),pointer :: ksq(:,:,:),s(:,:,:)
    real(8),intent(out) :: A(4,ng), B(4,ng)
    
    integer :: m
    !real(8) :: abd(4,0:1), k(2), x
    real(8) :: abd(4,0:2), k(2), x, abdm(4,0:2)   ! 2015_08_04 . scb
    real(8),dimension(4),target :: shk, chk, sm, km, mm
    ! 2015_08_04 . scb
    !data abd /4*0._8, 0.5_8, 0.625_8, -0.125_8, 0.625_8/
    data abd(1:4,0) /4*0._8   /
    data abd(1:4,1) /0.5_8, 0.625_8, -0.125_8, 0.625_8 /
    !data abd(1:4,2) /4*1.e30  /
    data abd(1:4,2) /1.e30, 0, 0, 1.e30  /   ! 2015_08_06 . scb
    
    integer :: idir
    ! added end
    
    ! 2015_08_05 . scb
    !x=(-1._8)**idx   ! = -1 at left boundary, +1 at right boundary    
    if(idx.eq.1) then
      x=-1._8
      idir=2
    elseif(idx.eq.2) then
      x=1._8
      idir=1
    endif    
    
    abdm=abd       
          
    !if(bcc.ne.0) then
    !  do m=1,ng
    !    abdm(1,1)=xsadfn(idir,i,m)*abd(1,1)
    !    abdm(3,1)=xsadfn(idir,i,m)*abd(3,1)
    !  enddo        
    !endif        
    ! added end

    do m=1,ng
      k(1)=sqrt(ksq(1,i,m)); k(2)=sqrt(ksq(2,i,m))   ! k=root(k2)
      km(:)=kmat(k(1),k(2))   ! K = [ k(1)   0  |  0  k(2) ]
      sm(:)=s(:,i,m)
      shk=sinhmat(k(1),k(2),1._8)    ! Sk = [ sinh(k(1))    0   |   0    sinh(k(2)) ]
      chk=coshmat(k(1),k(2),1._8)    ! Ck = [ cosh(k(1))    0   |   0    cosh(k(2)) ]
      mm=mmat(xsdn(i,m),xsd2n(i,m),h(i))   ! M = 2/h*[ D0   2*D0  |   0   D2 ]
      ! abd : albedo matrix
      ! abd(1:4,0) = [ 0 0 | 0 0 ]
      ! abd(1:4,1) = [ 0.5 0.625  |  -0.125 0.625 ] or [ 1/2  5/8  |  -1/8  5/8 ]
      ! A = -+ (M*S*K*Sk + ALB*S*Ck)
      ! B = ALB*S*Sk + M*S*K*Ck
! 2015_08_05 . scb      
      !A(:,m)=multi_matop(mm,sm,km,shk)+multi_matop(abd(:,bcc),sm,chk)
      !A(:,m)=x*A(:,m)
      !B(:,m)=multi_matop(abd(:,bcc),sm,shk)+multi_matop(mm,sm,km,chk)
      !      
      !if(bcc.ne.0) then
      !  abdm(1,1)=xsadfn(idir,i,m)*abd(1,1)
      !  abdm(3,1)=xsadfn(idir,i,m)*abd(3,1)
      !endif            
      A(:,m)=x*(multi_matop(mm,sm,km,shk)+multi_matop(abdm(:,bcc),sm,chk))
      B(:,m)=multi_matop(abdm(:,bcc),sm,shk)+multi_matop(mm,sm,km,chk)    
! added end      
    end do
  end subroutine

  subroutine UpdPsin(ng,ns,ne,phin,psin,psidn)
    use senm1d, only : xsnfn
    implicit none
    integer :: ng,ns,ne
    real(8),pointer :: phin(:,:,:,:),psin(:,:),psidn(:,:)

    integer :: i,j,m
    real(8) :: nufs
    
    do i=ns,ne
      do j=0,4
        psidn(j,i)=psin(j,i)
      end do
      !psin(1:4,i)=0.
      psin(1:4,i)=0.
    end do
    
    do i=ns,ne
      do m=1,ng
        nufs=xsnfn(i,m)
        do j=1,4
        !do j=1,4
          psin(j,i)=psin(j,i)+nufs*phin(1,j,i,m)
        end do
      end do
    end do
    
  end subroutine
  
  subroutine UpdPsilevel(ng,ns,ne,phin,psin,psidn)
    use senm1d, only : xsnfn
    implicit none
    integer :: ng,ns,ne
    real(8),pointer :: phin(:,:,:,:),psin(:,:),psidn(:,:)

    integer :: i,j,m
    real(8) :: nufs
    
    do i=ns,ne
      do j=0,4
        psidn(j,i)=psin(j,i)
      end do
      !psin(1:4,i)=0.
      psin(0:4,i)=0.
    end do
    
    do i=ns,ne
      do m=1,ng
        nufs=xsnfn(i,m)
        do j=0,4
        !do j=1,40
          psin(j,i)=psin(j,i)+nufs*phin(1,j,i,m)
        end do
      end do
    end do
    
  end subroutine
  
  function SetInterfaceCurrentRHS(mm,sm,c)
    implicit none
    integer,parameter :: l=1, r=2
    real(8) :: sm(4,2), mm(4,2)
    real(8) :: c(0:4,2,2)
    real(8) ::SetInterfaceCurrentRHS(2)
    real(8) :: mjs(2,2), pjs(2,2)
    
    mjs(:,l)=c(1,:,l)+3._8*c(2,:,l)+6._8*c(3,:,l)+10._8*c(4,:,l)
    pjs(:,l)=-submatvecop(submatop(mm(:,l),sm(:,l)),mjs(:,l))
    
    mjs(:,r)=c(1,:,r)-3._8*c(2,:,r)+6._8*c(3,:,r)-10._8*c(4,:,r)
    pjs(:,r)=-submatvecop(submatop(mm(:,r),sm(:,r)),mjs(:,r))
    
    SetInterfaceCurrentRHS=pjs(:,l)-pjs(:,r)
    
  end function
  
  !function SetInterfaceMmtRHS(sm,c)
  function SetInterfaceMmtRHS(sm,adf,c)    ! 2015_08_05 . scb for SP3 SENM ADF
    implicit none
    integer,parameter :: l=1, r=2
    real(8) :: sm(4,2)
    real(8) :: c(0:4,2,2)
    real(8) :: SetInterfaceMmtRHS(2)
    real(8) :: adf(2)        ! 2015_08_05 . scb for SP3 SENM ADF
    
    real(8) :: mphis(2,2), pphis(2,2)
    
    mphis(:,l)=c(0,:,l)+c(1,:,l)+c(2,:,l)+c(3,:,l)+c(4,:,l)
    pphis(:,l)=submatvecop(sm(:,l),mphis(:,l))
    pphis(1,l)=adf(l)*pphis(1,l)    ! 2015_08_05 . scb for SP3 SENM ADF
    
    mphis(:,r)=c(0,:,r)-c(1,:,r)+c(2,:,r)-c(3,:,r)+c(4,:,r)
    pphis(:,r)=submatvecop(sm(:,r),mphis(:,r))
    pphis(1,r)=adf(r)*pphis(1,r)    ! 2015_08_05 . scb for SP3 SENM ADF
    
    SetInterfaceMmtRHS=pphis(:,l)-pphis(:,r)
  end function
  
  !function SetBoundaryRHS(sm,mm,c,bc,idx)
  function SetBoundaryRHS(sm,mm,adf,c,bc,idx)      ! 2015_08_05 . scb for SP3 SENM ADF
    implicit none
    real(8) :: SetBoundaryRHS(2)
    real(8),dimension(4) :: sm, mm
    real(8) :: c(0:4,2)
    integer :: idx, bc
    !real(8) :: abd(4,0:1), lp_s(0:4,2), dlp_s(0:4,2)
    real(8) :: abd(4,0:2), lp_s(0:4,2), dlp_s(0:4,2)   ! 2015_08_04 . scb
    real(8) :: mphis(2), mjs(2), pphis(2), pjs(2)
    integer :: iod
    real(8) :: x,adf  ! 2015_08_04 . scb
    
    data abd(:,0) /4*0._8/
    data abd(:,1) /0.5_8, 0.625_8, -0.125_8, 0.625_8/
    !data abd(:,2) /4*1.e30/   ! 2015_08_04 . scb
    data abd(1:4,2) /1.e30, 0, 0, 1.e30  /   ! 2015_08_06 . scb
  
    data lp_s(:,1) /1._8, -1._8, 1._8, -1._8, 1._8/ ! Legendre Polynomial at x=-1
    data lp_s(:,2) /1._8, 1._8, 1._8, 1._8, 1._8/   ! Legendre Polynomial at x=1
    data dlp_s(:,1) /0._8, 1._8, -3._8, 6._8, -10._8/  ! Derivatives of Legendre Polynomial at x=-1
    data dlp_s(:,2) /0._8, 1._8, 3._8, 6._8, 10._8/   ! Derivatives of Legendre Polynomial at x=1
  
    mphis=0._8; mjs=0._8
    do iod=0,4
      mphis=mphis+(lp_s(iod,idx)*c(iod,:))   ! c(0:4, 1:2)
      mjs=mjs+(dlp_s(iod,idx)*c(iod,:))
    end do
    pphis=submatvecop(sm,mphis)
    pjs=-submatvecop(submatop(mm,sm),mjs)
    ! 2015_08_04 . scb
    !SetBoundaryRHS=pjs+(-1._8)**(idx+1)*submatvecop(abd(:,bc),pphis)    
    x=1._8
    if(bc.eq.1) then
      !abd(1,1)=adf*abd(1,1)
      !abd(3,1)=adf*abd(3,1)
    endif    
    if(idx.eq.2)  x=-1._8        
    SetBoundaryRHS=pjs+x*submatvecop(abd(:,bc),pphis)
    ! added end
  end function
 
  !subroutine UpdSrcq(m,ng,ns,ne,eigv,phin,psin,tlkg,qt,mqh)
  !subroutine UpdSrcq(m,ng,ns,ne,eigv,phin,psin,tlkg,tspc,qt,mqh)  ! 2013_08_12 .scb
  subroutine UpdSrcq(m,ng,ns,ne,eigv,phin,psin,tlkg,tspc,tdmd,qt,mqh)  ! 2014)10_06 .scb
    use senm1d, only : xssn, xschin
    implicit none
    integer :: m, ng, ns, ne
    real(8) :: eigv
    real(8),pointer :: phin(:,:,:,:), psin(:,:), tlkg(:,:,:,:), mqh(:,:,:,:), qt(:,:,:)
    real(8),pointer :: tspc(:,:,:,:)   ! 2013_08_12 . scb
    real(8),pointer :: tdmd(:,:,:,:)   ! 2014_10_06 . scb
    
    integer :: i, ms, iod
    !real(8) :: ss(0:4), qtmat(4), reigv, q(0:4,2,(ne-ns+1))
    real(8) :: ss(0:4), qtmat(4), reigv, q(0:4,2,ns:ne)   ! 2013_07_24 . scb
    
    reigv=1./eigv
    q=0._8
    do i=ns,ne
      ! scattering source update
      ss(0:4)=0._8
      do ms=1,ng
        ss(0:4)=ss(0:4)+xssn(i,m)%from(ms)*phin(1,0:4,i,ms)
      end do
      q(0:4,1,i)=ss(0:4)+reigv*xschin(i,m)*psin(0:4,i)
      ! transverse leakage
! 2013_08_12 . scb      
      !q(0:2,1,i)=q(0:2,1,i)-tlkg(1,0:2,i,m)
      !q(0:2,2,i)=-(2./5.)*tlkg(1,0:2,i,m)-(3./5.)*tlkg(2,0:2,i,m)
      !q(0:2,1,i)=q(0:2,1,i)-tlkg(1,0:2,i,m)+tspc(1,0:2,i,m)
      !q(0:2,2,i)=-(2./5.)*tlkg(1,0:2,i,m)-(3./5.)*tlkg(2,0:2,i,m)+tspc(2,0:2,i,m)
! 2014_10_06 . scb      
      q(0:2,1,i)=q(0:2,1,i)-tlkg(1,0:2,i,m)+tspc(1,0:2,i,m)+3.*tdmd(1,0:2,i,m)
      q(0:2,2,i)=-(2./5.)*tlkg(1,0:2,i,m)-(3./5.)*tlkg(2,0:2,i,m)+tspc(2,0:2,i,m)+1.2*tdmd(1,0:2,i,m)+1.4*tdmd(2,0:2,i,m)
! added end      
! added end      
      qtmat=qt(:,i,m)   ! qt = S_inv * D_inv
      do iod=0,4
        mqh(iod,:,i,m)=submatvecop(qtmat,q(iod,:,i))   ! q_tilda
      end do
    end do
  end subroutine
  
  subroutine UpdPsol(m,ng,ns,ne,ksq,mqh,mpsol)
    implicit none
    integer :: m, ng, ns, ne
    real(8),pointer :: ksq(:,:,:),mqh(:,:,:,:),mpsol(:,:,:,:)
    
    integer :: i, j
    real(8) :: k2(2), k4(2), k6(2), c(0:4,2)
    
    do i=ns,ne
      do j=1,2
        k2(j)=ksq(j,i,m)      
        k4(j)=k2(j)*k2(j)
        k6(j)=k4(j)*k2(j)
      end do
      do j=1,2
        c(0,j)=(mqh(0,j,i,m)*k4(j)+3.*mqh(2,j,i,m)*k2(j)+5.*mqh(4,j,i,m)*(21.+2.*k2(j)))/k6(j)
        c(1,j)=(15.*mqh(3,j,i,m)+mqh(1,j,i,m)*k2(j))/k4(j)
        c(2,j)=(mqh(2,j,i,m)*k2(j)+35.*mqh(4,j,i,m))/k4(j)
        c(3,j)=mqh(3,j,i,m)/k2(j)
        c(4,j)=mqh(4,j,i,m)/k2(j)
      end do
      mpsol(:,:,i,m)=c
    end do
    
  end subroutine
  
  subroutine SetRhs(m,ng,ns,ne,bc,ksq,s,mpsol,rhs)
    use senm1d, only : xsdn, xsd2n, h
    use senm1d, only : xsadfn   ! 2015_08_05 . scb for SP3 SENM ADF
    implicit none
    integer :: m, ng, ns, ne,bc(2)
    real(8),pointer :: ksq(:,:,:),s(:,:,:),mpsol(:,:,:,:),rhs(:,:)

    integer :: i,irow
! 2013_07_24 . scb    
    !real(8),dimension(2,(ne-ns+1)) :: k
    !real(8),dimension(4,(ne-ns+1)) :: sm,km,mm
    real(8),dimension(2,ns:ne) :: k
    real(8),dimension(4,ns:ne) :: sm,km,mm    
    real(8),dimension(2,ns:ne) :: adf    ! 2015_08_05 . scb for SP3 SENM ADF
! added end    
    real(8) :: temp(2), psol(0:4,2,2)

    do i=ns,ne
      k(1,i)=sqrt(ksq(1,i,m));    k(2,i)=sqrt(ksq(2,i,m))
      km(:,i)=kmat(k(1,i),k(2,i))
      sm(:,i)=s(:,i,m)
      mm(:,i)=mmat(xsdn(i,m),xsd2n(i,m),h(i))
      adf(:,i)=xsadfn(:,i,m)    ! 2015_08_05 . scb for SP3 SENM ADF
    end do
  
    irow=1
    !temp=SetBoundaryRHS(sm(:,ns), mm(:,ns), mpsol(:,:,ns,m), bc(1), 1)
    temp=SetBoundaryRHS(sm(:,ns), mm(:,ns), adf(2,ns), mpsol(:,:,ns,m), bc(1), 1)    ! 2015_08_05 . scb for SP3 SENM ADF
    rhs(:,irow)=temp
    do i=ns,ne-1
      psol(:,:,1)=mpsol(:,:,i,m)
      psol(:,:,2)=mpsol(:,:,i+1,m)
      irow=irow+1
      !temp=SetInterfaceMmtRHS(sm(:,i:i+1), psol)    ! 2015_08_05 . scb for SP3 SENM ADF
      temp=SetInterfaceMmtRHS(sm(:,i:i+1), [adf(1,i),adf(2,i+1)],psol)    ! 2015_08_05 . scb for SP3 SENM ADF
      rhs(:,irow)=temp
      irow=irow+1
      temp=SetInterfaceCurrentRHS(mm(:,i:i+1),sm(:,i:i+1),psol)
      rhs(:,irow)=temp
    end do
    irow=irow+1
    !temp=SetBoundaryRHS(sm(:,ne),mm(:,ne),mpsol(:,:,ne,m),bc(2),2)
    temp=SetBoundaryRHS(sm(:,ne),mm(:,ne),adf(1,ne),mpsol(:,:,ne,m),bc(2),2)    ! 2015_08_05 . scb for SP3 SENM ADF
    rhs(:,irow)=temp
    
  end subroutine

  subroutine UpdFlx(m,ns,ne,ksq,s,mpsol,Coeff,phin,mA,mB)
    implicit none
    integer :: m,ns,ne
    real(8),pointer :: ksq(:,:,:),s(:,:,:),mpsol(:,:,:,:),Coeff(:,:,:),phin(:,:,:,:), &
                       mA(:,:,:),mB(:,:,:)
    
    integer :: i,iod
    real(8) :: k(2), sm(4), mflx4th(0:4,2), flx4th(0:4,2)
    do i=ns,ne
      k(1)=sqrt(ksq(1,i,m))
      k(2)=sqrt(ksq(2,i,m))
      sm(:)=s(:,i,m)
      mA(:,i,m)=Coeff(:,2*i-1,m)
      mB(:,i,m)=Coeff(:,2*i,m)
    ! 4th order flux
      mflx4th(:,1)=flx_4th_expansion(mA(1,i,m),mB(1,i,m),mpsol(:,1,i,m),k(1))
      mflx4th(:,2)=flx_4th_expansion(mA(2,i,m),mB(2,i,m),mpsol(:,2,i,m),k(2))
      do iod=0,4
        flx4th(iod,:)=submatvecop(sm,mflx4th(iod,:))
      end do
      !do iod=0,4
      do iod=0,4
        phin(:,iod,i,m)=flx4th(iod,:)
      end do
      continue
    end do
  end subroutine
  
  function flx_4th_expansion(A,B,C,k)
    implicit none
    real(8) :: A, B, C(0:4), k
    real(8) :: k2, k4, sinhk, coshk, flx_4th_expansion(0:4)
    
    k2=k*k
    k4=k2*k2
       
    sinhk=sinh(k)
    coshk=cosh(k)
    
    flx_4th_expansion=0._8
    
    flx_4th_expansion(0)=c(0)+sinhk*A/k
    flx_4th_expansion(1)=c(1)+3._8*(coshk-sinhk/k)*B/k
    flx_4th_expansion(2)=c(2)-5._8*(3._8*coshk/k2-(1+3._8/k2)*sinhk/k)*A
    flx_4th_expansion(3)=c(3)+7._8/k*((1._8+15._8/k2)*coshk-(6._8+15._8/k2)*sinhk/k)*B
    flx_4th_expansion(4)=c(4)-9._8*(5._8/k2*(2._8+21._8/k2)*coshk-(1._8+45._8/k2+105._8/k4)*sinhk/k)*A

  end function

  subroutine UpdCurrent(ng,ns,ne,bc,ksq,s,mA,mB,mpsol,jnet1d)
    use senm1d, only : km,mm,sm,xsdn,xsd2n,h
    implicit none
    integer :: ng,ns,ne,bc(2)
    real(8) :: ksq(:,:,:),s(:,:,:),mA(:,:,:),mB(:,:,:),mpsol(:,:,:,:),jnet1d(:,:,:,:)
    
    integer,parameter :: l=1, r=2
    integer :: i,m
    real(8) :: k(2), jnet(2)
    
    do i=ns,ne
      do m=1,ng
        k(:)=sqrt(ksq(:,i,m))
        km(:,i,m)=kmat(k(1),k(2))
        sm(:,i,m)=s(:,i,m)
        mm(:,i,m)=mmat(xsdn(i,m),xsd2n(i,m),h(i))
      end do
    end do
    
    ! left-boundary
    if (bc(1).eq.0) then
      jnet1d(:,l,1,:)=0._8
    else
      do m=1,ng
        jnet=cal_j(mA(:,1,m),mB(:,1,m),mpsol(:,:,1,m),km(:,1,m),sm(:,1,m),mm(:,1,m),-1._8) ! original (-) sign
        jnet1d(:,l,1,m)=-jnet ! original (-) sign
      end do
    end if
    ! inner node
    do i=ns,ne-1
      do m=1,ng
        jnet=-cal_j(mA(:,i,m),mB(:,i,m),mpsol(:,:,i,m),km(:,i,m),sm(:,i,m),mm(:,i,m),1._8)
        jnet1d(:,r,i,m)=-jnet
        jnet1d(:,l,i+1,m)=jnet
      end do
    end do
    ! right boundary
    if (bc(2).eq.0) then
      jnet1d(:,r,ne,:)=0._8
    else
      do m=1,ng
        jnet=-cal_j(mA(:,ne,m),mB(:,ne,m),mpsol(:,:,ne,m),km(:,ne,m),sm(:,ne,m),mm(:,ne,m),1._8)
        jnet1d(:,r,ne,m)=-jnet
      end do
    end if
  end subroutine
  
  function cal_j(A,B,C,KM,SM,MM,x)
    implicit none
    real(8) :: A(2), B(2), C(0:4,2), KM(4), SM(4), MM(4), k(2)
    real(8) :: x, cal_j(2), MJH(2), MJP(2), SHK(4), CHK(4), MSK(4)
    integer ::iod
    
    k(1)=KM(1); k(2)=KM(4)
    MSK=multi_matop(MM,SM,KM)
    MJP=0._8
    do iod=0,4
      MJP=MJP+DLP(iod,x)*c(iod,:)
    end do
    SHK=sinhmat(k(1),k(2),x)
    CHK=coshmat(k(1),k(2),x)
    MJH=submatvecop(SHK,A)+submatvecop(CHK,B)
    cal_j=-submatvecop(MSK,MJH)-submatvecop(submatop(MM,SM),MJP)
  end function
  
  subroutine ChkPsiShp(ng,ns,ne,psin,psidn,errsum)
    implicit none
    
    integer :: ng, ns, ne
    real(8),pointer :: psin(:,:),psidn(:,:)
    real(8) :: errsum
    
    integer :: i, j
    real(8) :: sum,delta
    
    errsum=0._8
    do i=ns,ne
      sum=0._8
      do j=0,4
        sum=sum+(1._8/(2*j+1))*abs(psin(j,i)-psidn(j,i))
      end do
      if (psin(0,i).eq.0) then;  delta=0._8
      else; delta=sqrt(sum)/psin(0,i)
      end if
      errsum=errsum+delta
    end do
    errsum=errsum/(ne-ns+1)
  end subroutine
  

end module