module senmop
  real(8), private, pointer, dimension(:,:) :: temp_sol,temp_ya
  integer, private, pointer, dimension(:) :: temp_RowElmtIdx
  
  contains

  subroutine initsenmop(dim1)
    integer :: dim1
    allocate(temp_sol(2,dim1))
    allocate(temp_ya(2,dim1))
    allocate(temp_RowElmtIdx(dim1))
    temp_sol=0._8
    temp_ya=0._8
    temp_RowElmtIdx=0
  endsubroutine
  
  subroutine SetSimTrans(ng,qt,ksq,s)
    use matop_sp3, only : matinv,submatop,multi_matop,seteigs
    use senm1d, only : hinv,ns,ne,xsrn,xstrn,xsdn,xsd2n
    implicit none
    integer :: ng
    real(8),pointer,dimension(:,:,:) :: qt,ksq,s
    
    integer :: i, m
    real(8) :: d0, d2, sigr, sigt, hsqinv4
    real(8) :: A(4), D(4), DINV(4), AHAT(4), SINV(4), G(4), LAMBDA(2), EIGVEC(4)
    
    do i=ns,ne
      hsqinv4=4._8*(hinv(i)**2)   ! 4/h^2
      do m=1,ng
        d0=xsdn(i,m)*hsqinv4      ! 4*D0/h^2
        d2=xsd2n(i,m)*hsqinv4     ! 4*D2/h^2
        sigr=xsrn(i,m)            ! XSr
        sigt=xstrn(i,m)           ! XSt (or XStr)
        ! set up D matrix
        D=0;  A=0
        D(1)=d0;        D(2)=2.*d0
        D(3)=2./5.*d0;  D(4)=(4.*d0+3.*d2)/5.
        ! D = 4/h^2 * [ D0   2*D0   |    2/5*D0    4/5*D0 + 3/5*D2 ]
        A(1)=sigr;      A(4)=sigt
        DINV=matinv(D)
        AHAT=submatop(DINV,A)   ! A_hat = D_inv * A
        call seteigs(AHAT,LAMBDA,EIGVEC)   ! k2(1:2), S
        SINV=matinv(EIGVEC)
        G=multi_matop(SINV,AHAT,EIGVEC)   ! S_inv * A_hat * S
        qt(:,i,m)=submatop(SINV,DINV)    ! qt = S_inv * D_inv
        ksq(:,i,m)=LAMBDA
        s(:,i,m)=EIGVEC
        continue
      end do
    end do
  end subroutine
  
  subroutine Senmmat_LU(mg,N)
    use matop_sp3, only : matinv,submatop
    use senm1d, only : nElmt,ElmtLoc,DiagIdx,Elmt,LowerElmtIdx,nLowerElmt,LU
    implicit none
    integer,intent(in) :: mg, N
    integer :: i, j, k
    integer :: m, irow, icol, irow0, icol0, idx
    integer,pointer :: RowElmtIdx(:)
    real(8) :: Mult(4), InvDiag(4), temp(4)
    
    RowElmtIdx=>temp_RowElmtIdx
    
    do i=1,N
      LU(:,:,i,mg)=Elmt(:,:,i,mg)
    end do
    
    do i=1,N-1
      irow=i
      m=nLowerElmt(irow,mg); idx=DiagIdx(irow,mg)
      InvDiag=matinv(LU(:,idx,irow,mg))
      do j=1,m
        RowElmtIdx=0
        irow0=LowerElmtIdx(j,irow,mg)
        do k=1,nElmt(irow0,mg) 
          icol0=ElmtLoc(k,irow0,mg)
          RowElmtIdx(icol0)=k
        end do
        idx=RowElmtIdx(irow)
        Mult=submatop(LU(:,idx,irow0,mg), InvDiag)
        LU(:,idx,irow0,mg)=Mult
        continue
        do k=DiagIdx(irow,mg)+1,nElmt(irow,mg)
          icol=ElmtLoc(k,irow,mg)
          if (RowElmtIdx(icol).eq.0) cycle
          temp=submatop(mult,LU(:,k,irow,mg))
          Idx=RowElmtIdx(icol)
          LU(:,idx,irow0,mg)=LU(:,idx,irow0,mg)-temp
        end do
      end do
    end do
    
    do i=1,N
      irow=i; idx=DiagIdx(irow,mg)
      InvDiag=matinv(LU(:,idx,irow,mg))
      LU(:,idx,irow,mg)=InvDiag
    end do
    
  end subroutine

  subroutine SolSenmMatLU(ig,N,rhs,Coeff)
    use matop_sp3, only : submatvecop
    use senm1d, only : DiagIdx,ElmtLoc,nElmt,LU,Elmt
    implicit none
    
    integer :: ig,N
    real(8),pointer :: rhs(:,:), Coeff(:,:,:)
    
    integer :: i,j,k,l,m
    real(8) :: lmnt(2)
    real(8),pointer,dimension(:,:) :: sol, ya

    sol=>temp_sol
    ya=>temp_ya
    ! forward
    do i=1,N
      ya(:,i)=0._8
    end do
    ya(:,1)=rhs(:,1)
    do i=2,N
      ya(:,i)=rhs(:,i)
      k=DiagIdx(i,ig)
      do j=1,k-1
        l=ElmtLoc(j,i,ig)
        ya(:,i)=ya(:,i)-submatvecop(LU(:,j,i,ig),ya(:,l))
      end do
    end do
    ! backward
    k=DiagIdx(N,ig)
    sol(:,N)=submatvecop(LU(:,k,N,ig),ya(:,N))
    Coeff(:,N,ig)=submatvecop(LU(:,k,N,ig),ya(:,N))
    do i=N-1,1,-1
      lmnt=ya(:,i)
      k=DiagIdx(i,ig)
      m=nElmt(i,ig)
      do j=k+1,m
        l=ElmtLoc(j,i,ig)
        lmnt=lmnt-submatvecop(LU(:,j,i,ig),sol(:,l))
      end do
      sol(:,i)=submatvecop(LU(:,k,i,ig),lmnt)
      Coeff(:,i,ig)=submatvecop(LU(:,k,i,ig),lmnt)
    end do
    ! check solution
    do i=1,N
      m=nElmt(i,ig)
      lmnt=0
      do j=1,m
        l=ElmtLoc(j,i,ig)
        lmnt=lmnt+submatvecop(Elmt(:,j,i,ig),sol(:,l))
      end do
      ya(:,i)=lmnt
    end do
    do i=1,N
      ya(:,i)=ya(:,i)-rhs(:,i)
    end do
    
  end subroutine
  
  function kmat(k1,k2)
    implicit none
    real(8) :: k1,k2,kmat(4)
    kmat=0
    kmat(1)=k1; kmat(4)=k2
  end function

  function sinhmat(k1,k2,x)
    implicit none
    real(8) :: k1, k2, x, sinhmat(4)
    sinhmat=0
    sinhmat(1)=sinh(k1*x)
    sinhmat(4)=sinh(k2*x)
  end function

  function coshmat(k1,k2,x)
    implicit none
    real(8) :: k1, k2, x, coshmat(4)
    coshmat=0
    coshmat(1)=cosh(k1*x)
    coshmat(4)=cosh(k2*x)
  end function

  function mmat(d0,d2,hl)
    implicit none
    real(8) :: d0, d2, hl,hlinv,mmat(4)
    mmat=0
    hlinv=2._8/hl
    mmat(1)=hlinv*d0
    mmat(2)=2._8*hlinv*d0
    mmat(4)=hlinv*d2
  end function
  
  function DLP(I,X)
    implicit none
    integer :: i,j
    real(8) :: x,y,DLP,LPCoeff(0:4,0:4)
    data LPCoeff /  0._8,    0._8,    0._8,    0._8,    0._8, &
                    1._8,    0._8,    0._8,    0._8,    0._8, &
                    0._8,    3._8,    0._8,    0._8,    0._8, &
                  -1.5_8,    0._8,   7.5_8,    0._8,    0._8, &
                    0._8,  -7.5_8,    0._8,  17.5_8,    0._8 / 
    y=0
    do j=0,i
      y=y+(x**j)*LPCoeff(j,i)
    end do
    DLP=y
  end function
  
  end module



  !  
