module matop_sp3
  use const
  use geom, only : nxy,nz
  implicit none
  interface scalprod
    module procedure scalprod1g
    module procedure scalprodmg
  end interface
  
  interface multi_matop
    module procedure multi_matop2
    module procedure multi_matop3
    module procedure multi_matop4
    module procedure multi_matop5
  end interface
  
contains

function vecprod(A,B,n)
  real(8) :: A(:,:), B(:,:)
  integer :: l, n
  real(8) :: vecprod
  vecprod=0._8
  do l=1,n
    vecprod=vecprod+A(l,1)*B(l,1)
  end do
end function

function matinv(A)
  implicit none
  real(8) :: A(4), matinv(4), det
  det=1._8/(A(1)*A(4)-A(2)*A(3))
  matinv(1)=A(4)
  matinv(2)=-A(2)
  matinv(3)=-A(3)
  matinv(4)=A(1)
  matinv=det*matinv
end function

function submatop(A,B)
  implicit none
  real(8) :: A(4), B(4), submatop(4)
  submatop(1)=A(1)*B(1)+A(2)*B(3)
  submatop(3)=A(3)*B(1)+A(4)*B(3)
  submatop(2)=A(1)*B(2)+A(2)*B(4)
  submatop(4)=A(3)*B(2)+A(4)*B(4)
end function

function submatvecop(A,B)
  implicit none
  real(8) ::A(4), B(2), submatvecop(2)
  submatvecop(1)=A(1)*B(1)+A(2)*B(2)
  submatvecop(2)=A(3)*B(1)+A(4)*B(2)
end function

subroutine set_eigc(MAT,LAMBDA,EIGVEC)
  implicit none
  real(8) :: MAT(4), LAMBDA(2), EIGVEC(4)
  real(8) :: MAT2(4), a, b, c, d, temp
  equivalence(MAT2(1), a)
  equivalence(MAT2(2), b)
  equivalence(MAT2(3), c)
  equivalence(MAT2(4), d)
  
  MAT2=MAT
  temp=sqrt(a*a+4._8*b*c-2._8*a*d+d*d)
  LAMBDA(1)=0.5_8*(a+d-temp)
  LAMBDA(2)=0.5_8*(a+d+temp)
  EIGVEC(1)=-0.5_8*(-a+d+temp)/c
  EIGVEC(2)=0.5_8*(a-d+temp)/c
  EIGVEC(3)=1
  EIGVEC(4)=1
  
end subroutine

function multi_matop2(A1,A2)
  implicit none
  real(8) :: A1(4), A2(4)
  real(8) :: multi_matop2(4)
  multi_matop2=submatop(A1,A2)
end function

function multi_matop3(A1,A2,A3)
  implicit none
  real(8) :: A1(4), A2(4), A3(4)
  real(8) :: multi_matop3(4), temp(4)
  temp=submatop(A2,A3)
  multi_matop3=submatop(A1,temp)
end function

function multi_matop4(A1,A2,A3,A4)
  implicit none
  real(8) :: A1(4), A2(4), A3(4), A4(4)
  real(8) :: multi_matop4(4), temp1(4), temp2(4)
  temp1=submatop(A3,A4)
  temp2=submatop(A2,temp1)
  multi_matop4=submatop(A1,temp2)
end function

function multi_matop5(A1,A2,A3,A4,A5)
  implicit none
  real(8) :: A1(4), A2(4), A3(4), A4(4), A5(4)
  real(8) :: multi_matop5(4), temp1(4), temp2(4)
  temp1=submatop(A4,A5)
  temp2=submatop(A3,temp1)
  temp1=submatop(A2,temp2)
  multi_matop5=submatop(A1,temp1)
end function

subroutine seteigs(mat,lambda,eigvec)
  implicit none
  real(8) :: mat(4), lambda(2), eigvec(4)
  real(8) :: mat2(4), a, b, c, d
  real(8) :: temp
  equivalence(mat2(1),a)
  equivalence(mat2(2),b)
  equivalence(mat2(3),c)
  equivalence(mat2(4),d)
  mat2=mat
  eigvec=1._8
  temp=sqrt(a*a+4._8*b*c-2._8*a*d+d*d)
  lambda(1)=0.5_8*(a+d-temp)
  lambda(2)=0.5_8*(a+d+temp)
  eigvec(1)=-0.5_8*(-a+d+temp)/c
  eigvec(2)=0.5_8*(a-d+temp)/c
end subroutine

      function scalprod1g(x,y)
    real(8),pointer :: x(:,:),y(:,:)
    integer :: k, l, m
    real(8) :: scalprod1g
        
    scalprod1g=0
    do k=1,nz
      do l=1,nxy
        scalprod1g=scalprod1g+x(l,k)*y(l,k)
      end do
    end do

    return
  end function

  function scalprodmg(ng,x,y)
    real(8),pointer :: x(:,:,:),y(:,:,:)
    integer :: ng
    integer :: k, l, m
    real(8) :: scalprodmg
        
    scalprodmg=0
    do k=1,nz
      do l=1,nxy
        do m=1,ng
          scalprodmg=scalprodmg+x(m,l,k)*y(m,l,k)
        end do
      end do
    end do

    return
  end function    
end module