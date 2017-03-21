module matop
    implicit none
    
    contains

    !add outmat=mat1+mat2
    subroutine addmat(mat1,mat2,outmat,n)
        integer      :: n
        integer                 :: i,j
        real :: mat1(n,n), mat2(n,n), outmat(n,n)

        do i=1,n
        do j=1,n
            outmat(j,i)=mat1(j,i)+mat2(j,i)
        enddo
        enddo

        return
    end subroutine

    !subtract outmat=mat1-mat2
    subroutine subtractmat(mat1,mat2,outmat,n)
        integer      :: n
        integer                 :: i,j
        real :: mat1(n,n), mat2(n,n), outmat(n,n)

        do i=1,n
        do j=1,n
            outmat(j,i)=mat1(j,i)-mat2(j,i)
        enddo
        enddo

        return
    end subroutine


    !multiply outmat=mat1 x mat2
    subroutine matxmat(mat1, mat2, outmat, n)
        integer      :: n
        integer                 :: i,j,k
        real :: mat1(n,n), mat2(n,n), outmat(n,n)

        do i=1,n
        do j=1,n
            outmat(j,i)=0
            do k=1,n
                outmat(j,i)=outmat(j,i)+mat1(k,i)*mat2(j,k)
            enddo
        enddo
        enddo

        return
    end subroutine

    !multiply outmat=mat x diagmat
    subroutine matxdiag(mat, diagmat, outmat, n)
        integer      :: n
        integer                 :: i,j
        real :: mat(n,n), diagmat(n), outmat(n,n)

        do i=1,n
        do j=1,n
            outmat(j,i)=mat(j,i)*diagmat(j)
        enddo
        enddo

        return
    end subroutine

    !multiply outmat=diagmat x mat
    subroutine diagxmat(diagmat, mat, outmat, n)
        integer      :: n
        integer                 :: i,j
        real :: mat(n,n), diagmat(n), outmat(n,n)

        do i=1,n
        do j=1,n
            outmat(j,i)=diagmat(i)*mat(j,i)
        enddo
        enddo

        return
    end subroutine

    !multiply outvec=mat x vec
    subroutine matxvec(mat, vec, outvec, n)
        integer      :: n
        integer                 :: i,j
        real :: mat(n,n), vec(n), outvec(n)

        do i=1,n
            outvec(i)=0
            do j=1,n
            outvec(i)=outvec(i)+mat(j,i)*vec(j)
            enddo
        enddo

        return
    end subroutine 

    ! make l,u
    subroutine makelu(lu, pivot, n)
        integer      :: n
        real :: lu(n,n)
        integer                 :: pivot(n)
        integer                 :: i,j,k
        

        do i=1,n
            pivot(i)=i
        enddo

        do k=1,n
            ! if 0, do pivot
            if(lu(k,k) .eq. 0 ) then
                do i=k+1,n
                    if(lu(k,i) .ne. 0) then
                        pivot(k)=i
                        pivot(i)=k
                        exit
                    endif
                enddo
            endif

            lu(k,pivot(k))=1/lu(k,pivot(k))
            do i=k+1,n
                lu(k,pivot(i))=lu(k,pivot(i))*lu(k,pivot(k))
                do j=k+1,n
                    lu(j,pivot(i))=lu(j,pivot(i))-lu(k,pivot(i))*lu(j,pivot(k))
                enddo
            enddo
        enddo

        return
    end subroutine


    ! outvec= mat\vec.
    subroutine invmatxvec(matlu,vec,pivot,outvec,n)
        integer      :: n
        real, dimension(:,:)    :: matlu(n,n), vec(n), outvec(n)
        real, dimension(:)      :: y(n)
        integer                 :: pivot(n)
        integer                 :: i,j,k        

        ! forward
        do i=1,n
            y(pivot(i))=vec(pivot(i));
            do j=1,i-1
                y(pivot(i))=y(pivot(i))-matlu(j,pivot(i))*y(pivot(i))
            enddo
        enddo

        ! backward
        do i=n,1,-1
            outvec(pivot(i))=y(pivot(i))
            do j=n,i-1
                outvec(pivot(i))=outvec(pivot(i))-matlu(j,pivot(i))*outvec(pivot(i))
            enddo
            outvec(pivot(i))=outvec(pivot(i))*matlu(i,pivot(i))
        enddo

    end subroutine


    subroutine invmatxvec1(mat,vec,sol,n)
        integer      :: n
        real, dimension(:,:)    :: mat(n,n),vec(n),outvec(n),sol(n)
        integer                 :: i,j,k,it
        integer                 :: pivot(n),ipmax
        real                    :: fm, pmax,sum
        logical                 :: nonzero

        ! initialize pivot
        do i=1,n
            pivot(i)=i
        enddo

        ! forward elimination
        do k=1,n-1
            pmax=abs(mat(k,pivot(k)))
            ipmax=k
            do i=k+1,n
                if(abs(mat(k,pivot(i))) .gt. pmax) then
                    pmax=abs(mat(k,pivot(i)))
                    ipmax=i
                endif
            enddo

            it=pivot(ipmax)
            pivot(ipmax)=pivot(k)
            pivot(k)=it

            do i=k+1,n
                fm=mat(k,pivot(i))/mat(k,pivot(k))
                if(fm.eq.0) cycle
                do j=k,n
                    mat(j,pivot(i))=mat(j,pivot(i))-fm*mat(j,pivot(k))
                enddo
                vec(pivot(i))=vec(pivot(i))-fm*vec(pivot(k))
            enddo

        enddo


        !backward elimination
        outvec(pivot(n))=vec(pivot(n))/mat(n,pivot(n))
        do i=n-1,1,-1
            sum=0
            do j=i+1,n
                sum=sum+mat(j,pivot(i))*outvec(pivot(j))
            enddo
            outvec(pivot(i))=(vec(pivot(i))-sum)/mat(i,pivot(i))
        enddo

        do i=1,n
            sol(i)=outvec(pivot(i))
        enddo

        return
    end subroutine

! 2012_09_28 . scb    
    subroutine invmat(mat, inv, n)
        integer                 :: n
        real                    :: mat(n,n), inv(n,n) 
        
        integer                 :: i, j, k
        real                    :: fm
        real                    :: augmat(2*n,n)
            
!       augment input matrix with an identity matrix
        augmat(:,:) = 0
        do i=1,n
            do j=1,n
                augmat(j,i) = mat(j,i)
            enddo
            augmat(n+i,i) = 1
        end do
                
!       forward elimination
        do k=1,n-1
            augmat(k,k)=1/augmat(k,k)
            do i=k+1,n
                fm=augmat(k,i)*augmat(k,k)
                if(fm.eq.0) cycle
                augmat(k,i)=0
                do j=k+1,2*n
                    augmat(j,i)=augmat(j,i)-fm*augmat(j,k)
                enddo
            enddo
        enddo
        augmat(n,n)=1/augmat(n,n)
        
!       backward elimination
        do k=n,2,-1
            do i=k-1,1,-1
                fm=augmat(k,i)*augmat(k,k)
                augmat(k,i)=0
                do j=k+1,2*n
                    augmat(j,i)=augmat(j,i)-fm*augmat(j,k)
                enddo
            enddo
        enddo

!       make diagonal to 1
        do i=1,n
            do j=n+1,2*n                         
               augmat(j,i) = augmat(j,i)*augmat(i,i)
            enddo
            augmat(i,i)=1.
        enddo
            
        do i=1,n
            do j=1,n
                inv(j,i) = augmat(j+n,i)
            enddo
        enddo
    end subroutine 
! added end    
end module
