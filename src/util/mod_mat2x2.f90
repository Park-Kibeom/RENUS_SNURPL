module mat2x2
    implicit none

    contains
        subroutine matxmat2x2(a,b,c)
            real          :: a(2,2),b(2,2),c(2,2)

            c(1,1)=a(1,1)*b(1,1)+a(2,1)*b(1,2)
            c(2,1)=a(1,1)*b(2,1)+a(2,1)*b(2,2)
            c(1,2)=a(1,2)*b(1,1)+a(2,2)*b(1,2)
            c(2,2)=a(1,2)*b(2,1)+a(2,2)*b(2,2)             
        end subroutine 

        subroutine matxvec2x2(a,b,x)
            real          :: a(2,2),b(2),x(2)

            x(1)=a(1,1)*b(1)+a(2,1)*b(2)
            x(2)=a(1,2)*b(1)+a(2,2)*b(2)
        end subroutine 

        subroutine addmat2x2(a,b,c)
            real          :: a(2,2),b(2,2),c(2,2)

            c(1,1)=a(1,1) + b(1,1)
            c(2,1)=a(2,1) + b(2,1)
            c(1,2)=a(1,2) + b(1,2)
            c(2,2)=a(2,2) + b(2,2)
            
        end subroutine     

        subroutine invmat2x2(a,ainv)
            real          :: a(2,2),ainv(2,2)
            real          :: rdet,tmp

            rdet = 1/(a(1,1)*a(2,2)-a(2,1)*a(1,2))
            tmp  = a(1,1)
            ainv(1,1)=rdet*a(2,2)
            ainv(2,1)=-rdet*a(2,1)
            ainv(1,2)=-rdet*a(1,2)
            ainv(2,2)=rdet*tmp
        end subroutine 

        subroutine solmat2x2(a,b,x)
            real          :: a(2,2),b(2),x(2)
            real          :: ainv(2,2)

            call invmat2x2(a,ainv)
            call matxvec2x2(ainv,b,x)
        end subroutine 

end module