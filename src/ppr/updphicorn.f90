subroutine updphicorn(k)
    use ppr
    use geom,   only : ng
    integer,intent(in) :: k

    real :: ref(ng)
    real ,pointer,save,dimension(:,:,:) :: phicornl

    if(.not.associated(phicornl)) then
        allocate(phicornl(ncorn,nz,ng))
    endif     

    do m=1,ng
    do lc=1,ncorn

        inw=lcn(1,lc)
        ine=lcn(2,lc)
        isw=lcn(4,lc)
        ise=lcn(3,lc)

!       distinguish the existence nodes around each corner point
        itype=0
        if(inw/=0) itype=itype+8
        if(ine/=0) itype=itype+4
        if(isw/=0) itype=itype+2
        if(ise/=0) itype=itype+1

        select case(itype)
        case(1) !nw=0 ne=0 sw=0 se=1
        phicorn(lc,k,m)=phicorn(lc,k,m)                             &
                        -(jcornx(1,ise,k,m)+jcorny(1,ise,k,m))
        case(2) !nw=0 ne=0 sw=1 se=0
        phicorn(lc,k,m)=phicorn(lc,k,m)                             &
                        +(jcornx(2,isw,k,m)-jcorny(2,isw,k,m))
        case(3) !nw=0 ne=0 sw=1 se=1
        phicorn(lc,k,m)=phicorn(lc,k,m)                             &
                        +.5*(jcornx(2,isw,k,m)-jcornx(1,ise,k,m)    &
                        -jcorny(2,isw,k,m)-jcorny(1,ise,k,m))
        case(4) !nw=0 ne=1 sw=0 se=0
        phicorn(lc,k,m)=phicorn(lc,k,m)                             &
                        +(-jcornx(4,ine,k,m)+jcorny(4,ine,k,m))
        case(5) !nw=0 ne=1 sw=0 se=1
        phicorn(lc,k,m)=phicorn(lc,k,m)                             &
                        +.5*(-jcornx(4,ine,k,m)-jcornx(1,ise,k,m)   &
                        +jcorny(4,ine,k,m)-jcorny(1,ise,k,m))
        case(8) !nw=1 ne=0 sw=0 se=0
        phicorn(lc,k,m)=phicorn(lc,k,m)                             &
                        +(jcornx(3,inw,k,m)+jcorny(3,inw,k,m))
        case(10) !nw=1 ne=0 sw=1 se=0      
        phicorn(lc,k,m)=phicorn(lc,k,m)                             &
                        +.5*(jcornx(3,inw,k,m)+jcornx(2,isw,k,m)    &
                        +jcorny(3,inw,k,m)-jcorny(2,isw,k,m))
        case(12) !nw=1 ne=1 sw=0 se=0
        phicorn(lc,k,m)=phicorn(lc,k,m)                             &
                        +.5*(jcornx(3,inw,k,m)-jcornx(4,ine,k,m)    &
                        +jcorny(3,inw,k,m)+jcorny(4,ine,k,m))
        case(15) !nw=1 ne=1 sw=1 se=1
        phicorn(lc,k,m)=phicorn(lc,k,m)                             &
                        +.25*(jcornx(3,inw,k,m)-jcornx(4,ine,k,m)  &
                            +jcornx(2,isw,k,m)-jcornx(1,ise,k,m)    &
                            +jcorny(3,inw,k,m)-jcorny(2,isw,k,m)    &
                            +jcorny(4,ine,k,m)-jcorny(1,ise,k,m))
        end select

!        if(itype .ne. 15) phicorn(lc,k,m)=0 

        if(phicorn(lc,k,m).lt.zero) then
            phicorn(lc,k,m)=phicorn0(lc,k,m)
        endif
    enddo ! of lc
    enddo ! of group

    return
end subroutine
