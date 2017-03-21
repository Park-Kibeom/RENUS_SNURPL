module ppr
    use const

    logical                             :: usemss  = TRUE
    logical                             :: term15  = TRUE
    integer                             :: nmaxppr = 20

    integer                             :: ncorn

    real                                :: sqrt2    = 1.41421356,     &
                                           rsqrt2   = 0.70710678
    
    integer                             :: npin

    real,pointer,dimension(:,:,:)       :: kappa
    real,pointer,dimension(:,:,:,:)     :: qf2d,    &   !(0:12,nxy,nz,ng)
                                           qc2d,    &   !(0:12,nxy,nz,ng)
                                           pc2d,    &   !(0:12,nxy,nz,ng)
                                           hc2d,    &   !(1:8,nxy,nz,ng)
                                           jcornx,  &   !(1:4,nxy,nz,ng)
                                           jcorny       !(1:4,nxy,nz,ng)

!   coefficient of least square fitting method
    real,pointer,dimension(:,:,:) ::       clsqf01, &   !(nxy,nz,ng)
                                           clsqf02, &   
                                           clsqf11, &   
                                           clsqf12, &   
                                           clsqf21, &   
                                           clsqf22, &   
                                           clsqf31, &   
                                           clsqf32, &   
                                           clsqf41, &   
                                           clsqf42, &   
                                           clsqfx1y1, &   
                                           clsqf1221, &   
                                           clsqf1331, &   
                                           clsqfx2y2
                                           
!   coefficient of general solutions
    real,pointer,dimension(:,:,:) ::       cpc02,   &   !(nxy,nz,ng)
                                           cpc04,   &
                                           cpc022,  &
                                           cpc11,   &
                                           cpc12,   &
                                           cpc21,   &
                                           cpc22,   &
                                           chc6,    &
                                           chc13j,  &
                                           chc13p,  &
                                           chc57j,  &
                                           chc57p,  &
                                           chc8j,   &
                                           chc8p,   &
                                           chc24j,  &
                                           chc24a

!   coefficiets of corner partial current
    real,pointer,dimension(:,:,:,:) ::     cpjxh1,  & !(1:4,nxy,nz,ng)
                                           cpjxh2,  &
                                           cpjxh5,  &
                                           cpjxh6,  &
                                           cpjxh7,  &
                                           cpjxh8,  &
                                           cpjxp6,  &
                                           cpjxp7,  &
                                           cpjxp8,  &
                                           cpjxp9,  &
                                           cpjxp11, &
                                           cpjxp12, &
                                           cpjxp13, &
                                           cpjxp14, &
                                           cpjyh3,  &
                                           cpjyh4,  &
                                           cpjyh5,  &
                                           cpjyh6,  &
                                           cpjyh7,  &
                                           cpjyh8,  &
                                           cpjyp2,  &
                                           cpjyp3,  &
                                           cpjyp4,  &
                                           cpjyp9,  &
                                           cpjyp10, &
                                           cpjyp12, &
                                           cpjyp13, &
                                           cpjyp14

    real,pointer,dimension(:,:,:) ::       phicorn, &   !(ncorn,nz,ng)
                                           phicorn0,&   !(ncorn,nz,ng)
                                           powvalr, &   !(npin,npin,la)
                                           powval,  &   !(npin,npin,ng)
                                           avgjnetz     !(ng,nxy,nz)

    real,pointer,dimension(:,:) ::         pppeak,  &   !(nxya,0:nz)
                                           powval1      !(npin,npin)

    real,pointer,dimension(:) ::           powtot       !(ng)

    real,pointer,dimension(:,:,:,:,:) ::   trlzcff      !(9,ng,nxy,nz)
     
    integer,pointer,dimension(:) ::        nodec(:,:),  &   !(0:nx+2,0:ny+2)
                                           lcc(:,:),    &   !(8,ncorn)            !1-w,2-n,3-e,4-s,5-nw,6-ne,7-se,8-sw
                                           lcn(:,:)         !(4,ncorn)            !1-nw,2-ne,3-se,4-sw

    integer,pointer,dimension(:) ::        nxcs,    &   !(ny+1)
                                           nxce,    &   !(ny+1)
                                           lctox,   &   !(0:ncorn)
                                           lctoy,   &   !(0:ncorn)
                                           lcnw,    &   !(nxy), corners belonging to a node
                                           lcsw,    &   !(nxy)
                                           lcne,    &   !(nxy)
                                           lcse         !(nxy)

    interface
        subroutine mallocppr(npinl)
            integer             :: npinl
        end subroutine
        subroutine initppr(iftran)
            logical             :: iftran
        end subroutine
        subroutine caltrlz(iftran)
            logical             :: iftran
        end subroutine
        subroutine calhomo()
        end subroutine
    end interface
end module