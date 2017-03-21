module misc

    include 'global.h'
    include 'files.h'
    include 'xesm.h'
    include 'geom.h'
    include 'pow.h'
    include 'thgeom.inc'
    include 'ffdm.h'
    include 'xsec.h'
    include 'thop.inc'

    contains

    function misc_fluxlevel(flux) result(fnorm)
    ! compute absolute level for given flux

        integer :: k, l, m
        integer :: la, iat, vol

        real(8) :: fnorm
        real,pointer,dimension(:,:,:) :: flux !(1:ng,0:nxy,0:nz)
        real(8) :: totpow, powfac

        totpow=0.d0
        do k=kfbeg,kfend
          do l=1,nxy
            la=ltola(l)
            iat=iassytyp(la)
            if(.not.iffuela(iat)) cycle

            vol=volnode(l,k)
            do m=1,ng
              totpow=totpow+flux(m,l,k)*xskpf(m,l,k)*vol
            enddo
          enddo
        enddo
      
        powfac=totpow/volfuel*hac*pfa*pfa*1.e6 !cubic m to cubic cm
        if(hex) powfac=powfac*0.86602540378444 !txk: hex assy area correction
        fnorm=(plevel00*powfa)/powfac

        return 
    end function


    subroutine misc_rate(flux, mxs, rate)
    ! calculates rate for given flux and macro cross sections
    ! rate(:,:) = ( mxs(:,:,:) * flux(:,:,:) )

        integer :: k, l, m, la, iat

        real,pointer,dimension(:,:,:) :: flux
        real,pointer,dimension(:,:,:) :: mxs
        real,pointer,dimension(:,:)   :: rate

        rate(:,:) = 0.d0

        do k=kfbeg,kfend
            do l=1,nxy
                la=ltola(l)
                iat=iassytyp(la)
                if(.not.iffuela(iat)) cycle

                do m=1,ng
                    rate(l,k)=rate(l,k)+mxs(m,l,k)*flux(m,l,k)
                enddo

            enddo
        enddo  
    end subroutine misc_rate


end module misc