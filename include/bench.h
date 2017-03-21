      integer, parameter   :: IDM=1, IPPM=2,ITF=3,ITM=4, NPROP=4
      integer :: ntotpnt
      integer, pointer     :: npnt(:,:)   !NPROP, ncomp

      real   , dimension(:,:,:), pointer :: pnts(:,:,:),  & !ipnt,NPROP, ncomp
                                            sigbchtr,  &
                                            sigbcha,  &
                                            sigbchnf,  &
                                            sigbchkf,  &
                                            sigbchadf,  &
                                            sigbchs(:,:,:,:)

      common /benchi1/ntotpnt
      common /benchi/ npnt
      common /benchr/ pnts,sigbchtr,sigbcha,sigbchnf,sigbchkf,sigbchadf
      common /benchr4/sigbchs