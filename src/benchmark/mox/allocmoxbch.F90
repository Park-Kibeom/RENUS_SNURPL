    subroutine allocmoxbch
      use param
      use allocs

      include 'global.h'
      include 'geom.h'
      include 'xsec.h'
      include 'moxbch.inc'
      
      call dmalloc(npnt, NPROP, ncomp)
      call dmalloc(pnts,3,NPROP,ncomp)
      ntotpnt=27

      call dmalloc(sigbchtr, ntotpnt, ng, ncomp)
      call dmalloc(sigbcha, ntotpnt, ng, ncomp)
      call dmalloc(sigbchnf, ntotpnt, ng, ncomp)
      call dmalloc(sigbchkf, ntotpnt, ng, ncomp)
      call dmalloc(sigbchadf, ntotpnt, ng, ncomp)
      call dmalloc(sigbchs, ntotpnt, ng, ng, ncomp)
      call dmalloc(sigchi, ng, ncomp)

    end subroutine