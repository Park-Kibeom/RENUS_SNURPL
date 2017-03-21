Subroutine ExternalSource_Map

use Mod_FixedSource

include 'global.h'
include 'geom.h'

Integer :: ig,iy,ix,iz,iza,l,la,ja,ia,iSrctyp1,i


Do iy=1,ny
  Do ix=nxs(iy),nxe(iy)
    l=nodel(ix,iy)
    la=ltola(l)
    ja=latoja(la)
    ia=latoia(la)
    iSrctyp1=iSrctyp(la)
    ExtSrcMap(ix,iy) = iSrctyp(la)
  EndDo
EndDo

Do iSrctyp1=1,nSrctyp
  i=1
  Do iza=1,nza
    Do iz=1,nmeshz(iza)
      Srcdenz(i,iSrctyp1)=Srcdenza(iza,iSrctyp1)
      i=i+1
    EndDo
  EndDo
EndDo

return
End Subroutine