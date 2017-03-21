Module Mod_FixedSource
  Logical :: iffixed
  Integer :: nSrctyp
  Real :: resid
  Integer,Allocatable,Dimension(:) :: iSrctyp
  Integer,Allocatable,Dimension(:,:) :: ExtSrcMap
  Real,Allocatable,Dimension(:,:) :: Srcdenza,Srcdenz
End Module