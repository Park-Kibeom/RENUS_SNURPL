Subroutine Edit_Depl(id,burnn,eigv)

use param

include 'global.h'
include 'files.h'
include 'cntl.h'
include 'thlink.inc'

Logical :: First=.true.
Integer :: id
Real :: eigv,burnn

If(First) then
  Write(5689,*) id, 0., eigv
  First=.false.
Else
  Write(5689,*) id, burnn, eigv
Endif

End Subroutine