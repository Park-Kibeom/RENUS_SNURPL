! cards.h - input cards,speical characters and character variables for i/o
      character*1 DOT,BANG,BLANK,SLASH,AST
      parameter (DOT='.',BANG='!',BLANK=' ',SLASH='/',AST='*')
      parameter (mxncol=512,mxnfield=128)
      parameter (mxnccard=50,mxngcard=50,mxnxcard=50,mxnpcard=50)
      parameter (mxnfcard=50,mxntrcard=50,mxnthcard=50,mxn1dcard=50)
      parameter (mxnplcard=50)
      character*1 probe,sc(mxncol),sc80(80)
      character*512  oneline
      character*10 cardname,blockname
	  character*7 form7  
      character*79 hbar
      logical iffile
      equivalence (probe,oneline)
      equivalence (sc(1),oneline)
      character*10 ccard(mxnccard),gcard(mxngcard),xcard(mxnxcard),  &
                   pcard(mxnpcard),fcard(mxnfcard),trcard(mxntrcard),  &
                   thcard(mxnthcard),onedcard(mxn1dcard),  &
                   plcard(mxnplcard)                         
      logical      ifccard(mxnccard),ifgcard(mxngcard),ifxcard(mxnxcard),  &
                   ifpcard(mxnpcard),iffcard(mxnfcard),  &
                   iftrcard(mxntrcard),ifthcard(mxnthcard),  &
                   if1dcard(mxn1dcard),ifplcard(mxnplcard)                                           !onedk
      common /ncards/nccard,ngcard,nxcard,npcard,nfcard,ntrcard,nthcard,  &
                     n1dcard,nplcard                                                     !onedk
      common /acards/ccard,gcard,xcard,pcard,fcard,trcard,thcard,  &
                     onedcard,plcard                                                    !onedk
      common /lcards/ifccard,ifgcard,ifxcard,ifpcard,iffcard,iftrcard,  &
                     ifthcard,if1dcard,ifplcard                                           !onedk
      common /cconst/hbar
      logical ldum
      common /dumfield/idum(mxnfield),fdum(mxnfield),  &
                       ldum(mxnfield),nrepeat(mxnfield) !20 fields at maximum per line
