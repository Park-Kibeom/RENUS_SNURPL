    subroutine initbd
!
! initialize common blocks that are not specified by block data
      use param
!
      include 'global.h'
      include 'cards.h'
      include 'thfuel.inc'
!
      nccard=0
      ngcard=0
      nxcard=0
      npcard=0
      nfcard=0
      ntrcard=0
      nthcard=0
      n1dcard=0
      nplcard=0
! 
      ccard=BLANK
      gcard=BLANK
      xcard=BLANK
      pcard=BLANK
      fcard=BLANK
      trcard=BLANK
      thcard=BLANK
      onedcard=BLANK
      plcard=BLANK
!
! xsec cards
      i=0
      i=i+1
      xcard(i)='GROUP_SPEC'
      i=i+1
      xcard(i)='COMP_NUM'
      nxcard=i
! 
! geom cards
      i=0
      i=i+1
      gcard(1)='RAD_CONF'
      i=i+1
      gcard(2)='GRID_Z'
      i=i+1
      gcard(3)='NEUTMESH_Z'
      i=i+1 
      gcard(4)='ASSY_TYPE'
      ngcard=i
!
!thfuel.inc
      akfuel(0)   = 1.05
      akfuel(1)   = 0.
      akfuel(2)   = 0.
      akfuel(3)   = 0.
      akfuel(4)   = 2150.
      akfuel(5)   = 73.15
      akclad(0)   = 7.51
      akclad(1)   = 2.09e-2
      akclad(2)   = -1.45e-5
      akclad(3)   = 7.67e-9

      return
    end
