SUBROUTINE NUTYP

USE MASTERXSL
USE DEPLMOD

nu2id(1)  = 'U234'
nu2id(2)  = 'U235'
nu2id(3)  = 'U236'
nu2id(4)  = 'NP37'
nu2id(5)  = 'U238'
nu2id(6)  = 'PU48'
nu2id(7)  = 'NP39'
nu2id(8)  = 'PU49'
nu2id(9)  = 'PU40'
nu2id(10) = 'PU41'
nu2id(11) = 'PU42'
nu2id(12) = 'AM43'
nu2id(13) = 'TRU '
!nu2id(13) = 'RESI'
nu2id(14) = 'POI '
!nu2id(14) = 'POIS'
nu2id(15) = 'PM49'
nu2id(16) = 'SM49'
nu2id(17) = 'I135'
nu2id(18) = 'XE35'
nu2id(19) = 'LFP '
!nu2id(19) = 'FISP'
nu2id(20) = 'SB10'
nu2id(21) = 'H2O '
nu2id(22) = 'STM '
!nu2id(22) = 'STRM'
nu2id(23) = 'MXS '
!nu2id(23) = 'MACX'
nu2id(24) = 'CR1 '
nu2id(25) = 'CR2 '
nu2id(26) = 'CR3 '
!nu2id(24) = 'CRD1'
!nu2id(25) = 'CRD2'
!nu2id(26) = 'CRD3'
!nu2id(27) = 'TMOD'
!nu2id(27) = 'DETR'
nu2id(27) = 'DTR'

amass(1)  = 234.0000
amass(2)  = 235.0439
amass(3)  = 236.046
amass(4)  = 237.0480
amass(5)  = 238.0508
amass(6)  = 238.0500
amass(7)  = 239.0000
amass(8)  = 239.0520
amass(9)  = 240.0540
amass(10) = 241.0570
amass(11) = 242.0590
amass(12) = 243.0610

NNUCL  = 26
MNUCL  = 22
NCOMP2 = 500

END SUBROUTINE NUTYP