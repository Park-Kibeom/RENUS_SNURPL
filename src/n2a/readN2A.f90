module readN2A        ! 2016. 9. 5. ~   jjh
  
  use MASTERXSL, only : NBURN, NDER, NUCNUM, cnum1, IVER, basexs, ntype2, KAPPA_T, &
                        burnstep, derstep, boronstep, tfuelstep, tmodstep, modstep, nboron, ntfuel, ntmod, ndmod, &
                        ppmxs, tmxs, tfxs, dmxs, LFM_coeff
  implicit none

  integer, parameter :: charL=256, nbf=8
  integer, parameter :: ngroup=2
  integer, parameter :: nufis=1, fs=2, cap=3, tr=4, sct=5, ab=6, kap=7, MACX_ID=23
  integer, parameter :: brn=1, tf=2, tm=3, dm=4
  
  character*charL :: inputfilename
  logical :: Ladf2, Linterpol, Lautomsk
  real(nbf) :: boron, tmod, dmod, tfuel, pressure, void
!  real(nbf), allocatable :: burnupstep(:), branchstep(:)
!  real(nbf), allocatable :: boronstep1(:), tfuelstep1(:), tmodstep1(:), dmodstep1(:)       ! 2016. 12.20. jjh
  INTEGER::MAXBURN,MAXDER,MAXDMOD          ! 2016. 12.23. jjh
  INTEGER::MAXBORON,MAXTFUEL,MAXTMOD       ! 2016. 12.23. jjh
  INTEGER::MAXNB

  integer :: nassm, isonum(NUCNUM)
  character(3) :: isoXS(NUCNUM)
  integer, allocatable :: assmnum(:)
  character*charL, allocatable :: assmlist(:)

!  type XS_type
!    real(nbf), allocatable :: BaseXS(:,:)          ! (NBURN, group)
!    real(nbf), allocatable :: dXS(:,:,:,:)  ! (NDER(ICOMP), group, nbranch, branchtype(1:boron, 2:tfuel, 3:tmod, 4:dmod))
!  end type XS_type
    
  type Nuclei_XS
    character(4)           :: IsoName
    integer                :: IsoNum
    character(3)           :: XStype
    real(nbf), allocatable :: ND(:)      ! (NBURN)
!    real(nbf), allocatable :: ND_Branch(:,:)      ! (NDER(ICOMP), IVER*nbranch)
!    type(XS_type) :: Xsec(NTYPE2)
  end type Nuclei_XS
  
  type assmtype
    character*charL        :: CaseID
    real(nbf), allocatable :: ADF(:,:)      ! (NBURN, ngroup)
    real(nbf), allocatable :: dADF(:,:,:,:)      ! (NDER(ICOMP), ngroup, nbranch, IVER)
    type(nuclei_XS)        :: Nuclei(NUCNUM)
  end type assmtype

  type(assmtype), allocatable :: assm(:)  
   
contains

  subroutine readN2Aout(indev, icomp)
    implicit none

    logical :: First_read1=.true., First_read2=.true.
    integer :: indev, icomp
    character*charL :: fname, temp, dum, dum2
    integer :: a,b, IO=11
    integer :: nline=0, i, j, k,g,m, t, s, xs, br, n
    real :: keff(61)
    
    
    isonum( 1)=0    ; isonum( 2)=92235; isonum( 3)=92236; isonum( 4)=93237; isonum( 5)=92238
    isonum( 6)=94238; isonum( 7)=93239; isonum( 8)=94239; isonum( 9)=94240; isonum(10)=94241
    isonum(11)=94242; isonum(12)=95243; isonum(13)=0    ; isonum(14)=0    ; isonum(15)=61649
    isonum(16)=62649; isonum(17)=53635; isonum(18)=54635; isonum(19)=0    ; isonum(20)=5000
    isonum(21)=0    ; isonum(22)=0    ; isonum(23)=0    ; isonum(24)=0    ; isonum(25)=0   
    isonum(26)=0    ; isonum(27)=0 !   ; isonum(28)=0    
    

    isoXS( 1)='MIC'; isoXS( 2)='MIC'; isoXS( 3)='MIC'; isoXS( 4)='MIC'; isoXS( 5)='MIC'
    isoXS( 6)='MIC'; isoXS( 7)='MIC'; isoXS( 8)='MIC'; isoXS( 9)='MIC'; isoXS(10)='MIC'
    isoXS(11)='MIC'; isoXS(12)='MIC'; isoXS(13)='MAC'; isoXS(14)='MAC'; isoXS(15)='MIC'
    isoXS(16)='MIC'; isoXS(17)='MIC'; isoXS(18)='MIC'; isoXS(19)='MAC'; isoXS(20)='MIC'
    isoXS(21)='MIC'; isoXS(22)='MAC'; isoXS(23)='MAC'; isoXS(24)='MAC'; isoXS(25)='MAC'   
    isoXS(26)='MAC'; isoXS(27)='MIC'
!    isoXS(26)='MAC'; isoXS(27)='MAC'; isoXS(28)='MIC'
    
    !isonum( 1)=90232 ; isonum( 2)=91233 ; isonum( 3)=92233 ; isonum( 4)=92234
    !isonum( 5)=92235 ; isonum( 6)=92236 ; isonum( 7)=92238 ; isonum( 8)=93237
    !isonum( 9)=93239 ; isonum(10)=94248 ; isonum(11)=94249 ; isonum(12)=94240
    !isonum(13)=94241 ; isonum(14)=94242 ; isonum(15)=95241 ; isonum(16)=95242
    !isonum(17)=95243 ; isonum(18)=96242 ; isonum(19)=96244 ; isonum(20)=0
    !isonum(21)=61649 ; isonum(22)=62649 ; isonum(23)=53635 ; isonum(24)=54635
    !isonum(25)=0     ; isonum(26)=0     ; isonum(26)=0     ; isonum(28)=5000
    !isonum(29)=0     ; isonum(30)=0     ; isonum(31)=0     ; isonum(32)=0     
    !isonum(33)=0     ; isonum(34)=0     
    !
    !isoXS( 1)='MIC' ; isoXS( 2)='MIC' ; isoXS( 3)='MIC' ; isoXS( 4)='MIC'
    !isoXS( 5)='MIC' ; isoXS( 6)='MIC' ; isoXS( 7)='MIC' ; isoXS( 8)='MIC'
    !isoXS( 9)='MIC' ; isoXS(10)='MIC' ; isoXS(11)='MIC' ; isoXS(12)='MIC' 
    !isoXS(13)='MIC' ; isoXS(14)='MIC' ; isoXS(15)='MIC' ; isoXS(16)='MIC'
    !isoXS(17)='MIC' ; isoXS(18)='MIC' ; isoXS(19)='MIC' ; isoXS(20)='MAC'
    !isoXS(21)='MIC' ; isoXS(22)='MIC' ; isoXS(23)='MIC' ; isoXS(24)='MIC'
    !isoXS(25)='MIC' ; isoXS(26)='MAC' ; isoXS(27)='MAC' ; isoXS(28)='MIC'
    !isoXS(29)='MAC' ; isoXS(30)='MAC' ; isoXS(31)='MAC' ; isoXS(32)='MAC'
    !isoXS(33)='MAC' ; isoXS(34)='MIC'   
    
   nassm=0
   fname=trim(inputfilename)
   open(io, file=fname, readonly)
   do
     read(io, '(a256)') temp
     if (temp(1:10) /='          ') then
       read(temp, *) dum
       if (dum(1:10)=='comp_index') then
!         read(temp, *) dum, nassm
         nassm=nassm+1
      end if
     end if
     
     if (temp(1:1)=='.') exit       
   end do
   close(io)
   
   if(first_read1==.true.) then
     allocate(assmlist(cnum1),assmnum(nassm))
     assmnum=0
     first_read1=.false.
   
   open(io, file=fname, readonly)
   t=0
   do
     read(io, '(a256)') temp
     if (temp(1:10) /='          ') then
       read(temp, *) dum
       if (dum(1:10)=='comp_index') then
          t=t+1
          read(temp, *) dum, assmnum(t), assmlist(assmnum(t))
       end if
     end if
     if (temp(1:1)=='.') exit       
   end do
   close(io)   
   
  do i=1, nassm
    
    open(io, file=assmlist(assmnum(i)), readonly)  
    do 
      read(io, '(a256)' ) temp
      call Line_upper(temp)
      if (temp(1:1)=='.') then
        exit
      else if (temp(2:14)=='ASSEMBLY TYPE') then
        read(temp,*) dum, dum, dum,assm(assmnum(i))%caseID
      else if (temp(3:13)=='BURNUP STEP') then
        read(temp,*) dum, dum, dum, nburn(assmnum(i)) 
      else if (temp(3:13)=='BRANCH STEP') then
        read(temp,*) dum, dum, dum, nder(assmnum(i)) 
      else if (temp(3:12)=='BORON STEP') then
        read(temp,*) dum, dum, dum, nboron(assmnum(i))
      else if (temp(3:12)=='TFUEL STEP') then  
        read(temp,*) dum, dum, dum, ntfuel(assmnum(i))
      else if (temp(3:11)=='TMOD STEP') then  
        read(temp,*) dum, dum, dum, ntmod(assmnum(i))        
      else if (temp(3:11)=='DMOD STEP') then  
        read(temp,*) dum, dum, dum, ndmod(assmnum(i))
      else if (temp(2:10)=='LFM_COEFF') then    
        read(io, *) dum    !nufis=1, fs=2, cap=3, tr=4, sct=5, ab=6, kap=7
        read(io, *) dum, LFM_coeff(assmnum(i), 1, tr, 1:2)
        read(io, *) dum, LFM_coeff(assmnum(i), 2, tr, 1:2)
        read(io, *) dum, LFM_coeff(assmnum(i), 1, ab, 1:2)
        read(io, *) dum, LFM_coeff(assmnum(i), 2, ab, 1:2)
        read(io, *) dum, LFM_coeff(assmnum(i), 1, nufis, 1:2)
        read(io, *) dum, LFM_coeff(assmnum(i), 2, nufis, 1:2)
        read(io, *) dum, LFM_coeff(assmnum(i), 1, sct, 1:2)
        read(io, *) dum, LFM_coeff(assmnum(i), 2, sct, 1:2)
        read(io, *) dum, LFM_coeff(assmnum(i), 1, kap, 1:2)
        read(io, *) dum, LFM_coeff(assmnum(i), 2, kap, 1:2)
        LFM_coeff(assmnum(i),:,:,:)=LFM_coeff(assmnum(i),:,:,:)/100
        !exit
      end if
    end do
    close(io)
  end do 

MAXBURN=NBURN(1); MAXDER=NDER(1); MAXDMOD=NDMOD(1)
maxboron=nboron(1); maxtfuel=ntfuel(1); maxtmod=ntmod(1)
DO I=1,cnum1
    IF(NBURN(I).GT.MAXBURN) MAXBURN=NBURN(I)
    IF(NDER(I).GT.MAXDER) MAXDER=NDER(I)
    IF(Nboron(I).GT.MAXboron) MAXboron=Nboron(I)
    IF(Ntfuel(I).GT.MAXtfuel) MAXtfuel=Ntfuel(I)
    IF(Ntmod(I).GT.MAXtmod) MAXtmod=Ntmod(I)
    IF(NDMOD(I).GT.MAXDMOD) MAXDMOD=NDMOD(I)    
END DO
 end if   
   
!  io=indev
!  t=icomp     
    
  ! read N2A's output
  open(io, file=assmlist(icomp))
  read(io, *) dum
  read(io, *) dum
  read(io, *) dum
  read(io, *) dum
  read(io, *) NBURN(icomp), nder(icomp), ndmod(icomp), Ladf2, Linterpol, Lautomsk    ! 2016. 12.20. jjh
  read(io, *) dum
  read(io, *) boron, tmod, dmod, tfuel, pressure                                          ! 2016. 12.20. jjh
  read(io, '(1/)')
  
  maxnb=max(maxboron, maxtfuel, maxtmod, maxdmod)
  allocate(assm(icomp)%ADF(maxburn,ngroup), assm(icomp)%dADF(maxder,ngroup,maxnb, IVER))  
  assm(icomp)%ADF=1._8; assm(icomp)%dADF=0._8
  
  b=maxnb*maxder
    if (First_read2==.true.) then
    ALLOCATE(BASEXS(CNUM1,NTYPE2,MAXBURN,NUCNUM,ngroup))
    ALLOCATE(KAPPA_T(cnum1,NUCNUM,ngroup)) 
    ALLOCATE(BORONSTEP(CNUM1,MAXBORON), TFUELSTEP(CNUM1,MAXTFUEL), TMODSTEP(CNUM1,MAXTMOD))
    ALLOCATE(PPMXS(CNUM1,NTYPE2,MAXDER,NUCNUM,ngroup,MAXBORON))
    ALLOCATE(TFXS(CNUM1,NTYPE2,MAXDER,NUCNUM,ngroup,MAXTFUEL))
    ALLOCATE(TMXS(CNUM1,NTYPE2,MAXDER,NUCNUM,ngroup,MAXTMOD))
    ALLOCATE(DMXS(CNUM1,NTYPE2,MAXDER,NUCNUM,ngroup,MAXDMOD))
    
    BASEXS(:,:,:,:,:) = 0.d0;  KAPPA_T(:,:,:)    = 0.d0
    boronstep(:,:)    = 0.d0;  tfuelstep(:,:)    = 0.d0 ; tmodstep(:,:)     = 0.d0
    PPMXS(:,:,:,:,:,:)  = 0.d0;TFXS(:,:,:,:,:,:)   = 0.d0;TMXS(:,:,:,:,:,:)   = 0.d0;DMXS(:,:,:,:,:,:) = 0.d0
    First_read2=.false.
  end if
  
  
  do k=1, NUCNUM
    allocate(assm(icomp)%nuclei(k)%ND(maxburn)) !, assm(icomp)%nuclei(k)%ND_branch(b,IVER))
    assm(icomp)%nuclei(k)%ND=0._8
!    assm(icomp)%nuclei(k)%ND_branch=0._8 

!    do j=1, NTYPE2   ! NTYPE2=7
!      allocate(assm(icomp)%nuclei(k)%xsec(j)%basexs(maxburn, ngroup))
!      allocate(assm(icomp)%nuclei(k)%xsec(j)%dxs(maxder, ngroup, maxnb, IVER))
!      assm(icomp)%nuclei(k)%xsec(j)%basexs=0._8; assm(icomp)%nuclei(k)%xsec(j)%dXS=0._8
!    end do
  end do
    
  
  do i=1, NBURN(icomp)
    read(io, *) Burnstep(icomp,i)
!    Burnstep(icomp,i)=burnupstep(i)
  end do
  read(io, '(1/)')
  
  
  do i=1, NDER(ICOMP)
    read(io, *) Derstep(icomp,i)
!    Derstep(icomp,i)=branchstep(i)
  end do
  read(io, '(1/)')
  
  ! 2016. 12.20. jjh 
  do i=1, nboron(icomp)
    read(io, *) boronstep(icomp,i)
  end do
  read(io, '(1/)')  
  
  do i=1, ntfuel(icomp)
    read(io, *) tfuelstep(icomp,i)
  end do
  read(io, '(1/)')  
  
  do i=1, ntmod(icomp)
    read(io, *) tmodstep(icomp,i)
  end do
  read(io, '(1/)')  
 ! end add
  
  
  do i=1, ntmod(icomp)
    read(io, *) modstep(icomp,i)
  end do
  read(io, '(1/)')
  
  read (io, *) dum, g
    do i=1, NBURN(icomp)
      read(io, *) assm(icomp)%ADF(i,1:g)
    end do
  read(io, '(1/)')
  
  do s=1, nboron(icomp)
    read (io, *) dum, dum, dum, m
    read (io, *) dum, g
    do i=1, NDER(ICOMP)
      read(io, *) assm(icomp)%dADF(i,1:g, s,brn)
    end do
    read(io, *) 
  end do
  read(io, *)
  
  do s=1, ntfuel(icomp)
    read (io, *) dum, dum, dum, m
    read (io, *) dum, g
    do i=1, NDER(ICOMP)
      read(io, *) assm(icomp)%dADF(i,1:g, s,tf)
    end do
    read(io, *)    
  end do
  read(io, *)
  
  do s=1, ntmod(icomp)
    read (io, *) dum, dum, dum, m
    read (io, *) dum, g
    do i=1, NDER(ICOMP)
      read(io, *) assm(icomp)%dADF(i,1:g, s,tm)
    end do
    read(io, *)  
  end do
  read(io, *)
  
  do s=1, ntmod(icomp)
    read (io, *) dum, dum, dum, m
    read (io, *) dum, g
    do i=1, NDER(ICOMP)
      read(io, *) assm(icomp)%dADF(i,1:g, s,dm)
    end do
  end do
  read(io, '(2/)') 

  do k=1, NUCNUM
    read(io, '(a4)') assm(icomp)%nuclei(k)%isoname
    assm(icomp)%nuclei(k)%isonum=isonum(k)
    assm(icomp)%nuclei(k)%XStype=isoXS(k)
    do j=1, ngroup
      read(io, *) dum, dum, g
      read(io, *) dum
      do i=1, NBURN(icomp)
        read(io, *) a, assm(icomp)%nuclei(k)%ND(i), basexs(icomp,1:NTYPE2-1,i,k,g), KAPPA_T(icomp,k,g)
!        read(io, *) a, assm(icomp)%nuclei(k)%ND(i), basexs(icomp,nufis,i,k,g), basexs(icomp,fs,i,k,g), basexs(icomp,cap,i,k,g), &
!                    basexs(icomp,tr,i,k,g), basexs(icomp,sct,i,k,g), basexs(icomp,ab,i,k,g), KAPPA_T(icomp,k,g)
!        read(io, *) a, assm(icomp)%nuclei(k)%ND(i), assm(icomp)%nuclei(k)%xsec(nufis)%basexs(i, g), assm(icomp)%nuclei(k)%xsec(fs)%basexs(i, g), &
!                    assm(icomp)%nuclei(k)%xsec(cap)%basexs(i, g), assm(icomp)%nuclei(k)%xsec(tr)%basexs(i, g), assm(icomp)%nuclei(k)%xsec(sct)%basexs(i, g), &
!                    assm(icomp)%nuclei(k)%xsec(ab)%basexs(i, g), assm(icomp)%nuclei(k)%xsec(kap)%basexs(i, g)
      end do
      
      do s=1, nboron(icomp)
        read(io, *) dum, dum, m 
        read(io, *) dum 
        do i=1, NDER(ICOMP)
          read(io, *) a,dum, ppmxs(icomp,1:ntype2,i,k,g,m)
!          read(io, *) a, assm(icomp)%nuclei(k)%ND_branch(i+(s-1)*nboron(icomp),brn), assm(icomp)%nuclei(k)%xsec(nufis)%dXS(i,g,s,brn), assm(icomp)%nuclei(k)%xsec(fs)%dXS(i,g,s,brn), &
!                      assm(icomp)%nuclei(k)%xsec(cap)%dXS(i,g,s,brn), assm(icomp)%nuclei(k)%xsec(tr)%dXS(i,g,s,brn), assm(icomp)%nuclei(k)%xsec(sct)%dXS(i,g,s,brn), &
!                      assm(icomp)%nuclei(k)%xsec(ab)%dXS(i,g,s,brn), assm(icomp)%nuclei(k)%xsec(kap)%dXS(i,g,s,brn)
        end do  
      end do  ! s : nboron(icomp)
      
      do s=1, ntfuel(icomp)
        read(io, *) dum, dum, m 
        read(io, *) dum 
        do i=1, NDER(ICOMP)
          read(io, *) a, dum, tfxs(icomp,1:ntype2,i,k,g,m)
!          read(io, *) a, assm(icomp)%nuclei(k)%ND_branch(i+(s-1)*ntfuel(icomp),tf), assm(icomp)%nuclei(k)%xsec(nufis)%dXS(i,g,s,tf), assm(icomp)%nuclei(k)%xsec(fs)%dXS(i,g,s,tf), &
!                      assm(icomp)%nuclei(k)%xsec(cap)%dXS(i,g,s,tf), assm(icomp)%nuclei(k)%xsec(tr)%dXS(i,g,s,tf), assm(icomp)%nuclei(k)%xsec(sct)%dXS(i,g,s,tf), &
!                      assm(icomp)%nuclei(k)%xsec(ab)%dXS(i,g,s,tf), assm(icomp)%nuclei(k)%xsec(kap)%dXS(i,g,s,tf)
        end do  
      end do  ! s : ntfuel(icomp)
        
      do s=1, ntmod(icomp)
        read(io, *) dum, dum, m 
        read(io, *) dum         
        do i=1, NDER(ICOMP)
          read(io, *) a, dum, tmxs(icomp,1:ntype2,i,k,g,m)
!          read(io, *) a, assm(icomp)%nuclei(k)%ND_branch(i+(s-1)*ntmod(icomp),tm), assm(icomp)%nuclei(k)%xsec(nufis)%dXS(i,g,s,tm), assm(icomp)%nuclei(k)%xsec(fs)%dXS(i,g,s,tm), &
!                      assm(icomp)%nuclei(k)%xsec(cap)%dXS(i,g,s,tm), assm(icomp)%nuclei(k)%xsec(tr)%dXS(i,g,s,tm), assm(icomp)%nuclei(k)%xsec(sct)%dXS(i,g,s,tm), &
!                      assm(icomp)%nuclei(k)%xsec(ab)%dXS(i,g,s,tm), assm(icomp)%nuclei(k)%xsec(kap)%dXS(i,g,s,tm)
        end do  
      end do  ! s : ntmod(icomp)        
      
      do s=1, ndmod(icomp)
        read(io, *) dum, dum, m 
        read(io, *) dum         
        do i=1, NDER(ICOMP)
          read(io, *) a, dum, dmxs(icomp,1:ntype2,i,k,g,m)         
!          read(io, *) a, assm(icomp)%nuclei(k)%ND_branch(i+(s-1)*ndmod(icomp),dm), assm(icomp)%nuclei(k)%xsec(nufis)%dXS(i,g,s,dm), assm(icomp)%nuclei(k)%xsec(fs)%dXS(i,g,s,dm), &
!                      assm(icomp)%nuclei(k)%xsec(cap)%dXS(i,g,s,dm), assm(icomp)%nuclei(k)%xsec(tr)%dXS(i,g,s,dm), assm(icomp)%nuclei(k)%xsec(sct)%dXS(i,g,s,dm), &
!                      assm(icomp)%nuclei(k)%xsec(ab)%dXS(i,g,s,dm), assm(icomp)%nuclei(k)%xsec(kap)%dXS(i,g,s,dm)
        end do  
      end do  ! s : ndmod(icomp)        
    end do  ! j : ngroup
    read(io, *) dum
  end do   !  k : NUCNUM
 
  close(io)
 
  call N2A_SCANXSL(indev, icomp)
  
  
  !=================================================
    ! keff
  
    !  nufis=1, fs=2, cap=3, tr=4, sct=5, abs=6, kap=7
!    t=icomp
!    do k =1, NBURN(1)   
!       keff(k) = (assm(icomp)%nuclei(MACX_ID)%xsec(nufis)%basexs(k,1) &
!                  * (assm(icomp)%nuclei(MACX_ID)%xsec(ab)%basexs(k,2) + assm(icomp)%nuclei(MACX_ID)%xsec(sct)%basexs(k,2)) &
!                + assm(icomp)%nuclei(MACX_ID)%xsec(nufis)%basexs(k,2) * assm(icomp)%nuclei(MACX_ID)%xsec(sct)%basexs(k,1) ) &
!                / (assm(icomp)%nuclei(MACX_ID)%xsec(ab)%basexs(k,1) &
!                  * (assm(icomp)%nuclei(MACX_ID)%xsec(ab)%basexs(k,2) + assm(icomp)%nuclei(MACX_ID)%xsec(sct)%basexs(k,2)) &
!                   + assm(icomp)%nuclei(MACX_ID)%xsec(ab)%basexs(k,2) * assm(icomp)%nuclei(MACX_ID)%xsec(sct)%basexs(k,1))
!    end do    
  
  
  
  
  
  
    
  end subroutine readN2Aout

  
  
end module readN2A