! 2012_09_07 . scb
  !subroutine initbdf(unit_tm, bdforder, ng, nxy, nz)   ! 2013_05_22 . scb
  subroutine initbdf(unit_tm, bdforder, ng, nxy, nz, ifsp3)   ! 2013_05_22 . scb
    
    use allocs
    use bdf
    !use geomhex,  only : ifhexsp3   ! 2012_12_06 . scb
      
    !include 'geomh.h'   ! 2013_05_22 . scb
      
    real :: unit_tm
    integer :: bdforder, ng, nxy, nz
    logical :: ifsp3   ! 2013_05_22 . scb
      
    integer,parameter :: ng2=2
      
    unittm=unit_tm
    bdforder_mod=bdforder
    bdforder0=bdforder    ! 2014_09_04 . scb
!need power update? => is used for step size control method
    bdfordersave=bdforder
           
    !call dmalloc0(bdfarray%bdfcoef,0,bdforder)        
    call dmalloc0(bdfarray%bdfcoef,0,5)     ! 2014_12_29 . scb
    bdfcoef => bdfarray%bdfcoef
! 2013_05_16 . scb      
    !call dmalloc0(bdfarray%phibdf,1,ng2,0,nxy,0,nz,1,bdforder)
    !phibdf => bdfarray%phibdf                
    !if (ng.ne.ng2) then
    !  call dmalloc0(bdfarray%phifbdf,1,ng,0,nxy,0,nz,1,bdforder)
    !  phifbdf => bdfarray%phifbdf
    !endif
    !call dmalloc0(bdfarray%phifbdf,1,ng,0,nxy,0,nz,1,bdforder)
    call dmalloc0(bdfarray%phifbdf,1,ng,0,nxy,0,nz,1,5)   ! 2014_12_29 . scb
    phifbdf => bdfarray%phifbdf              
    if (ng.ne.ng2) then
      !call dmalloc0(bdfarray%phibdf,1,ng2,0,nxy,0,nz,1,bdforder)
      call dmalloc0(bdfarray%phibdf,1,ng2,0,nxy,0,nz,1,5)    ! 2014_12_29 . scb
      phibdf => bdfarray%phibdf
    else
      phibdf => phifbdf
    endif
! added end      
      
! 2012_12_06 . scb      
    !if(ifhexsp3) then
    if(ifsp3) then  ! 2013_08_12 . scb
      !call dmalloc0(bdfarray%phibdf2,1,ng2,0,nxy,0,nz,1,bdforder)
      !phibdf2 => bdfarray%phibdf2          
      !if (ng.ne.ng2) then
      !  call dmalloc0(bdfarray%phifbdf2,1,ng,0,nxy,0,nz,1,bdforder)
      !  phifbdf2 => bdfarray%phifbdf2
      !endif        
      !call dmalloc0(bdfarray%phifbdf2,1,ng,0,nxy,0,nz,1,bdforder)
      call dmalloc0(bdfarray%phifbdf2,1,ng,0,nxy,0,nz,1,5)   ! 2014_12_29 . scb
      phifbdf2 => bdfarray%phifbdf2          
      if (ng.ne.ng2) then
        !call dmalloc0(bdfarray%phibdf2,1,ng,0,nxy,0,nz,1,bdforder)
        call dmalloc0(bdfarray%phibdf2,1,ng,0,nxy,0,nz,1,5)   ! 2014_12_29 . scb
        phibdf2 => bdfarray%phibdf2
      else
        phibdf2 => phifbdf2
      endif    
    endif
! added end      
    
    !call dmalloc(deltmarray,bdforder)
    call dmalloc(deltmarray,5)   ! 2014_12_29 . scb
      
! 2014_09_18 . scb   
    !flagsrcdt=.true.
    if(flagsrcdt) then 
      if(ifsp3) then
        
      else
        
        !call dmalloc0(srcdtbdf,1,ng,1,nxy,1,nz,1,3,0,bdforder)
        call dmalloc0(srcdtbdf,1,ng,1,nxy,1,nz,1,3,0,5)   ! 2014_12_29 . scb
      
        srcdt => srcdtbdf(:,:,:,:,0)
        call dmalloc(spsrcdt,ng,nxy,nz,3)
      endif
      
    elseif(flagsrcdt2) then
      if(ifsp3) then
        !call dmalloc0(dtcff2,0,2,1,2,1,ng,1,nxy,1,nz,1,3,0,bdforder)   ! 2014_10_06 . scb
        call dmalloc0(dtcff2,0,2,1,2,1,ng,1,nxy,1,nz,1,3,0,5)   ! 2014_12_29 . scb
      else
        !call dmalloc0(dtcff,0,2,1,ng,1,nxy,1,nz,1,3,0,bdforder)
        call dmalloc0(dtcff,0,2,1,ng,1,nxy,1,nz,1,3,0,5)   ! 2014_12_29 . scb
      endif      
    endif    
! added end      
      
  end subroutine