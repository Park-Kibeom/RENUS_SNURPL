! file name =  ff.h power form functions 
      logical :: printff
      real,pointer,dimension(:,:,:,:) ::     &      !(npin-x,npin-y,ng,nffset)
                                        ff,  &      ! ff for unrodded set
                                        ffr         ! ffr for unrodded set
      common /ffc/ ff, ffr

      integer,pointer,dimension(:) ::         iffset ! (ncomp)
      common /ffc2/                         npin  ,  & ! number of pins in a row
                                            nffset,  & ! number of form function set
                                            iffset     ! ff set number for each comp.
      common /ffl/ printff