      common /ipertbs/ npbank, & ! number of banks that changes positions
                       maxntm    ! max number of perturb times of banks
      
      integer, pointer, dimension(:) :: idpbank, &   !(npbanks), id of banks that changes positions
                                        ntpbank      !(npbanks), number of perturb times of each banks
      common /ipertb/ idpbank, ntpbank 
     
      real, pointer, dimension(:,:,:) :: tbank    !(nmaxtmbank,2,npbank, perturb times of each bank
      common /dperrb/ tbank