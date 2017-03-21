module itrinfo
    use const
    implicit none
    integer                 :: nnodal, ncmfd2g, ncmfdmg 
    real                    :: tnodal, tcmfd, tother
    real                    :: tcmfd2g       ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB              
    
    integer :: iterinner2g, iterinnermg  ! 2012_09_26 . scb
    
    contains
    
    subroutine resetitrinfo()
      nnodal  = 0
      ncmfd2g = 0
      ncmfdmg = 0
      tnodal  = 0.0
      tcmfd   = 0.0
      tcmfd2g = 0.0       ! added in ARTOS ver. 0.2 . 2012_07_20 by SCB              
      tother  = 0.0
    end subroutine
    
end module