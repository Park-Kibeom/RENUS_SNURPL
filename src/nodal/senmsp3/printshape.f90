! 2015_01_14 . scb
  subroutine printshape
  
    use geom, only : ng
    use senm1d
  
    real :: k=0.d0, sinhk(2), coshk(2)
  
    do ih=ns,ne
      do ig=1,ng
        write(114,*) ig, 1, ih
        write(114,*) h(ih)
        write(114,*) phin(:,0,ih,ig)
        
        write(114,*) sqrt(ksq(1,ih,ig)), sqrt(ksq(2,ih,ig))
        
        write(114,*) s(1:2, ih, ig)
        write(114,*) s(3:4, ih, ig)
      
        write(114,*) mA(:,ih,ig)
        write(114,*) mB(:,ih,ig)
        do io=0,4
          write(114,*) mpsol(io,1,ih,ig), mpsol(io,2,ih,ig)
        enddo   
        
        write(114,*) xsnfn(ih,ig), xsdn(ih,ig)
      enddo
    enddo
  
    close(114)
    
  end subroutine

  