  subroutine dbgsp3
  
    use geom,     only : ng, nxy, nz
    use sp3senm,  only : phishp
    use sfam,     only : phif
    
    open(unit=810,file='dbgsp3_cmfd',status='unknown')
    open(unit=811,file='dbgsp3_xdir',status='unknown')
    open(unit=812,file='dbgsp3_ydir',status='unknown')
    open(unit=813,file='dbgsp3_zdir',status='unknown')
    
    do ig=1,ng
      do iz=1,nz
        write(810,815) phif(ig,1:nxy,iz)
        write(811,815) phishp(1,0,1,1:nxy,iz,ig)
        write(812,815) phishp(1,0,2,1:nxy,iz,ig)
        write(813,815) phishp(1,0,3,1:nxy,iz,ig)
      enddo
      write(810,*) 
      write(811,*) 
      write(812,*)
      write(813,*)
    enddo
    
815 format(50e20.12)    
  
  end subroutine