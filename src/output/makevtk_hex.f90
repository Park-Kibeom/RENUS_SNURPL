! 2014_01_17 . scb
    subroutine makevtk_hex(final)
    
      use allocs
      use vtk
      
      include 'global.h'
      include 'files.h'     
      !include 'geomh.h'
      
      integer :: nx1, ny1, nz1
      integer,save :: istep=0
      logical,save :: first=.true.
      character :: str3*3, str5*5
      character(20) :: dataname
      
      logical :: final
      
      if(first) then
        first=.false.
#ifndef CVF
        res = makedirqq('vtk/') 
#endif  
        
        call dmalloc0(iovtk,0,ng)
        iovtk(0)=109
        do ig=1,ng
          iovtk(ig)=iovtk(ig-1)+1
          write(str3,'(i3.3)') ig
#ifndef CVF          
          open(unit=iovtk(ig),file='vtk/'//trim(caseid)//'_group'//str3//'.vtk',status='unknown') 
#else
          open(unit=iovtk(ig),file=trim(caseid)//'_group'//str3//'.vtk',status='unknown')   ! 2014_08_08 . scb
#endif
          
          write(iovtk(ig),'(a)') '# vtk DataFile Version 3.0'
          
          write(iovtk(ig),'(a)') trim(caseid)//' / group '//str3//' flux'
      
          write(iovtk(ig),'(a)') 'ASCII'
          write(iovtk(ig),'(a)') 'DATASET UNSTRUCTURED_GRID'          
        enddo       
#ifndef CVF        
        open(unit=iovtk(0),file='vtk/'//trim(caseid)//'_psi.vtk',status='unknown') 
#else   
        open(unit=iovtk(0),file=trim(caseid)//'_psi.vtk',status='unknown')  ! 2014_08_08 . scb
#endif
          
        write(iovtk(0),'(a)') '# vtk DataFile Version 3.0'
          
        write(iovtk(0),'(a)') trim(caseid)//' / fission source'
      
        write(iovtk(0),'(a)') 'ASCII'
        write(iovtk(0),'(a)') 'DATASET UNSTRUCTURED_GRID'            
              
        call dmalloc0(pr,0,6,1,nxy)
        call dmalloc0(pz,0,nz)
      
      !  ix=0
      !  do ixa=1,nxa
      !    htemp=hxa(ixa)/nmeshx(ixa)
      !    do imx=1,nmeshx(ixa)
      !      ix1=ix+1
      !      px(ix1)=px(ix)+htemp  
      !      ix=ix1
      !    enddo
      !  enddo
      !
      !  iy=0
      !  do iya=1,nya
      !    htemp=hya(iya)/nmeshy(iya)
      !    do imy=1,nmeshx(iya)
      !      iy1=iy+1
      !      py(iy1)=py(iy)-htemp  
      !      iy=iy1
      !    enddo
      !  enddo
      !
      !  iz=0
      !  do iza=1,nza
      !    htemp=hza(iza)/nmeshz(iza)
      !    do imz=1,nmeshz(iza)
      !      iz1=iz+1
      !      pz(iz1)=pz(iz)+htemp  
      !      iz=iz1
      !    enddo
      !  enddo
      !
      !  call dmalloc0(indp,0,nx,0,ny,0,nz)    ! point index information
      !
      !  do ig=0,ng
      !    write(iovtk(ig),'(a6,i10,a10)') 'POINTS', np, 'float'      
      !  enddo
      !  
      !  ip=0
      !  do iz=0,nz
      !    do iy=0,ny
      !      do ix=0,nx
      !        do ig=0,ng
      !          write(iovtk(ig),'(3f15.6)') px(ix),py(iy),pz(iz)
      !        enddo
      !        
      !        indp(ix,iy,iz)=ip
      !        ip=ip+1
      !      enddo
      !    enddo
      !  enddo      
      !
      !  if(np.ne.ip)  stop 'there is some error in np dismatch !!'
      !
      !  allocate(ncinfo(nc))
      !  call dmalloc(indc,nx,ny,nz)
      !  ic=0
      !  do iz=1,nz
      !    iz1=iz-1
      !    do iy=1,ny
      !      iy1=iy-1
      !      do ix=1,nx
      !        ix1=ix-1
      !        ic=ic+1
      !        indc(ix,iy,iz)=ic
      !      
      !        call dmalloc(ncinfo(ic)%ip,8)
      !        ncinfo(ic)%ix=ix
      !        ncinfo(ic)%iy=iy
      !        ncinfo(ic)%iz=iz
      !        ncinfo(ic)%ip(1)=indp(ix1,iy,iz1)
      !        ncinfo(ic)%ip(2)=indp(ix,iy,iz1)
      !        ncinfo(ic)%ip(3)=indp(ix1,iy1,iz1)
      !        ncinfo(ic)%ip(4)=indp(ix,iy1,iz1)
      !        ncinfo(ic)%ip(5)=indp(ix1,iy,iz)
      !        ncinfo(ic)%ip(6)=indp(ix,iy,iz)
      !        ncinfo(ic)%ip(7)=indp(ix1,iy1,iz)
      !        ncinfo(ic)%ip(8)=indp(ix,iy1,iz)
      !      enddo
      !    enddo
      !  enddo           
      !
      !  if(nc.ne.ic)  stop 'there is some error in nc dismatch !!'
      !
      !  nxyz=nxy*nz
      !
      !  do ig=0,ng
      !    write(iovtk(ig),'(a5,i10,I10)') 'CELLS', nxyz, 9*nxyz
      !  enddo
      !  
      !  do iz=1,nz
      !    do ixy=1,nxy
      !      ix=ltoi(ixy)
      !      iy=ltoj(ixy)
      !    
      !      ic=indc(ix,iy,iz)
      !      do ig=0,ng
      !        write(iovtk(ig),'(9i8)') 8,ncinfo(ic)%ip
      !      enddo            
      !    enddo
      !  enddo    
      !
      !  do ig=0,ng
      !    write(iovtk(ig),*)
      !    write(iovtk(ig),'(a10,i15)') 'CELL_TYPES', nxyz
      !    do ic=1,nxyz
      !      write(iovtk(ig),'(i10)') 11
      !    enddo      
      !
      !    write(iovtk(ig),*)
      !    write(iovtk(ig),'(a9,i16)') 'CELL_DATA', nxyz
      !  enddo
      endif
      !
      !write(str5,'(i5.5)') istep*ivtkfreq
      !if(final) str5='FINAL'
      !write(dataname,'(a9, a5)') 'TIMESTEP_',str5
      !
      !do ig=1,ng                
      !  write(iovtk(ig),'(a7,8x,a20,a10,i5)') 'SCALARS', dataname, 'double', 1
      !  write(iovtk(ig),'(a12,a13)') 'LOOKUP_TABLE', 'default'
      !  do iz=1,nz
      !    do ixy=1,nxy
      !      write(iovtk(ig),'(e25.8)') phif(ig,ixy,iz)*flxlevel
      !    enddo
      !  enddo
      !enddo
      !
      !write(iovtk(0),'(a7,8x,a20,a10,i5)') 'SCALARS', dataname, 'double', 1
      !write(iovtk(0),'(a12,a13)') 'LOOKUP_TABLE', 'default'
      !do iz=1,nz
      !  do ixy=1,nxy
      !    write(iovtk(0),'(e25.8)') psi(ixy,iz)*flxlevel
      !  enddo
      !enddo      
      !            
      !istep=istep+1
  
    end subroutine
! added end    