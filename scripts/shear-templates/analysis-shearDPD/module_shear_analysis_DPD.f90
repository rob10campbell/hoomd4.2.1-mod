!############################################
subroutine corrected_temperature(nparticles,nlayerY,pos,vel,mass,typeid,Lbox,pxx,pxy,pxz,kT)
implicit none
integer, intent(in) :: nparticles, nlayerY
real(8), intent(in) :: pos(0:nparticles-1,0:2)
real(8), intent(in) :: vel(0:nparticles-1,0:2)
real(8), intent(in) :: Lbox(0:2)
real(8), intent(in) :: mass(0:nparticles-1)
integer, intent(in) :: typeid(0:nparticles-1)
!f2py depend(nparticles) pos,vel,mass,typeid
real(8), intent(out) :: pxx, pxy, pxz, kT
real(8) :: bin_width, inv_bin_width, vel_x
real(8) :: vel_bin(0:nlayerY-1)
integer :: hist_bin(0:nlayerY-1)
integer :: i, k
pxx=0.d0; pxy=0.d0; pxz=0.d0; kT=0.d0
hist_bin = 0
vel_bin = 0.d0
bin_width = Lbox(1)/nlayerY
inv_bin_width = 1.d0/bin_width

! get the solvent velocity profile
do i=0,nparticles-1
  if(typeid(i) == 0) then
    k = floor((pos(i,1)+0.5d0*Lbox(1))*inv_bin_width)
    if(k >= nlayerY) then
      k = nlayerY -1
    endif
    hist_bin(k) = hist_bin(k) + 1
    vel_bin(k) = vel_bin(k) + vel(i,0)
  endif
enddo
do i=0,nlayerY-1
  if(hist_bin(i) == 0) then
    hist_bin(i) =1
  endif
enddo
vel_bin = vel_bin/hist_bin

! remove effect of solvents
do i=0,nparticles-1
  k = floor((pos(i,1)+0.5d0*Lbox(1))*inv_bin_width)
  if(k >= nlayerY) then
    k = nlayerY -1
  endif
  vel_x = vel(i,0) - vel_bin(k)
  pxx = pxx + mass(i) * (vel_x * vel_x)
  pxy = pxy + mass(i) * (vel_x * vel(i,1))
  pxz = pxz + mass(i) * (vel_x * vel(i,2))
enddo

! recalculate temperature
kT = pxx/3.d0/nparticles
! normalize pressure by the simulation box volume
pxx = pxx/Lbox(1)**3
pxy = pxy/Lbox(1)**3
pxz = pxz/Lbox(1)**3
end subroutine corrected_temperature
!#################################################


!#################################################
subroutine vel_profile(nparticles,nlayerY,pos,vel,typeid,Lbox)
implicit none
integer, intent(in) :: nparticles, nlayerY
real(8), intent(in) :: pos(0:nparticles-1,0:2)
real(8), intent(in) :: vel(0:nparticles-1,0:2)
real(8), intent(in) :: Lbox(0:2)
integer, intent(in) :: typeid(0:nparticles-1)
!f2py depend(nparticles) pos,vel,typeid
real(8) :: bin_width, inv_bin_width
real(8) :: vel_bin(0:nlayerY-1)
integer :: hist_bin(0:nlayerY-1)
integer :: i, k
bin_width = Lbox(1)/nlayerY
inv_bin_width = 1.d0/bin_width
hist_bin = 0
vel_bin = 0.d0
do i=0,nparticles-1
  if(typeid(i) == 0) then
    k = floor((pos(i,1)+0.5d0*Lbox(1))*inv_bin_width)
    if(k >= nlayerY) then
      k = nlayerY -1
    endif
    hist_bin(k) = hist_bin(k) + 1
    vel_bin(k) = vel_bin(k) + vel(i,0)
  endif
enddo
do i=0,nlayerY-1
  if(hist_bin(i) == 0) then
    hist_bin(i) = 1
  endif
enddo
vel_bin = vel_bin/hist_bin
open(unit=14,file='vel_profile.txt',action='write')
do i=0,nlayerY-1
  write(14,'(f20.5,x,f20.5)') -0.5*Lbox(1)+(i+0.5d0)*bin_width,vel_bin(i)
enddo
close(14)
end subroutine vel_profile
!#################################################


!###########################################
subroutine fabric_tensor(nparticles,pos,radius,Lbox,cut_off,m_xy,Rxy_avg)
implicit none
integer, intent(in) :: nparticles
real(8), intent(in) :: pos(0:nparticles-1,0:2)
real(8), intent(in) :: radius(0:nparticles-1)
!f2py depend(nparticles) pos,radius
real(8), intent(in) :: Lbox(0:2)
real(8), intent(in) :: cut_off, m_xy
real(8), intent(out) :: Rxy_avg
real(8) :: dr(0:2), r1, r2, inv_Lbox(0:2), img
integer :: i, j!, coord_num
inv_Lbox = 1.d0/Lbox
!coord_num = 0
Rxy_avg = 0.d0
do i=0,nparticles-2
  do j=i+1,nparticles-1
    dr = pos(i,:)-pos(j,:)
    dr(2) = dr(2) - Lbox(2) * dnint(dr(2)*inv_Lbox(2))
    img = Lbox(1) * dnint(dr(1)*inv_Lbox(1))
    dr(1) = dr(1) - img
    dr(0) = dr(0) - img * m_xy
    dr(0) = dr(0) - Lbox(0) * dnint(dr(0)*inv_Lbox(0))
    r2 = sum(dr*dr)
    r1 = dsqrt(r2) - radius(i) - radius(j)
    if(r1 <= cut_off) then
      Rxy_avg = Rxy_avg + dr(0)*dr(1)/r2
      !coord_num = coord_num + 1
    endif
  enddo
enddo
Rxy_avg = Rxy_avg/nparticles
end subroutine fabric_tensor
!###########################################


!#########################################
subroutine density_profile(natom,nlayers,pos,typeid,Lbox,hist_z)
implicit none
integer, intent(in) :: natom, nlayers
real(8), intent(in) :: pos(0:natom-1)
real(8), intent(in) :: Lbox(0:2)
integer, intent(in) :: typeid(0:natom-1)
!f2py depend(natom) pos,vel,typeid
!f2py depend(nlayers) hist_z
real(8) :: bin_width, inv_bin_width
integer, intent(out) :: hist_z(0:nlayers-1)
integer :: i, k
bin_width = Lbox(2)/nlayers
inv_bin_width = 1.d0/bin_width
hist_z = 0
do i=0,natom-1
        if(typeid(i) == 0) then
                k = floor((pos(i)+0.5d0*Lbox(2))*inv_bin_width)
                if(k >= nlayers) then
                        k = nlayers -1
                endif
                hist_z(k) = hist_z(k) + 1
        endif
enddo
end subroutine density_profile

