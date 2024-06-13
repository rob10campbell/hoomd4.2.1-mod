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
