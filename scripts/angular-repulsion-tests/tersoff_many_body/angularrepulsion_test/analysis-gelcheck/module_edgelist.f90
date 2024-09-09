subroutine edgelist_calc(nframes,ncolloids,radii,allpos,lbox,rcut,edge_output)
implicit none

! INPUTS
integer, intent(in) :: nframes                                ! number of colloids
integer, intent(in) :: ncolloids                              ! number of colloids
real(8), intent(in) :: radii(0:ncolloids-1)                   ! colloid radii
real(8), intent(in) :: allpos(0:nframes-1,0:ncolloids-1,0:2)  ! colloid positions
!f2py depend(ncolloids) radii
!f2py depend(nframes,ncolloids) allpos
real(8), intent(in) :: lbox(0:2)                              ! simulation box size
real(8), intent(in) :: rcut                                   ! surface-surface cut-off for a bond
character*20, intent(in) :: edge_output                       ! output filepath for edgelist

! OUTPUTS
!N/A

! INTERNAL VARIABLES
real(8) :: rij(0:2), hij                 ! center-center (rij) and surface-surface (hij) interaction distance
integer :: f                             ! frame index, f
integer :: i,j                           ! colloid indices i, j
character(len=10) :: file_id_frame       ! variable for frame number as string
character(len=50) :: filename_edge       ! variable for frame-specific filepath/name


! for each frame
do f=0,nframes-1
  ! write the frame number into a string for naming files
  write(file_id_frame, '(i0)') f

  ! construct the edgelist filename
  filename_edge = trim(adjustl(edge_output)) // trim(adjustl(file_id_frame)) // '.csv'
  ! create the file for this frame
  open(unit=14,file=trim(filename_edge),action='write')
  ! write the header
  write(14,fmt="(*(g0:','))") "i","j"

  ! caluclate the edge list, considering for all possible pairs
  do i=0,ncolloids-2
    do j=i+1,ncolloids-1
     ! calculate the center-center interparticle distance
     rij(0) = allpos(f,i,0) - allpos(f,j,0)
     rij(1) = allpos(f,i,1) - allpos(f,j,1)
     rij(2) = allpos(f,i,2) - allpos(f,j,2)
     ! adjust for interactions across the periodic boundaries
     ! NOTE: assuming that there is NO shear deforming the box
     rij = rij - lbox * dnint(rij/lbox)
     ! conver to surface surface distance
     hij = dsqrt(sum(rij*rij)) - (radii(i) + radii(j))
     ! if it is a bond
     if(hij<=rcut) then
       ! record in the edge list
       write(14,fmt="(*(g0:','))") i,j
     endif
   enddo
  enddo

  ! close the file
  close(14)

enddo

end subroutine edgelist_calc
