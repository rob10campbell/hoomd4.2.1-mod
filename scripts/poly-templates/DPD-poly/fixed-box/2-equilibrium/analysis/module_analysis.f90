! Fortran module for analyzing a colloid sim GSD file
! NOTE: requires compilation with compile-module
! NOTE: is run from a matching sim-analysis.py Python script
! NOTE: ASSUMES 1 colloid type and 1 solvent type
! (Rob Campbell)

! Calculates:
! * average coordination number (Z) for each frame
! * mean squared displacement (MSD)
! * pair correlation function (PCF) g(r1, r2)
! * radial distribution function (RDF) g(r)

  ! NOTE: f2py only understands limited number of KIND parameters
  ! use real(kind=8) not real64

  ! NOTE: fortran indexes from 1 by default, switch to base 0 to match C/Python
  ! for indexing from 0, 0:ncolloid-1 is of length ncolloid, 0:2 is length 3

  ! NOTE: must tell f2py about any variable dependencies or you will get Value Errors


!###########################################
subroutine coordination_number(Lbox,ncolloid,R_C1,pos,m_xy,cut_off,Zavg)
! calculates average coordination number (Z), AKA contact number

  implicit none

  !!! INPUTS
  ! the simulation box size (L_X, L_Y, L_Z)
  real(8), intent(in) :: Lbox(0:2)

  ! number of colloid particles
  integer, intent(in) :: ncolloid 
  ! the radius of each colloid particle
  real(8), intent(in) :: R_C1(0:ncolloid-1)
  ! xyz position of all colloid particles
  real(8), intent(in) :: pos(0:ncolloid-1,0:2)

  ! tell f2py about dependencies (required as commented out)
  !f2py depend(ncolloid) pos,R_C1
 
  ! the cut-off distance for defining a bond, and
  ! the xy tilt factor (box deformation from shear in x-direction)
  real(8), intent(in) :: m_xy, cut_off

  !!! OUTPUTS
  ! average contact number
  real(8), intent(out) :: Zavg

  !!! INTERNAL VALUES
  ! the center-center interparticle distance (dr), 
  ! the surface-surface interparticle distance (h_ij), 
  ! 1/box-size, and the modification for
  ! for particles interacting across a boundary (img)   
  real(8) :: dr(0:2), h_ij, inv_Lbox(0:2), img
  ! particle tags (i, j) and coordination number (Z)
  integer :: i, j, Z

  !!! SET INITIAL VALUES
  ! calculate 1/box-size
  inv_Lbox = 1.d0/Lbox
  ! set the initial contact number to zero
  Z = 0

  !!! calculate Z for each colloid particle
  ! for all colloid particles
  do i=0,ncolloid-1
        ! check it's interaction with each other colloid particle
        do j=0,ncolloid-1
                ! calculate the center-center interparticle distance
                dr = pos(i,:)-pos(j,:) 
                ! and adjust the calculation for particles interacting
                ! across a boundary (x=0, y=1, z=2), including box tilt
                ! from any shear in the x-direction (m_xy)
                dr(2) = dr(2) - Lbox(2) * dnint(dr(2)*inv_Lbox(2))
                img = Lbox(1) * dnint(dr(1)*inv_Lbox(1))
                dr(1) = dr(1) - img
                dr(0) = dr(0) - img * m_xy
                dr(0) = dr(0) - Lbox(0) * dnint(dr(0)*inv_Lbox(0)) 
                ! update coordinates and convert to 
                ! surface-surface interparticle distance
                h_ij = dsqrt(sum(dr*dr)) - R_C1(i) - R_C1(j)
                ! if the 2nd particle is within the cut-off distance
                if(h_ij <= cut_off) then
                        ! increase the coordination number by 1
                        Z = Z + 1
                endif
        enddo
  enddo

  ! caluclate the average coordination number
  Z = Z - ncolloid
  Zavg = float(Z)/ncolloid

end subroutine coordination_number


!###########################################
subroutine msd_calculation(nframe,ncolloid,allpos,t1)
! calculate the mean-squared displacement of all colloid particles

  implicit none

  !!! INPUTS
  ! number of frames, total number of colloids
  integer, intent(in) :: nframe, ncolloid
  ! xyz positions of all colloids in all frames
  real(8), intent(in) :: allpos(0:nframe-1,0:ncolloid-1,0:2) 
  ! timestep converstion factor (t1 = period * dt_Integration)
  integer, intent(in) :: t1

  ! tell f2py about dependencies (required as commented out)
  !f2py depend(nframe,ncolloid) allpos
  !f2py depend(ncolloid) typeid

  !!! OUTPUTS
  ! N/A (output file created)

  !!! INTERNAL VALUES
  ! the center-center change in position (dpos) of each
  ! particle, and it's value squared (dpos2)
  real(8) :: dpos(0:ncolloid-1,0:2), dpos2(0:ncolloid-1)
  ! frame tag (i)
  integer :: i
  ! mean squared displacement, it's corrected sample variance,
  ! and it's (corrected) sample standard deviation
  real(8) :: msd_value, msd_corrsmplvar, msd_smplstddev

  !!! SET INITIAL VALUES
  ! N/A

  !!! create a file to record the MSD and standard deviation
  open(unit=14,file='msd.txt',action='write') 
  ! label the output data values
  write(14,'(A,x,A,x,A)') 'DPD-time', 'msd-colloids', 'msd_smplstddev' 

  !!! calculate the msd for each frame
  do i=1,nframe-1
        ! get the change in position relative to frame 0 for all colloids
        dpos = allpos(i,:,:) - allpos(0,:,:)
        ! get the square of the change in position
        !   dpos = (particles, positions (x,y,z))
        !   dpos*dpos = (particles, squared positions (x*x, y*y, z*z)
        !   sum(,dim=2) sums the 2nd dimension (x*x+y*y+z*z)
        !   -> creates a 1D array of dpos2 for all particles
        dpos2 = sum(dpos*dpos, dim=2)
        ! sum all the squares of relative change in position, and take the mean
        msd_value = sum(dpos2)/ncolloid
        ! calculate the sample standard deviation (1/N-1) of the MSD
        msd_corrsmplvar = (sum(dpos2*dpos2)/dble(ncolloid) - msd_value**2) / dble(ncolloid-1)
        msd_smplstddev = dsqrt(msd_corrsmplvar)

        ! write data to the output file
        write(14,'(I5,x,f20.5,x,f20.5)') i*t1, msd_value, msd_smplstddev
  enddo

  ! close output file
  close(14)

end subroutine msd_calculation


!###########################################
! PCF TO BE ADDED


!###########################################
subroutine rdf_calc(Lbox,m_xy,nparticle,typeid,pos,bin_width,nlayers)
! calculate the radial distribution function (RDF) g(r), the density of (other)
! particles per unit volume relative to the position of a chosen particle (r1)
! see the sim-analysis.py script for more details

  implicit none

  !!! INPUTS
  ! the simulation box size (L_X, L_Y, L_Z), and the
  ! xy tilt factor (box deformation from shear in x-direction)
  real(8), intent(in) :: Lbox(0:2), m_xy
  ! total number of particles 
  integer, intent(in) :: nparticle
  ! the typeid of all particles
  integer, intent(in) :: typeid(0:nparticle-1)
  ! the xyz position of all particles
  real(8), intent(in) :: pos(0:nparticle-1,0:2)

  ! tell f2py about dependencies (required as commented out)
  !f2py depend(nparticle) typeid,pos 
 
  ! and the size of each bin/layer/slab in the search box
  real(8), intent(in) :: bin_width
  ! the number of bins/layers/slabs in the search box
  integer, intent(in) :: nlayers

  !!! OUTPUTS
  ! N/A (output files created)

  !!! INTERNAL VALUES
  ! 1/bin-size
  real(8) :: inv_bin_width 
  ! the center-center interparticle distance (dr) and 
  ! corrected center-center interparticle distance (r1),
  ! 1/box-size, and the modification for particles interacting
  ! across a boundary (img)
  real(8) :: dr(0:2), r1, inv_Lbox(0:2), img
  ! particle tags (i, j) and bin tag (k)
  integer :: i, j, k
  ! lists for counting the number of other particles at each bin
  ! hist_bin_cs = colloid-solvent; hist_bin_cc = colloid-colloid
  integer :: hist_bin_cs(0:nlayers-1), hist_bin_cc(0:nlayers-1)
  ! the number pi, the current and previous volumes
  real(8) :: pi, curr_volume, prev_volume
  ! the particle density calculated for a bin/layer
  ! (colloid-solvent and colloid-colloid)
  real(8) :: cs_density, cc_density 

  !!! SET INITIAL VALUES
  ! calculate 1/box-size
  inv_Lbox = 1.d0/Lbox
  ! calculate 1/bin-size
  inv_bin_width = 1.d0/bin_width
  ! calculate pi (pi = arctan(1) * 4)
  pi = 4.d0*datan(1.d0)
  ! initialize history bins at 0
  hist_bin_cs = 0
  hist_bin_cc = 0

  !!! count the number of particles in each bin/layer
  ! for N-1 particles
  do i=0,nparticle-2
        ! interacting eith the remaining unbinned particles
        do j=i+1,nparticle-1

                ! if COLLOID-SOLVENT interaction
                if(typeid(i) .ne. typeid(j)) then
                        ! calculate the center-center interparticle distance
                        dr = pos(i,:) - pos(j,:)
                        ! adjust dr for particles interacting across a boundary:
                        ! z-boundary
                        dr(2) = dr(2) - Lbox(2) * dnint(dr(2)*inv_Lbox(2))
                        ! y-boundary
                        img = Lbox(1) * dnint(dr(1)*inv_Lbox(1))
                        dr(1) = dr(1) - img
                        ! x-boundary
                        dr(0) = dr(0) - img * m_xy
                        dr(0) = dr(0) - Lbox(0) * dnint(dr(0)*inv_Lbox(0))
                        ! set the corrected center-center interparticle distance
                        r1 = dsqrt(sum(dr*dr))
                        ! find which bin/layer the interparticle distance is
                        ! located in (bins range 0:nlayers-1)
                        k = floor(r1*inv_bin_width)
                        ! if the interaction is inside the search box
                        if(k < nlayers) then
                                ! increase the number of particles interactions
                                ! in that bin by 2 (once for each particle)
                                hist_bin_cs(k) = hist_bin_cs(k) + 2
                        endif

                ! if a COLLOID-COLLOID interaction
                elseif((typeid(i) .eq. 1) .and. (typeid(j) .eq. 1)) then
                        ! calculate the center-center interparticle distance
                        dr = pos(i,:) - pos(j,:)
                        ! adjust dr for particles interacting across a boundary:
                        ! z-boundary
                        dr(2) = dr(2) - Lbox(2) * dnint(dr(2)*inv_Lbox(2))
                        ! y-boundary
                        img = Lbox(1) * dnint(dr(1)*inv_Lbox(1))
                        dr(1) = dr(1) - img
                        ! x-boundary
                        dr(0) = dr(0) - img * m_xy
                        dr(0) = dr(0) - Lbox(0) * dnint(dr(0)*inv_Lbox(0))
                        ! set the corrected center-center interparticle distance
                        r1 = dsqrt(sum(dr*dr))
                        ! find which bin/layer the interparticle distance is
                        ! located in (bins range 0:nlayers-1)
                        k = floor(r1*inv_bin_width)
                        ! if the interaction is inside the search box
                        if(k < nlayers) then
                                ! increase the number of particles interactions
                                ! in that bin by 2 (once for each particle)
                                hist_bin_cc(k) = hist_bin_cc(k) + 2
                        endif
                endif
        enddo
  enddo
                       
  ! create 2 files to record rdf data
  ! (one for colloid-solvent and one for colloid-colloid)
  open(unit=14,file='rdf-CS.txt',action='write')
  write(14,'(A,x,A)') 'position', 'density-per-unit-volume'
  open(unit=28,file='rdf-CC.txt',action='write')
  write(28, '(A,x,A)') 'position', 'density-per-unit-volume'

  ! for each bin/layer
  do i=0,nlayers-1
        ! calculate the volume for the current and previous bin/layer
        curr_volume = (4/3) * pi * ((i+1)*bin_width)**3
        prev_volume = (4/3) * pi * (i*bin_width)**3

        ! calculate the COLLOID-SOLVENT density per unit volume
        cs_density = hist_bin_cs(i) / (curr_volume - prev_volume)
        ! write COLLOID-SOLVENT data
        write(14,'(f20.5,x,f20.5)') (i+0.5d0)*bin_width, cs_density

        ! calculate the COLLOID-COLLOID density per unit volume
        cc_density = hist_bin_cc(i) / (curr_volume - prev_volume)
        ! write the COLLOID-COLLOID data
        write(28,'(f20.5,x,f20.5)') (i+0.5d0)*bin_width, cc_density

  enddo

  ! close files
  close(14)
  close(28)

end subroutine rdf_calc
