! Fortran module for analyzing a bimodal colloid sim GSD file
! NOTE: requires compilation with compile-module
! NOTE: is run from a matching sim-bimodal-analysis.py Python script
! NOTE: ASSUMES 2 colloid types and 1 solvent type
! (Rob Campbell)

! Calculates:
! * average coordination number (<Z>) for each frame
! * mean squared displacement (MSD) for each frame
! * g(r) of the final simulation frame (gelled state)
! * network edgelist (all bonded colloid pairs) for each frame
! * number density fluctuations ("index of dispersion") per frame

  ! NOTE: f2py only understands limited number of KIND parameters
  ! use real(kind=8) not real64

  ! NOTE: fortran indexes from 1 by default, switch to base 0 to match C/Python
  ! for indexing from 0, 0:ncolloid-1 is of length ncolloid, 0:2 is length 3

  ! NOTE: must tell f2py about any variable dependencies or you will get Value Errors


!###########################################
subroutine coordination_number(Lbox,ncolloid1,ncolloid2,ncolloids,typeid,R_C,pos,m_xy,cut_off,Zavg,Z_12avg,Z_21avg,Z_11avg,Z_22avg)
! calculates average coordination number (Z), AKA contact number

  implicit none

  !!! INPUTS
  ! the simulation box size (L_X, L_Y, L_Z)
  real(8), intent(in) :: Lbox(0:2)

  ! total number of colloids and number of each type
  integer, intent(in) :: ncolloid1, ncolloid2, ncolloids
  ! the typeid of all particles
  integer, intent(in) :: typeid(0:ncolloids-1)
  ! the radius of each colloid particle
  real(8), intent(in) :: R_C(0:ncolloids-1)
  ! xyz position of all colloid particles
  real(8), intent(in) :: pos(0:ncolloids-1,0:2)

  ! tell f2py about dependencies (required as commented out)
  !f2py depend(ncolloids) typeid,R_C,pos
 
  ! the cut-off distance for defining a bond, and
  ! the xy tilt factor (box deformation from shear in x-direction)
  real(8), intent(in) :: m_xy, cut_off

  !!! OUTPUTS
  ! type-independent average contact number
  real(8), intent(out) :: Zavg
  ! type-dependent average contact numbers
  real(8), intent(out) :: Z_12avg, Z_21avg, Z_11avg, Z_22avg

  !!! INTERNAL VALUES
  ! the center-center interparticle distance (dr), 
  ! the surface-surface interparticle distance (h_ij), 
  ! 1/box-size, and the modification for
  ! for particles interacting across a boundary (img)   
  real(8) :: dr(0:2), h_ij, inv_Lbox(0:2), img
  ! particle tags (i, j), type-independent coordination number (Z)
  integer :: i, j, Z
  ! colloid-type specific coordination numbers (Z)
  integer :: Z_12, Z_21, Z_11, Z_22

  !!! SET INITIAL VALUES
  ! calculate 1/box-size
  inv_Lbox = 1.d0/Lbox
  ! set the initial type-independent coordination number to zero
  Z = 0
  ! set the initial type1-type2 coordination number to zero
  Z_12 = 0
  ! set the initial type2-type1 coordination number to zero
  Z_21 = 0
  ! set the initial type1-type1 coordination number to zero
  Z_11 = 0
  ! set the initial type2-type2 coordination number to zero
  Z_22 = 0

  !!! calculate type-independent Z for each colloid particle
  ! for all colloid particles
  do i=0,(ncolloids-1)
        ! check it's interaction with each other colloid particle
        do j=0,(ncolloids-1)
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
                h_ij = dsqrt(sum(dr*dr)) - R_C(i) - R_C(j)
                ! if the 2nd particle is within the cut-off distance
                if(h_ij <= cut_off) then
                        ! increase the coordination number by 1
                        Z = Z + 1
                endif
        enddo
  enddo

  ! caluclate the average coordination number
  Z = Z - ncolloids
  Zavg = float(Z)/ncolloids

  !!! calculate type-dependent Zs for each colloid particle
  ! for all colloid particles
  do i=0,(ncolloids-1) 
        ! check it's interaction with each other colloid particle
        do j=0,(ncolloids-1) 
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
                h_ij = dsqrt(sum(dr*dr)) - R_C(i) - R_C(j)
                ! if the 2nd particle is within the cut-off distance
                if(h_ij <= cut_off) then
                        ! if i=type 1 interacting with type2 --> Z_12
                        if((typeid(i) .eq. 1) .AND. (typeid(i) .ne. typeid(j))) then
                                ! increase the Z_12 count by 1
                                Z_12 = Z_12 + 1
                        ! if i=type 2 interacting with type1 --> Z_21
                        else if((typeid(i) .eq. 2) .AND. (typeid(i) .ne. typeid(j))) then
                                ! increase the Z_12 count by 1
                                Z_21 = Z_21 + 1
                        ! if type1-type1 interaction --> Z_11
                        else if((typeid(i) .eq. typeid(j)) .AND. (typeid(i) .eq. 1)) then
                                ! increase the Z_11 count by 1
                                Z_11 = Z_11 + 1
                        ! if type2-type2 interaction --> Z22
                        else if((typeid(i) .eq. typeid(j)) .AND. (typeid(i) .eq. 2)) then
                                ! increase the Z_22 count by 1
                                Z_22 = Z_22 + 1
                        endif
                endif
        enddo
  enddo

  ! caluclate the average coordination numbers
  Z_12 = Z_12
  Z_12avg = float(Z_12)/ncolloid1 
  Z_21 = Z_21
  Z_21avg = float(Z_21)/ncolloid2 
  Z_11 = Z_11 - ncolloid1
  Z_11avg = float(Z_11)/ncolloid1
  Z_22 = Z_22 - ncolloid2 
  Z_22avg = float(Z_22)/ncolloid2

end subroutine coordination_number




!###########################################
subroutine msd_calculation(nframe,ncolloids,ncolloid1,ncolloid2,nsolvents,allpos,allpos1,allpos2,allsolvpos)
! calculate the mean-squared displacement of all colloid particles

  implicit none

  !!! INPUTS
  ! number of frames, total number of each type of particle
  integer, intent(in) :: nframe, ncolloids, ncolloid1, ncolloid2, nsolvents
  ! xyz positions of all colloids in all frames
  real(8), intent(in) :: allpos(0:nframe-1,0:ncolloids-1,0:2) 
  ! xyz positions of all type 1 colloids in all frames
  real(8), intent(in) :: allpos1(0:nframe-1,0:ncolloid1-1,0:2) 
  ! xyz positions of all type 2 colloids in all frames
  real(8), intent(in) :: allpos2(0:nframe-1,0:ncolloid2-1,0:2) 
  ! xyz positions of all solvent particles in all frames
  real(8), intent(in) :: allsolvpos(0:nframe-1,0:nsolvents-1,0:2)
  ! timestep converstion factor (t1 = period * dt_Integration)
  !real(8), intent(in) :: t1

  ! tell f2py about dependencies (required as commented out)
  !f2py depend(nframe,ncolloids) allpos
  !f2py depend(nframe,ncolloid1) allpos1
  !f2py depend(nframe,ncolloid2) allpos2
  !f2py depend(nframe,nsolvents) allsolvpos

  !!! OUTPUTS
  ! N/A (output file created)

  !!! INTERNAL VALUES
  ! the center-center change in position of each colloid
  ! (dpos), and it's value squared (dpos_2)
  real(8) :: dcollpos(0:ncolloids-1,0:2), dcollpos_2(0:ncolloids-1)
  ! the center-center change in position of each type 1 colloid
  ! (dpos1), and it's value squared (dpos1_2)
  real(8) :: dcollpos1(0:ncolloid1-1,0:2), dcollpos1_2(0:ncolloid1-1)
  ! the center-center change in position of each type 2 colloid
  ! (dpos2), and it's value squared (dpos2_2)
  real(8) :: dcollpos2(0:ncolloid2-1,0:2), dcollpos2_2(0:ncolloid2-1)
  ! the center-center change in position (dpos) of each
  ! solvent particle, and it's value squared (dpos2)
  real(8) :: dsolvpos(0:nsolvents-1,0:2), dsolvpos2(0:nsolvents-1)
  ! frame tag (i)
  integer :: i
  ! type 1 mean squared displacement, it's corrected sample variance,
  ! and it's (corrected) sample standard deviation
  real(8) :: coll_msd_value, coll_msd_corrsmplvar, coll_msd_smplstddev
  ! type 1 mean squared displacement, it's corrected sample variance,
  ! and it's (corrected) sample standard deviation
  real(8) :: coll_msd1_value, coll_msd1_corrsmplvar, coll_msd1_smplstddev
  ! type 2 mean squared displacement, it's corrected sample variance,
  ! and it's (corrected) sample standard deviation
  real(8) :: coll_msd2_value, coll_msd2_corrsmplvar, coll_msd2_smplstddev
  ! mean squared displacement of solvents, it's corrected sample variance,
  ! and it's (corrected) sample standard deviation
  real(8) :: solv_msd_value, solv_msd_corrsmplvar, solv_smplstddev

  !!! SET INITIAL VALUES
  ! N/A

  !!! create a file to record the MSD and standard deviation
  open(unit=14,file='msd.txt',action='write') 
  ! label the output data values
  write(14,'(A,x,A,x,A,x,A,x,A,x,A,x,A)') 'DPD-time','msd_colloid_1','msd1_smplstddev', &
'msd_colloid_2','msd2_smplstddev','msd_solvent','msd_solv_smplstddev', &
'msd_allcolloids', 'msd_all_smplstddev'

  !!! calculate the msd for each frame
  do i=1,nframe-1
        ! for all colloids
        ! get the change in position relative to frame 0 for all type 1 colloids
        dcollpos = allpos(i,:,:) - allpos(0,:,:)
        ! get the square of the change in position
        !   dcollpos = (particles, positions (x,y,z))
        !   dcollpos*dcollpos = (particles, squared positions (x*x, y*y, z*z)
        !   sum(,dim=2) sums the 2nd dimension (x*x+y*y+z*z)
        !   -> creates a 1D array of dcollpos1_2 for all particles
        dcollpos_2 = sum(dcollpos*dcollpos, dim=2) 
        ! sum all the squares of relative change in position, and take the mean
        coll_msd_value = sum(dcollpos_2)/ncolloids
        ! calculate the sample standard deviation (1/N-1) of the MSD
        coll_msd_corrsmplvar = (sum(dcollpos_2*dcollpos_2)/dble(ncolloids) - coll_msd_value**2) / dble(ncolloids-1) 
        coll_msd_smplstddev = dsqrt(coll_msd_corrsmplvar)

        ! for type 1
        ! get the change in position relative to frame 0 for all type 1 colloids
        dcollpos1 = allpos1(i,:,:) - allpos1(0,:,:)
        ! get the square of the change in position
        !   dcollpos1 = (particles, positions (x,y,z))
        !   dcollpos1*dcollpos1 = (particles, squared positions (x*x, y*y, z*z)
        !   sum(,dim=2) sums the 2nd dimension (x*x+y*y+z*z)
        !   -> creates a 1D array of dcollpos1_2 for all particles
        dcollpos1_2 = sum(dcollpos1*dcollpos1, dim=2) 
        ! sum all the squares of relative change in position, and take the mean
        coll_msd1_value = sum(dcollpos1_2)/ncolloid1
        ! calculate the sample standard deviation (1/N-1) of the MSD
        coll_msd1_corrsmplvar = (sum(dcollpos1_2*dcollpos1_2)/dble(ncolloid1) - coll_msd1_value**2) / dble(ncolloid1-1) 
        coll_msd1_smplstddev = dsqrt(coll_msd1_corrsmplvar)

        ! for type 2
        ! get the change in position relative to frame 0 for all type 2 colloids
        dcollpos2 = allpos2(i,:,:) - allpos2(0,:,:)
        ! get the square of the change in position
        !   dcollpos2 = (particles, positions (x,y,z))
        !   dcollpos2*dcollpos2 = (particles, squared positions (x*x, y*y, z*z)
        !   sum(,dim=2) sums the 2nd dimension (x*x+y*y+z*z)
        !   -> creates a 1D array of dcollpos2_2 for all particles
        dcollpos2_2 = sum(dcollpos2*dcollpos2, dim=2) 
        ! sum all the squares of relative change in position, and take the mean
        coll_msd2_value = sum(dcollpos2_2)/ncolloid2
        ! calculate the sample standard deviation (1/N-1) of the MSD
        coll_msd2_corrsmplvar = (sum(dcollpos2_2*dcollpos2_2)/dble(ncolloid2) - coll_msd2_value**2) / dble(ncolloid2-1)
        coll_msd2_smplstddev = dsqrt(coll_msd2_corrsmplvar)

        !! for SOLVENTS
        ! get the change in position relative to frame 0 for all solvent
        ! particles
        dsolvpos = allsolvpos(i,:,:) - allsolvpos(0,:,:)
        ! get the square of the change in position
        !   dsolvpos = (solvents, positions (x,y,z))
        !   dsolvpos*dsolvpos = (solvents, squared positions (x*x, y*y, z*z)
        !   sum(,dim=2) sums the 2nd matrix dimension (x*x+y*y+z*z)
        !   -> creates a 1D array of dsolvpos2 for all solvent particles
        dsolvpos2 = sum(dsolvpos*dsolvpos, dim=2)
        ! sum all the squares of relative change in position, and take the mean
        solv_msd_value = sum(dsolvpos2)/nsolvents
        ! calculate the sample standard deviation (1/N-1) of the MSD
        solv_msd_corrsmplvar = (sum(dsolvpos2*dsolvpos2)/dble(nsolvents) - solv_msd_value**2) / dble(nsolvents-1)
        solv_smplstddev = dsqrt(solv_msd_corrsmplvar)

        ! write data to the output file
        write(14,'(I5,x,f20.5,x,f20.5,x,f20.5,x,f20.5,x,f20.5,x,f20.5,x,f20.5,x,f20.5)') i, coll_msd1_value, &
coll_msd1_smplstddev, coll_msd2_value, coll_msd2_smplstddev, solv_msd_value, solv_smplstddev, coll_msd_value, &
coll_msd_smplstddev
  enddo

  ! close output file
  close(14)

end subroutine msd_calculation




!###########################################
subroutine gofr_calc(Lbox,m_xy,nparticle,ncolloid1,ncolloid2,nsolvents,typeid,pos,bin_width,nlayers,cs_gauss,c1s_gauss)
! calculate g(r) the probability density per unit volume 
! see the sim-analysis.py script for more details

  implicit none

  !!! INPUTS
  ! the simulation box size (L_X, L_Y, L_Z), and the
  ! xy tilt factor (box deformation from shear in x-direction)
  real(8), intent(in) :: Lbox(0:2), m_xy
  ! total number of particles 
  integer, intent(in) :: nparticle
  ! total number of particles in each group
  integer, intent(in) :: ncolloid1, ncolloid2, nsolvents
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
  ! the Guassian distribution of colloid-solvent
  real(8), intent(in) :: cs_gauss, c1s_gauss

  !!! OUTPUTS
  ! N/A (output files created)

  !!! INTERNAL VALUES
  !! For counting particles in each bin !!
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
  ! colloid type independent:
  ! hist_bin_cs = colloid-solvent; hist_bin_cc = colloid-colloid
  integer :: hist_bin_cs(0:nlayers-1), hist_bin_cc(0:nlayers-1)
  ! colloid type dependent: (c1 = colloid 1, c2 = colloid 2)
  integer :: hist_bin_c1s(0:nlayers-1), hist_bin_c2s(0:nlayers-1)
  integer :: hist_bin_c1c1(0:nlayers-1), hist_bin_c2c2(0:nlayers-1)
  integer :: hist_bin_c1c2(0:nlayers-1)

  !! For the final calculation !!
  ! the number pi, the current and previous volumes
  real(8) :: pi, curr_volume, prev_volume
  ! temporary 1/binvolume for each bin
  real(8) :: inv_binvolume
  ! total box volume
  real(8) :: Lbox_vol
  ! adjusted Gaussian number density of each bimodal pair
  !real(8) :: cs_gauss, c1s_gauss 
  real(8) :: c2s_gauss, c1c2_gauss 
  ! 1/numberdensity for each like-particle group
  real(8) :: inv_cc_nd
  real(8) :: inv_c1c1_nd, inv_c2c2_nd
  ! 1/number for each like-particle group
  real(8) :: inv_ncc
  real(8) :: inv_nc1c1, inv_nc2c2
  ! temporary normalized counts for like-particle bins
  real(8) :: norm_cc_num !norm_cs_num, norm_cc_num
  real(8) :: norm_c1c1_num, norm_c2c2_num
  ! the particle density calculated for a bin/layer
  real(8) :: cs_density, cc_density 
  real(8) :: c1s_density, c2s_density 
  real(8) :: c1c1_density, c2c2_density, c1c2_density 

  !!! SET INITIAL VALUES
  ! calculate 1/box-size
  inv_Lbox = 1.d0/Lbox
  ! calculate 1/bin-size
  inv_bin_width = 1.d0/bin_width
  ! calculate pi (pi = arctan(1) * 4)
  pi = 4.d0*datan(1.d0)
  ! calculate the total volume
  Lbox_vol = Lbox(0)*Lbox(1)*Lbox(2)
  ! calculate a Gaussian number density for each bimodal pair
  ! (n1 * n2 * 2) and adjust for spherical shells
  ! ...but for some reason cs and c1s must be calculated outside the module
  !cs_gauss = ((ncolloid1+ncolloid2)*nsolvents*2.d0/Lbox_vol)*((4.d0/3.d0)*pi) 
  !c1s_gauss = (ncolloid1*nsolvents*2.d0/Lbox_vol)*((4.d0/3.d0)*pi)
  c2s_gauss = (ncolloid2*nsolvents*2.d0/Lbox_vol)*((4.d0/3.d0)*pi)
  c1c2_gauss = (ncolloid1*ncolloid2*2.d0/Lbox_vol)*((4.d0/3.d0)*pi)
  ! calculate 1/numberdensity for each like-particle group
  inv_cc_nd = 1.d0 / ((ncolloid1 + ncolloid2) / Lbox_vol)
  inv_c1c1_nd = 1.d0 / (ncolloid1 / Lbox_vol)
  inv_c2c2_nd = 1.d0 / (ncolloid2 / Lbox_vol)
  ! calculate 1/number for each like-particle group
  inv_ncc = 1.d0 / (ncolloid1 + ncolloid2)
  inv_nc1c1 = 1.d0 / ncolloid1
  inv_nc2c2 = 1.d0 / ncolloid2
  ! initialize history bins at 0
  hist_bin_cs = 0 
  hist_bin_c1s = 0
  hist_bin_c2s = 0
  hist_bin_cc = 0
  hist_bin_c1c1 = 0 
  hist_bin_c2c2 = 0
  hist_bin_c1c2 = 0

  !!! count the number of particles in each bin/layer
  ! for N-1 particles
  do i=0,nparticle-2
        ! interacting eith the remaining unbinned particles
        do j=i+1,nparticle-1

                ! if COLLOID-SOLVENT interaction
                if((typeid(i) .ne. typeid(j)) .AND. ((typeid(i) .eq. 0) .OR. (typeid(j) .eq. 0))) then 
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
                                ! if type 1 colloid
                                if((typeid(i) .eq. 1) .OR. (typeid(j) .eq. 1)) then
                                hist_bin_c1s(k) = hist_bin_c1s(k) + 2
                                ! if type 2 colloid
                                else if ((typeid(i) .eq. 2) .OR. (typeid(j) .eq. 2)) then
                                hist_bin_c2s(k) = hist_bin_c2s(k) + 2 
                                endif 
                        endif

                ! if a COLLOID-COLLOID interaction
                elseif((typeid(i) .ne. 0) .AND. (typeid(j) .ne. 0)) then
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
                                ! if TYPE1-TYPE1
                                if((typeid(i) .eq. 1) .AND. (typeid(j) .eq. 1)) then
                                        hist_bin_c1c1(k) = hist_bin_c1c1(k) + 2
                                ! if TYPE2-TYPE2
                                else if((typeid(i) .eq. 2) .AND. (typeid(j) .eq. 2)) then 
                                        hist_bin_c2c2(k) = hist_bin_c2c2(k) + 2
                                ! if TYPE1-TYPE2
                                else if(typeid(i) .ne. typeid(j)) then
                                        hist_bin_c1c2(k) = hist_bin_c1c2(k) + 2
                                endif
                        endif
                endif
        enddo
  enddo
                       
  ! create file to record gofr data
  open(unit=15,file='gofr.txt',action='write')
  write(15,'(A,x,A,x,A,x,A,x,A,x,A,x,A,x,A)') 'position', 'cs_density', &
'c1s_density', 'c2s_density', 'cc_density', 'c1c1_density', 'c2c2_density', &
'c1c2_density' 

  ! for each bin/layer
  do i=0,nlayers-1
        ! calculate the volume for the current and previous bin/layer
        curr_volume = (4.d0/3.d0) * pi * ((i+1)*bin_width)**3
        prev_volume = (4.d0/3.d0) * pi * (i*bin_width)**3
        inv_binvolume = 1.d0 / (curr_volume - prev_volume)

        ! for comparing un-like particles: 
        ! normalize with the Gaussian distribution of a bimodal system
        ! COLLOID-SOLVENT 
        cs_density = hist_bin_cs(i)/cs_gauss/((i+1)**3-i**3)/bin_width**3
        ! COLLOID1-SOLVENT
        c1s_density = hist_bin_c1s(i)/c1s_gauss/((i+1)**3-i**3)/bin_width**3 
        ! COLLOID2-SOLVENT
        c2s_density = hist_bin_c2s(i)/c2s_gauss/((i+1)**3-i**3)/bin_width**3 


        ! COLLOID-COLLOID
        ! for comparing like-particles:
        ! normalize by the total number of particles considered
        norm_cc_num = hist_bin_cc(i) * inv_ncc 
        ! calculate the COLLOID-COLLOID density per unit volume
        ! normalized by the density of a uniform distribution
        cc_density = norm_cc_num * inv_binvolume * inv_cc_nd 

        ! COLLOID1-COLLOID1
        ! for comparing like-particles:
        ! normalize by the total number of particles considered
        norm_c1c1_num = hist_bin_c1c1(i) * inv_nc1c1 
        ! calculate the COLLOID1-COLLOID1 density per unit volume
        ! normalized by the density of a uniform distribution
        c1c1_density = norm_c1c1_num * inv_binvolume * inv_c1c1_nd

        ! COLLOID2-COLLOID2
        ! for comparing like-particles:
        ! normalize by the total number of particles considered
        norm_c2c2_num = hist_bin_c2c2(i) * inv_nc2c2 
        ! calculate the COLLOID2-COLLOID2 density per unit volume
        ! normalized by the density of a uniform distribution
        c2c2_density = norm_c2c2_num * inv_binvolume * inv_c2c2_nd

        ! COLLOID1-COLLOID2
        ! for comparing un-like particles: 
        ! normalize with the Gaussian distribution of a bimodal system
        c1c2_density = hist_bin_c1c2(i)/c1c2_gauss/((i+1)**3-i**3)/bin_width**3  

        ! write the data
        write(15,'(f20.5,x,f20.5,x,f20.5,x,f20.5,x,f20.5,x,f20.5,x,f20.5,x,f20.5)') (i+0.5d0)*bin_width, & 
cs_density, c1s_density, c2s_density, cc_density, c1c1_density, c2c2_density, c1c2_density

  enddo

  ! close files
  close(15)

end subroutine gofr_calc





!###########################################
subroutine edgelist_calc(nframes,ncolloids,radii,allpos,lbox,rcut,edge_output)
implicit none

! INPUTS
integer, intent(in) :: nframes                              ! number of colloids
integer, intent(in) :: ncolloids                            ! number of colloids
real(8), intent(in) :: radii(0:ncolloids-1)                  ! colloid radii
real(8), intent(in) :: allpos(0:nframes-1,0:ncolloids-1,0:2) ! colloid positions
!f2py depend(ncolloids) radii
!f2py depend(nframes,ncolloids) allpos
real(8), intent(in) :: lbox(0:2)                     ! simulation box size
real(8), intent(in) :: rcut                        ! surface-surface cut-off for a bond
character*20, intent(in) :: edge_output            ! output filepath for edgelist

! OUTPUTS
!N/A

! INTERNAL VARIABLES
real(8) :: rij(0:2), hij                   ! center-center (rij) and surface-surface (hij) interaction distance
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
  do i=0,ncolloids-1
    do j=i+1,ncolloids
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





!###########################################
subroutine ndfluc(nframe,nparticle,Lbox,typeid,allpos,binwidth_x,binwidth_y,binwidth_z,nlayers_x,nlayers_y,nlayers_z)
! calculates number density fluctuation, AKA "index of dispersion"
! see the sim-analysis.py script for more details

  implicit none

  !!! INPUTS
  ! number of frames and number of particles
  integer, intent(in) :: nframe, nparticle
  ! the simulation box size (L_X, L_Y, L_Z)
  real(8), intent(in) :: Lbox(0:2)
  ! particle type IDs
  integer, intent(in) :: typeid(0:nparticle-1)
  ! the dimesions of each bin
  real(8), intent(in) :: binwidth_x, binwidth_y, binwidth_z  
  ! the number of layers in each dimensions
  integer, intent(in) :: nlayers_x, nlayers_y, nlayers_z
  ! xyz position of each particle
  real(8), intent(in) :: allpos(0:nframe-1,0:nparticle-1,0:2)
  ! timestep converstion factor (t1 = period * dt_Integration)
  !integer, intent(in) :: t1

  ! tell f2py about dependencies (required as commented out)
  !f2py depend(nparticle) typeid
  !f2py depend(nframe,nparticle) allpos

  !!! OUTPUTS
  ! N/A (output file created)

  !!! INTERNAL VALUES
  ! the volume of each bin
  real(8) :: binvol
  ! inverse bin widths (1/binwidth)
  real(8) :: inv_binwidth_x, inv_binwidth_y, inv_binwidth_z
  ! the number of bins
  integer :: nbins
  ! indices for looping through frames and particles
  integer :: i, j
  ! tags for bin position
  integer :: x_bin, y_bin, z_bin
  ! history bins for recording particle counts
  real(8) :: hist_bin_c(0:nlayers_x-1,0:nlayers_y-1,0:nlayers_z-1) 
  real(8) :: hist_bin_c1(0:nlayers_x-1,0:nlayers_y-1,0:nlayers_z-1) 
  real(8) :: hist_bin_c2(0:nlayers_x-1,0:nlayers_y-1,0:nlayers_z-1) 
  real(8) :: hist_bin_s(0:nlayers_x-1,0:nlayers_y-1,0:nlayers_z-1) 
  ! the number of binned particles
  real(8) :: sum_histbin_c, sum_histbin_c1, sum_histbin_c2, sum_histbin_s
  ! the means for each particle type
  real(8) :: mean_c, mean_c1, mean_c2, mean_s
  ! the squares of each particle, and their averages N^2 and <N^2>
  real(8) :: sq_histbin_c, sq_histbin_c1, sq_histbin_c2, sq_histbin_s
  real(8) :: NC_sq, NC1_sq, NC2_sq, NS_sq
  ! the variance of each particle type
  real(8) :: variance_sq_c, variance_sq_c1, variance_sq_c2, variance_sq_s
  ! the number density fluctuations for each particle type
  real(8) :: NDF_C, NDF_C1, NDF_C2, NDF_S

  !!! SET INITIAL VALUES
  ! calculate bin volume
  binvol = binwidth_x * binwidth_y * binwidth_z
  ! calculate 1/binwidths
  inv_binwidth_x = 1.d0/binwidth_x
  inv_binwidth_y = 1.d0/binwidth_y
  inv_binwidth_z = 1.d0/binwidth_z
  ! calculate number of bins
  nbins = nlayers_x * nlayers_y * nlayers_z


  !!! create a file to record the NDF values
  open(unit=14,file='ndfluc.txt',action='write')
  ! label the output data values
  write(14,'(A,x,A,x,A,x,A,x,A)') 'simframe','ndfluc_C','ndfluc_C1','ndfluc_C2', 'ndfluc_S'


  !!! calculate the NDF for each frame
  do i=0,nframe-1
        ! assign particles to bins
        do j=0,nparticle-1
                ! scale bin position for a box size -L/2 to L/2 (not 0 to L)
                x_bin = dnint(allpos(i,j,0)+0.5*Lbox(0))*inv_binwidth_x
                y_bin = dnint(allpos(i,j,1)+0.5*Lbox(1))*inv_binwidth_y
                z_bin = dnint(allpos(i,j,2)+0.5*Lbox(2))*inv_binwidth_z
                ! for particles outside the layers, assign them to the last bin
                if(x_bin >= nlayers_x) then
                        x_bin = nlayers_x-1
                else if(y_bin >= nlayers_y) then
                        y_bin = nlayers_y-1
                else if(z_bin >= nlayers_z) then
                        z_bin = nlayers_z-1
                endif
                ! assign to the history bin based on particle type
                ! COLLOID1
                if(typeid(j) .eq. 1) then
                        hist_bin_c1(x_bin,y_bin,z_bin) = hist_bin_c1(x_bin,y_bin,z_bin) + 1
                        hist_bin_c(x_bin,y_bin,z_bin) = hist_bin_c(x_bin,y_bin,z_bin) + 1
                ! COLLOID2
                else if(typeid(i) .eq. 2) then
                        hist_bin_c2(x_bin,y_bin,z_bin) = hist_bin_c2(x_bin,y_bin,z_bin) + 1
                        hist_bin_c(x_bin,y_bin,z_bin) = hist_bin_c(x_bin,y_bin,z_bin) + 1
                ! SOLVENT
                else if(typeid(i) .eq. 0) then
                        hist_bin_s(x_bin,y_bin,z_bin) = hist_bin_s(x_bin,y_bin,z_bin) + 1
                endif
        enddo
              
        ! calculate the mean, AKA the average number density per bin, <N>
        sum_histbin_c = sum(hist_bin_c/binvol) ! N_C
        mean_c = sum_histbin_c / nbins !<N_C>

        sum_histbin_c1 = sum(hist_bin_c1/binvol) ! N_C1
        mean_c1 = sum_histbin_c1 / nbins !<N_C1>

        sum_histbin_c2 = sum(hist_bin_c2/binvol) ! N_C2
        mean_c2 = sum_histbin_c2 / nbins !<N_C2>

        sum_histbin_s = sum(hist_bin_s/binvol) ! N_S
        mean_s = sum_histbin_s / nbins !<N_S>

        ! calculate the variance squared (<N^2> - <N>^2)
        sq_histbin_c = sum((hist_bin_c/binvol)**2) ! N_C^2
        NC_sq = sq_histbin_c / nbins ! <N_C^2>
        variance_sq_c = NC_sq - (mean_c*mean_c) ! <N_C^2> - <N_C>^2
 
        sq_histbin_c1 = sum((hist_bin_c1/binvol)**2) ! N_C1^2
        NC1_sq = sq_histbin_c1 / nbins ! <N_C1^2>
        variance_sq_c1 = NC1_sq - (mean_c1*mean_c1) ! <N_C1^2> - <N_C1>^2

        sq_histbin_c2 = sum((hist_bin_c2/binvol)**2) ! N_C2^2
        NC2_sq = sq_histbin_c2 / nbins ! <N_C2^2>
        variance_sq_c2 = NC2_sq - (mean_c2*mean_c2) ! <N_C2^2> - <N_C2>^2

        sq_histbin_s = sum((hist_bin_s/binvol)**2) ! N_S^2
        NS_sq = sq_histbin_s / nbins ! <N_S^2>
        variance_sq_s = NS_sq - (mean_s*mean_s) ! <N_S^2> - <N_S>^2

        ! calculate the number density fluctuation
        NDF_C = variance_sq_c / mean_c
        NDF_C1 = variance_sq_c1 / mean_c1
        NDF_C2 = variance_sq_c2 / mean_c2
        NDF_S = variance_sq_s / mean_s

        ! write data to the output file
        write(14,'(I5,x,f20.5,x,f20.5,x,f20.5,x,f20.5)') (i+1), NDF_C, NDF_C1, NDF_C2, NDF_S 

  enddo

  ! close output file
  close(14)

end subroutine ndfluc                 
