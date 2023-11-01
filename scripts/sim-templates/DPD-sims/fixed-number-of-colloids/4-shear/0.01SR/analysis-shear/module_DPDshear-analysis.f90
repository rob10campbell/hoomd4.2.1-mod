! Fortran module for analyzing a colloid sim.
! NOTE: requires compilation with compile-module
! NOTE: is run from a matching sim-analysis-shear.py Python script
! (Rob Campbell)

! Calculates:
! * corrected (kinetic) pressure and temperature

  ! NOTE: f2py only understands limited number of KIND parameters
  ! use real(kind=8) not real64

  ! NOTE: fortran indexes from 1 by default, switch to base 0 to match C/Python
  ! for indexing from 0, 0:nframe-1 is of length nframe, 0:2 is length 3


!##############################################
subroutine corrected_temperature(nlayerY,nparticle,typeid,Lbox,mass,pos,vel,pxx,pxy,pxz,kT)
! see description in sim-analysis.py

  implicit none

  !!! INPUTS
  ! total number of layers/bins (y-direction), total number of particles
  integer, intent(in) :: nlayerY, nparticle  

  ! tell f2py about dependencies (required as commented out)
  !f2py depend(nparticle) typeid,mass,pos,vel

  ! the typeid of all particles
  integer, intent(in) :: typeid(0:nparticle-1)
  ! size of the simulation box (L_X, L_Y, L_Z)
  real(8), intent(in) :: Lbox(0:2)
  ! the mass of all particles
  real(8), intent(in) :: mass(0:nparticle-1)
  ! xyz position of all particles
  real(8), intent(in) :: pos(0:nparticle-1,0:2) 
  ! velocity vector (xyz) for all particles 
  real(8), intent(in) :: vel(0:nparticle-1,0:2)

  !!! OUTPUTS
  ! the xx, xy, and xz compoenents of the pressure tensor
  ! and the kinetic temperature (kT)
  real(8), intent(out) :: pxx, pxy, pxz, kT

  !!! INTERNAL VALUES
  ! bin/layer size, 1/bin-size, x-component of velocity
  real(8) :: bin_width, inv_bin_width, vel_x
  ! a list for storing particle x-velocities in each bin/layer
  real(8) :: velx_bin(0:nlayerY-1)
  ! a list for storing the number of particles in each bin/layer
  integer :: hist_bin(0:nlayerY-1)
  ! particle tag (i) and bin tag (k)
  integer :: i, k
  ! simulation box volume
  real(8) :: Vbox
  
  !!! SET INITIAL VALUES
  pxx=0.d0; pxy=0.d0; pxz=0.d0; kT=0.d0
  hist_bin = 0
  velx_bin = 0.d0
  bin_width = Lbox(1)/nlayerY
  inv_bin_width = 1.d0/bin_width
  Vbox = 0.d0

  !!! first bin SOLVENTS only
  ! for all particles
  do i=0,nparticle-1
    ! if it is a solvent particle
    if(typeid(i) == 0) then
      ! scale the y-position by 1/2 the box size and 
      ! use that to find it's bin/layer (bin range 1:nlayerY-1)
      k = floor((pos(i,1)+0.5d0*Lbox(1))*inv_bin_width)
      ! wrap outlying particles into the last bin
      if(k >= nlayerY) then
        k = nlayerY - 1
      endif
      ! increase the particle count in the assigned bin by 1
      hist_bin(k) = hist_bin(k) + 1
      ! add the particles v_x to the matching x-velocity bin
      velx_bin(k) = velx_bin(k) + vel(i,0)
    endif
  enddo

  ! check all layers/bins and update to avoid discontinuities
  do i=0,nlayerY-1
    ! if there are no particles in that layer
    if(hist_bin(i) == 0) then
      ! add one particle to avoid discontinuities
      hist_bin(i) = 1
    endif
  enddo

  ! normalize each x-veloctiy bin by its number of particles
  ! (i.e. average SOLVENT velocity per layer)
  velx_bin = velx_bin/hist_bin


  !!! then for ALL PARTICLES
  ! correct the pressure (normalize x-velocity by bin/layer's solvent velocity)
  do i=0,nparticle-1
    ! scale the y-position by 1/2 the box size and 
    ! use that to find it's bin/layer (bin range 1:nlayerY-1)
    k = floor((pos(i,1)+0.5d0*Lbox(1))*inv_bin_width)
    ! wrap outlying particles into the last bin
    if(k >= nlayerY) then
      k = nlayerY - 1
    endif

    ! find the corrected x-component of particle velocity
    ! (relative to the bin's (average) solvent velocity)
    vel_x = vel(i,0) - velx_bin(k)

    ! use the corrected velocity to update the pressure
    ! (removes effect of solvent velocity profile; gives
    !    the pressure from colloid motion ONLY)
    pxx = pxx + mass(i) * (vel_x * vel_x - vel(i,0)*vel(i,0))
    pxy = pxy + mass(i) * (vel_x * vel(i,1) - vel(i,0)*vel(i,1))
    pxz = pxz + mass(i) * (vel_x * vel(i,2) - vel(i,0)*vel(i,2))

  enddo

  ! use the corrected pressure to calculate corrected kinetic temperature
  ! per dimension per particle
  kT = pxx/3.d0/nparticle

  ! find simulation box volume
  Vbox = Lbox(0)*Lbox(1)*Lbox(2)
  ! normalize corrected pressure by the simulation box volume
  pxx = pxx/Vbox
  pxy = pxy/Vbox
  pxz = pxz/Vbox

end subroutine corrected_temperature
