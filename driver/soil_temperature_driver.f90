!
! compile: 
!

program soil_temperature

use output
use routines, only: tsnosoi, thermoprop, noahmp_parameters, noahmp_options

  implicit none

!---------------------------------------------------------------------
!  inputs start
!---------------------------------------------------------------------

  real          :: dt
  integer       :: maxtime
  character*256 :: output_filename
  integer       :: output_freq
  integer       :: solution_method
  real          :: temperature_mean
  real          :: temperature_amplitude_daily
  real          :: temperature_amplitude_annual
  integer       :: isltyp
  integer       :: nsoil
  integer       :: nsnow
  real          :: tbot
  integer       :: structure_option
  real          :: soil_depth
  integer       :: bottom_temperature_option
  integer       :: soil_boundary_condition_option

  real, allocatable, dimension(:) :: zsoil   ! depth of layer-bottom from soil surface
  real, allocatable, dimension(:) :: dzsnso  ! snow/soil layer thickness [m]
  real, allocatable, dimension(:) :: sice    ! soil ice content [m3/m3]
  real, allocatable, dimension(:) :: sh2o    ! soil liquid water content [m3/m3]
  real, allocatable, dimension(:) :: stc     ! soil temperature [K]
  logical :: initial_theory                  ! initial with theoretical temperature
  logical :: initial_uniform                 ! initial all levels the same
  real    :: initial_sh2o_value              ! constant sh2o value
  real    :: initial_sice_value              ! constant sice value
  real    :: initial_stc_value               ! constant stc value

  !--------------------!
  !  soil parameters   !
  !--------------------!

  real, dimension(12) ::        qtz  ! quartz content
  real, dimension(12) ::     maxsmc  ! porosity
  real                :: csoil_data  ! soil heat capacity
  real                ::  zbot_data  ! deep soil temperature depth
     
  namelist / timing          / dt,maxtime,output_filename,solution_method,output_freq
  namelist / forcing         / temperature_mean, temperature_amplitude_daily, &
                               temperature_amplitude_annual
  namelist / structure       / isltyp,nsoil,nsnow,structure_option,soil_depth,tbot,&
                               initial_theory
  namelist / fixed_initial   / zsoil,dzsnso,sice,sh2o,stc
  namelist / uniform_initial / initial_uniform,initial_sh2o_value,&
                               initial_sice_value,initial_stc_value
  namelist / soil_parameters / qtz,maxsmc,csoil_data,zbot_data
  namelist / options         / bottom_temperature_option, &
                               soil_boundary_condition_option
 
!---------------------------------------------------------------------
!  inputs end
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  additional variables passed to thermoprop and tsnosoi
!---------------------------------------------------------------------

  integer                         :: isnow
  integer                         :: ice     ! not used
  real                            :: snowh
  real                            :: ssoil
  real                            :: tg
  real                            :: tg_previous  ! previous timestep tg
  real                            :: ur      ! not used
  real                            :: lat     ! not used
  real                            :: z0m     ! not used
  real                            :: zlvl    ! not used
  real                            :: sag     ! not used
  character(len=128)              :: errmsg
  integer                         :: errflg

  real, allocatable, dimension(:) :: smc         !total soil water content [m3/m3]
  real, allocatable, dimension(:) :: zsnso       !depth from snow surface (negative) [m]
  real, allocatable, dimension(:) :: df          !thermal conductivity [W/m/K]
  real, allocatable, dimension(:) :: hcpct       !heat capacity [J/m3/K]
  real, allocatable, dimension(:) :: snice       !snow layer ice
  real, allocatable, dimension(:) :: snliq       !snow layer liquid
  real, allocatable, dimension(:) :: snicev      !not needed
  real, allocatable, dimension(:) :: snliqv      !not needed
  real, allocatable, dimension(:) :: epore       !not needed
  real, allocatable, dimension(:) :: fact        !not needed

  integer                         :: ist   = 1   !surface type
  integer                         :: iloc  = 1   !grid index
  integer                         :: jloc  = 1   !grid index
  integer                         :: vegtyp = 14 !vegetation type

!---------------------------------------------------------------------
!  local variables
!---------------------------------------------------------------------

  integer :: itime, iz          ! some loop counters
  integer :: simulation_time    ! total time into simulation
  integer :: ntime      = 0     ! number of timesteps to run
  real, allocatable, dimension(:) :: theoretical_temperature   ! estimate T from diffusion solution
  real, allocatable, dimension(:) :: depth_midpoint            ! calculation node depth
  real, allocatable, dimension(:) :: depth_interface           ! calculation node depth
  real    :: damp_depth_daily
  real    :: damp_depth_annual
  real    :: period_daily  = 3600.0 * 24
  real    :: period_annual = 3600.0*24*365
  real, parameter :: pi = 3.14159265
  real, parameter :: fillvalue = 1.d30
  real    :: storage_before, storage_after, bottom_flux, energy_balance

!---------------------------------------------------------------------
!  parameters
!---------------------------------------------------------------------

  type (noahmp_parameters) :: parameters
  
!---------------------------------------------------------------------
!  end declarations
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  read input file
!---------------------------------------------------------------------

  open(30, file="namelist.input", form="formatted")
   read(30, timing)
   read(30, forcing)
   read(30, structure)
   read(30, uniform_initial)
   read(30, soil_parameters)
   read(30, options)
  close(30)

!---------------------------------------------------------------------
!  allocate for dynamic levels
!---------------------------------------------------------------------

  allocate (zsoil (       1:nsoil))
  allocate (dzsnso(-nsnow+1:nsoil))
  allocate (zsnso (-nsnow+1:nsoil))
  allocate (sice  (       1:nsoil))
  allocate (sh2o  (       1:nsoil))
  allocate (smc   (       1:nsoil))
  allocate (stc   (       1:nsoil))
  allocate (df    (-nsnow+1:nsoil))
  allocate (hcpct (-nsnow+1:nsoil))
  allocate (snice (-nsnow+1:0    ))
  allocate (snliq (-nsnow+1:0    ))
  allocate (snicev(-nsnow+1:0    ))
  allocate (snliqv(-nsnow+1:0    ))
  allocate (epore (-nsnow+1:0    ))
  allocate (fact  (-nsnow+1:nsoil))
  
  allocate (theoretical_temperature (1:nsoil))   ! estimate T from diffusion solution [K]
  allocate (depth_midpoint  (1:nsoil))           ! midpoint depth (negative) [m]
  allocate (depth_interface (1:nsoil))           ! interface depth (negative) [m]

  allocate (parameters%quartz(nsoil))
  allocate (parameters%smcmax(nsoil))

!---------------------------------------------------------------------
!  transfer parameters
!---------------------------------------------------------------------

  parameters%smcmax = maxsmc(isltyp)
  parameters%quartz = qtz(isltyp)

  parameters%csoil  = csoil_data
  parameters%zbot   =  zbot_data
  
  parameters%urban_flag = .false.
  
!---------------------------------------------------------------------
!  transfer options
!---------------------------------------------------------------------

  call noahmp_options(bottom_temperature_option, soil_boundary_condition_option)

!---------------------------------------------------------------------
!  read input file, part 2: initialize
!---------------------------------------------------------------------

  if(structure_option == 1) then       ! user-defined levels
    open(30, file="namelist.input", form="formatted")
     read(30, fixed_initial)
    close(30)
    if(initial_theory) print*, "prescribed temperature ICs will be replaced with theoretical"
  elseif(structure_option == 2 .or. structure_option == 3) then  ! fixed levels
    dzsnso = soil_depth / nsoil    ! soil layer thickness, positive
    do iz = 1, nsoil
      zsoil(iz) = -1. * sum(dzsnso(1:iz))  ! depth to interface, negative
    end do
    if(.not.initial_uniform) &
      stop "structure_option > 1 must have initial_uniform == .true."
  end if
  
  if(initial_uniform) then
    sh2o = initial_sh2o_value
    sice = initial_sice_value
    stc  = initial_stc_value
    if(initial_theory) print*, "uniform temperature ICs will be replaced with theoretical"
  end if

!---------------------------------------------------------------------
! initialize any other values
!---------------------------------------------------------------------

  smc       = sh2o + sice  ! initial volumetric soil water

  do iz = 1, nsoil
    if(smc(iz) > parameters%smcmax(iz)) then  ! adjust if oversaturated
      print*, "adjusting oversaturation in level: ",iz
      print*, "  old smc: ",smc(iz)
      print*, "  porosity: ",parameters%smcmax(iz)
      sh2o(iz) = sh2o(iz) * parameters%smcmax(iz) / smc(iz)
      sice(iz) = sice(iz) * parameters%smcmax(iz) / smc(iz)
      smc(iz)  = sh2o(iz) + sice(iz)
      print*, "  new smc: ",smc(iz)
    end if
  end do

  tg        = temperature_mean  ! surface temperature
  isnow     = 0                 !
  snowh     = 0.0               !
  ssoil     = fillvalue         !
  ur        = huge(0.0)         ! 
  lat       = huge(0.0)         ! 
  z0m       = huge(0.0)         ! 
  zlvl      = huge(0.0)         ! 
  sag       = huge(0.0)         !
  ice       = huge(0)           !
  df        = fillvalue
  hcpct     = fillvalue
  theoretical_temperature = fillvalue
  simulation_time         = 0.0
  bottom_flux             = fillvalue
  energy_balance          = fillvalue
  storage_before          = fillvalue
  storage_after           = fillvalue
  
  ntime      =  nint(maxtime * 3600.0 / dt)

  zsnso          = 0.0
  zsnso(1:nsoil) = zsoil          ! interface depths

  depth_interface = zsnso

  depth_midpoint(1) = 0.5 * zsoil(1)
  do iz = 2, nsoil
    depth_midpoint(iz) = 0.5*(zsoil(iz-1) + zsoil(iz))
  end do

  if(initial_theory) then
   ! get df and hcpct for theoretical temperature
    call thermoprop (parameters  ,nsoil   ,nsnow   ,isnow   ,ist     ,dzsnso  , & !in
                         dt      ,snowh   ,snice   ,snliq   ,                   & !in
                         smc     ,sh2o    ,tg      ,stc     ,ur      ,          & !in
                         lat     ,z0m     ,zlvl    ,vegtyp  ,                   & !in
                         df      ,hcpct   ,snicev  ,snliqv  ,epore   ,          & !out
                         fact    )                                                !out

    damp_depth_daily  = sqrt(period_daily*df(1)/hcpct(1)/pi)
    damp_depth_annual = sqrt(period_annual*df(1)/hcpct(1)/pi)

    if(structure_option == 1 .or. structure_option == 2) then
      stc = temperature_mean + &
            temperature_amplitude_daily  * exp(depth_midpoint/damp_depth_daily)  * sin(2*pi/period_daily*simulation_time  + depth_midpoint/damp_depth_daily) + &
            temperature_amplitude_annual * exp(depth_midpoint/damp_depth_annual) * sin(2*pi/period_annual*simulation_time + depth_midpoint/damp_depth_annual)
    elseif(structure_option == 3) then
      stc = temperature_mean + &
            temperature_amplitude_daily  * exp(depth_interface/damp_depth_daily)  * sin(2*pi/period_daily*simulation_time  + depth_interface/damp_depth_daily) + &
            temperature_amplitude_annual * exp(depth_interface/damp_depth_annual) * sin(2*pi/period_annual*simulation_time + depth_interface/damp_depth_annual)
    end if

  end if

!---------------------------------------------------------------------
! create output file and add initial values
!---------------------------------------------------------------------

  if(mod(ntime,output_freq) /= 0) stop "ntime must be divisible by output_freq"
  call initialize_output(output_filename, ntime+1, output_freq, nsoil)
  call add_to_output(0,output_freq,nsoil,tg,stc(1:nsoil),df(1:nsoil),hcpct(1:nsoil),ssoil, &
                       theoretical_temperature,energy_balance, bottom_flux)

!---------------------------------------------------------------------
! start the time loop
!---------------------------------------------------------------------

  do itime = 1, ntime
  
    simulation_time = itime * dt
    
  !---------------------------------------------------------------------
  ! calculate the surface temperature
  !---------------------------------------------------------------------

    tg = temperature_mean + &
         temperature_amplitude_daily  * sin(2*pi/period_daily*simulation_time) + &
         temperature_amplitude_annual * sin(2*pi/period_annual*simulation_time)
   
  !---------------------------------------------------------------------
  ! call the routines
  !---------------------------------------------------------------------

    call thermoprop (parameters  ,nsoil   ,nsnow   ,isnow   ,ist     ,dzsnso  , & !in
                         dt      ,snowh   ,snice   ,snliq   ,                   & !in
                         smc     ,sh2o    ,tg      ,stc     ,ur      ,          & !in
                         lat     ,z0m     ,zlvl    ,vegtyp  ,                   & !in
                         df      ,hcpct   ,snicev  ,snliqv  ,epore   ,          & !out
                         fact    )                                                !out

    storage_before = sum( hcpct(1:nsoil) * stc(1:nsoil) * dzsnso(1:nsoil) )
   
    damp_depth_daily  = sqrt(period_daily*df(1)/hcpct(1)/pi)
    damp_depth_annual = sqrt(period_annual*df(1)/hcpct(1)/pi)
    
    if(itime == 1) print *, "damping depth daily: ",damp_depth_daily
    if(itime == 1) print *, "damping depth annual: ",damp_depth_annual

    if(structure_option == 1 .or. structure_option == 2) then
      theoretical_temperature = temperature_mean + &
           temperature_amplitude_daily  * exp(depth_midpoint/damp_depth_daily)  * sin(2*pi/period_daily*simulation_time  + depth_midpoint/damp_depth_daily) + &
           temperature_amplitude_annual * exp(depth_midpoint/damp_depth_annual) * sin(2*pi/period_annual*simulation_time + depth_midpoint/damp_depth_annual)
    elseif(structure_option == 3) then
      theoretical_temperature = temperature_mean + &
           temperature_amplitude_daily  * exp(depth_interface/damp_depth_daily)  * sin(2*pi/period_daily*simulation_time  + depth_interface/damp_depth_daily) + &
           temperature_amplitude_annual * exp(depth_interface/damp_depth_annual) * sin(2*pi/period_annual*simulation_time + depth_interface/damp_depth_annual)
    end if

    if(solution_method == 0 .and. (structure_option == 1 .or. structure_option == 2)) then  ! use noahmp method

     ! top and bottom fluxes defined from previous timestep temperatures 

      ssoil = df(isnow+1)/(0.5*dzsnso(isnow+1)) * (tg - stc(isnow+1))
      if(bottom_temperature_option == 1) then
       bottom_flux = 0.0
      elseif(bottom_temperature_option == 2) then
       bottom_flux = df(nsoil)/(depth_midpoint(nsoil)-parameters%zbot) * (stc(nsoil)-tbot)
      end if
          
      call tsnosoi (parameters,ice     ,nsoil   ,nsnow   ,isnow   ,ist     , & !in
                  tbot      ,zsnso   ,ssoil   ,df      ,hcpct   ,          & !in
                  sag       ,dt      ,snowh   ,dzsnso  ,                   & !in
                  tg        ,iloc    ,jloc    ,                            & !in
                  stc       ,errmsg  ,errflg     )                           !inout

    elseif(solution_method == 1 .and. (structure_option == 1 .or. structure_option == 2)) then   ! use diffusion_implicit_midpoint subroutine

      call diffusion_implicit_midpoint(nsoil, zsnso, dt, tbot, parameters%zbot, &
                                         df, hcpct, tg, stc, bottom_temperature_option)
    
     ! top and bottom fluxes defined from updated temperatures 

      ssoil = df(isnow+1)/(0.5*dzsnso(isnow+1)) * (tg - stc(isnow+1))

      if(bottom_temperature_option == 1) then
        bottom_flux = 0.0
      elseif(bottom_temperature_option == 2) then
        bottom_flux = df(nsoil)/(depth_midpoint(nsoil)-parameters%zbot) * (stc(nsoil)-tbot)
      end if

    elseif(solution_method == 2 .and. structure_option == 3) then ! use diffusion_implicit_interface

      if(bottom_temperature_option /= 1) stop "only zero flux for the structure option 3"

      call diffusion_implicit_interface(nsoil, depth_midpoint, depth_interface, dt, df, hcpct, tg, stc)
    
      ssoil = -1.0 * df(1)/depth_interface(1) * (tg - stc(1))

      bottom_flux = 0.0

    elseif(solution_method == 3 .and. structure_option == 3) then ! use diffusion_implicit_interface

      if(bottom_temperature_option /= 1) stop "only zero flux for the structure option 3"

      call diffusion_cn_interface(nsoil, depth_midpoint, depth_interface, dt, df, hcpct, tg, tg_previous, stc)
    
      ssoil = -1.0 * df(1)/depth_interface(1) * (tg - stc(1))

      bottom_flux = 0.0

    else
    
      print *, "structure_option: ",structure_option
      print *, "solution_method: ",solution_method
      stop "no valid solution-structure combination"

    end if
    
    tg_previous = tg
    
    storage_after = sum( hcpct(1:nsoil) * stc(1:nsoil) * dzsnso(1:nsoil) )
    
    energy_balance = storage_after - storage_before - ssoil * dt + bottom_flux * dt

  !---------------------------------------------------------------------
  ! accumulate some fields and error checks
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  ! add to output file
  !---------------------------------------------------------------------

    if(mod(itime,output_freq) == 0) then

      call add_to_output(itime,output_freq,nsoil,tg,stc(1:nsoil),df(1:nsoil),hcpct(1:nsoil),ssoil, &
                         theoretical_temperature,energy_balance, bottom_flux)
    end if
    
  end do ! time loop

  call finalize_output()
   
end program

subroutine diffusion_implicit_midpoint(nsoil, zsnso, dt, tbot, zbot, df, hcpct, tg, stc, bottom_temperature_option)

  use routines, only: rosr12
  implicit none
  
  integer               , intent(in)    :: nsoil
  real, dimension(nsoil), intent(in)    :: zsnso  !< thickness of snow/soil layers [m]
  real                  , intent(in)    :: dt     !< time step [s]
  real                  , intent(in)    :: tbot   !< 
  real                  , intent(in)    :: zbot   !< 
  real, dimension(nsoil), intent(in)    :: df     !< thermal conductivity
  real, dimension(nsoil), intent(in)    :: hcpct  !< heat capacity (j/m3/k)
  real                  , intent(in)    :: tg     !< surface temperature (k)
  real, dimension(nsoil), intent(inout) :: stc    !< snow/soil/lake temp. (k)
  integer               , intent(in)    :: bottom_temperature_option  

  ! local variables
  
  real, dimension(0:nsoil)   :: depth_interface  ! depth of interface (negative)
  real, dimension(0:nsoil+1) :: depth_node       ! depth of node (negative)
  
  real, dimension(1:nsoil)   :: dz_interface     ! distance between interfaces
  real, dimension(0:nsoil)   :: dz_node          ! distance between nodes
  
  real, dimension(1:nsoil)   :: heat_capacity    ! heat capacity
  real, dimension(0:nsoil)   :: thermal_cond     ! thermal conductivity
  
  real, dimension(1:nsoil)   :: a,b,c,d,d2       ! tridiagonal terms
  
  integer :: iz
  
  a = huge(1.0)
  b = huge(1.0)
  c = huge(1.0)
  d = huge(1.0)
  d2 = huge(1.0)
  
  heat_capacity = hcpct
  thermal_cond(1:nsoil) = df(1:nsoil)
  thermal_cond(0) = df(1)
  
  depth_interface(0)       = 0.0
  depth_interface(1:nsoil) = zsnso(1:nsoil)

  depth_node(0)       = 0.0
  do iz = 1, nsoil
    depth_node(iz)    = 0.5 * (depth_interface(iz-1) + depth_interface(iz))
  end do
  depth_node(nsoil+1) = zbot
  
  do iz = 1, nsoil
    dz_interface(iz)    = depth_interface(iz-1) - depth_interface(iz)
  end do
  
  do iz = 0, nsoil
    dz_node(iz)         = depth_node(iz) - depth_node(iz+1)
  end do
  
  do iz = 1, nsoil
  
    if(iz == nsoil .and. bottom_temperature_option == 1) then
      b(iz) = 1 + dt/heat_capacity(iz)/dz_interface(iz) *  thermal_cond(iz-1)/dz_node(iz-1) 
    else
      b(iz) = 1 + dt/heat_capacity(iz)/dz_interface(iz) *  thermal_cond(iz-1)/dz_node(iz-1) + dt/heat_capacity(iz)/dz_interface(iz) * thermal_cond(iz)/dz_node(iz)
    end if
    
    if(iz .ne. nsoil) c(iz) = - dt/heat_capacity(iz)/dz_interface(iz) * thermal_cond(iz)/dz_node(iz)

    if(iz .ne. 1)     a(iz) = - dt/heat_capacity(iz)/dz_interface(iz) * thermal_cond(iz-1)/dz_node(iz-1)

    if(iz == 1) then
      d(iz) = stc(iz) + dt/heat_capacity(iz)/dz_interface(iz) *  thermal_cond(iz-1)/dz_node(iz-1) * tg
    elseif(iz == nsoil) then
      if(iz == nsoil .and. bottom_temperature_option == 1) then
        d(iz) = stc(iz)
      else
        d(iz) = stc(iz) + dt/heat_capacity(iz)/dz_interface(iz) *  thermal_cond(iz)/dz_node(iz) * tbot
      end if
    else
      d(iz) = stc(iz)
    end if

  end do
  
  call rosr12 (stc,a,b,c,d,d2,1,nsoil,0)

end subroutine diffusion_implicit_midpoint

subroutine diffusion_implicit_interface(nsoil, depth_midpoint, depth_interface, dt, df, hcpct, tg, stc)

  use routines, only: rosr12
  implicit none
  
  integer               , intent(in)    :: nsoil
  real, dimension(nsoil), intent(in)    :: depth_midpoint   ! depth to midpoints (negative) [m]
  real, dimension(nsoil), intent(in)    :: depth_interface  ! depth to interfaces (negative) [m]
  real                  , intent(in)    :: dt               ! time step [s]
  real, dimension(nsoil), intent(in)    :: df               ! thermal conductivity (W/m2/K)
  real, dimension(nsoil), intent(in)    :: hcpct            ! heat capacity (J/m3/K)
  real                  , intent(in)    :: tg               ! surface temperature (K)
  real, dimension(nsoil), intent(inout) :: stc              ! soil temperature (K)

  ! local variables
  
  real, dimension(1:nsoil)   :: dz_interface     ! distance between interfaces
  real, dimension(1:nsoil)   :: dz_midpoint      ! distance between midpoints
  
  real, dimension(1:nsoil)   :: heat_capacity    ! heat capacity
  real, dimension(1:nsoil)   :: thermal_cond     ! thermal conductivity
  
  real, dimension(1:nsoil)   :: a,b,c,d,d2       ! tridiagonal terms
  
  integer :: iz
  
  a = huge(1.0)
  b = huge(1.0)
  c = huge(1.0)
  d = huge(1.0)
  d2 = huge(1.0)
  
  heat_capacity = hcpct
  thermal_cond(1:nsoil) = df(1:nsoil)
  
  dz_interface(1) = 0 - depth_interface(1)
  do iz = 2, nsoil
    dz_interface(iz)    = depth_interface(iz-1) - depth_interface(iz)
  end do
  
  do iz = 1, nsoil-1
    dz_midpoint(iz)     = depth_midpoint(iz) - depth_midpoint(iz+1)
  end do
  dz_midpoint(nsoil)    = 2.0*(depth_midpoint(iz) - depth_interface(iz))
  
  do iz = 1, nsoil
  
    if(iz == nsoil) then
      b(iz) = 1 + dt/heat_capacity(iz)/dz_midpoint(iz) *  thermal_cond(iz)/dz_interface(iz) 
    else
      b(iz) = 1 + dt/heat_capacity(iz)/dz_midpoint(iz) *  thermal_cond(iz)/dz_interface(iz) + dt/heat_capacity(iz)/dz_midpoint(iz) * thermal_cond(iz+1)/dz_interface(iz+1)
    end if
    
    if(iz .ne. nsoil) c(iz) = - dt/heat_capacity(iz)/dz_midpoint(iz) * thermal_cond(iz+1)/dz_interface(iz+1)

    if(iz .ne. 1)     a(iz) = - dt/heat_capacity(iz)/dz_midpoint(iz) * thermal_cond(iz)/dz_interface(iz)

    if(iz == 1) then
      d(iz) = stc(iz) + dt/heat_capacity(iz)/dz_midpoint(iz) *  thermal_cond(iz)/dz_interface(iz) * tg
    else
      d(iz) = stc(iz)
    end if

  end do
  
  call rosr12 (stc,a,b,c,d,d2,1,nsoil,0)

end subroutine diffusion_implicit_interface

subroutine diffusion_cn_interface(nsoil, depth_midpoint, depth_interface, dt, df, hcpct, tg, tg_previous, stc)

  use routines, only: rosr12
  implicit none
  
  integer               , intent(in)    :: nsoil
  real, dimension(nsoil), intent(in)    :: depth_midpoint   ! depth to midpoints (negative) [m]
  real, dimension(nsoil), intent(in)    :: depth_interface  ! depth to interfaces (negative) [m]
  real                  , intent(in)    :: dt               ! time step [s]
  real, dimension(nsoil), intent(in)    :: df               ! thermal conductivity (W/m2/K)
  real, dimension(nsoil), intent(in)    :: hcpct            ! heat capacity (J/m3/K)
  real                  , intent(in)    :: tg               ! surface temperature (K)
  real                  , intent(in)    :: tg_previous      ! previous surface temperature (K)
  real, dimension(nsoil), intent(inout) :: stc              ! soil temperature (K)

  ! local variables
  
  real, dimension(1:nsoil)   :: dz_interface     ! distance between interfaces
  real, dimension(1:nsoil)   :: dz_midpoint      ! distance between midpoints
  
  real, dimension(1:nsoil)   :: heat_capacity    ! heat capacity
  real, dimension(1:nsoil)   :: thermal_cond     ! thermal conductivity
  
  real, dimension(1:nsoil)   :: a,b,c,d,d2       ! tridiagonal terms
  
  integer :: iz
  
  a = huge(1.0)
  b = huge(1.0)
  c = huge(1.0)
  d = huge(1.0)
  d2 = huge(1.0)
  
  heat_capacity = hcpct
  thermal_cond(1:nsoil) = df(1:nsoil)
  
  dz_interface(1) = 0 - depth_interface(1)
  do iz = 2, nsoil
    dz_interface(iz)    = depth_interface(iz-1) - depth_interface(iz)
  end do
  
  do iz = 1, nsoil-1
    dz_midpoint(iz)     = depth_midpoint(iz) - depth_midpoint(iz+1)
  end do
  dz_midpoint(nsoil)    = 2.0*(depth_midpoint(iz) - depth_interface(iz))
  
  do iz = 1, nsoil
  
    if(iz == nsoil) then
      b(iz) = 1 + dt/heat_capacity(iz)/dz_midpoint(iz)/2.0 *  thermal_cond(iz)/dz_interface(iz) 
    else
      b(iz) = 1 + dt/heat_capacity(iz)/dz_midpoint(iz)/2.0 *  thermal_cond(iz)/dz_interface(iz) + dt/heat_capacity(iz)/dz_midpoint(iz)/2.0 * thermal_cond(iz+1)/dz_interface(iz+1)
    end if
    
    if(iz .ne. nsoil) c(iz) = - dt/heat_capacity(iz)/dz_midpoint(iz)/2.0 * thermal_cond(iz+1)/dz_interface(iz+1)

    if(iz .ne. 1)     a(iz) = - dt/heat_capacity(iz)/dz_midpoint(iz)/2.0 * thermal_cond(iz)/dz_interface(iz)

    if(iz == 1) then
      d(iz) = stc(iz) + dt/heat_capacity(iz)/dz_midpoint(iz)/2.0 *  thermal_cond(iz)/dz_interface(iz) * tg &
        + dt/heat_capacity(iz)/dz_midpoint(iz)/2.0 *  thermal_cond(iz)/dz_interface(iz) * (tg_previous-stc(iz)) &
        + dt/heat_capacity(iz)/dz_midpoint(iz)/2.0 *  thermal_cond(iz+1)/dz_interface(iz+1) * (stc(iz)-stc(iz+1))
    elseif(iz == nsoil) then
      d(iz) = stc(iz) &
        + dt/heat_capacity(iz)/dz_midpoint(iz)/2.0 *  thermal_cond(iz)/dz_interface(iz) * (stc(iz-1)-stc(iz))
    else
      d(iz) = stc(iz) &
        + dt/heat_capacity(iz)/dz_midpoint(iz)/2.0 *  thermal_cond(iz)/dz_interface(iz) * (stc(iz-1)-stc(iz)) &
        + dt/heat_capacity(iz)/dz_midpoint(iz)/2.0 *  thermal_cond(iz+1)/dz_interface(iz+1) * (stc(iz)-stc(iz+1))
    end if

  end do
  
  call rosr12 (stc,a,b,c,d,d2,1,nsoil,0)

end subroutine diffusion_cn_interface
