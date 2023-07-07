module routines

  TYPE noahmp_parameters
     real, allocatable, dimension(:) :: smcmax       ! porosity (volumetric)
     real, allocatable, dimension(:) :: quartz       ! quartz content
     real                            :: zbot         ! 
     real                            :: csoil        ! 
     logical                         :: urban_flag   ! 
     
  END TYPE noahmp_parameters
  
  integer :: opt_tbot
  integer :: opt_stc
  
  integer, parameter :: kind_phys = 8

  real (kind=kind_phys), parameter :: tfrz   = 273.16    !< freezing/melting point (k)
  real (kind=kind_phys), parameter :: cwat   = 4.188e06  !< specific heat capacity of water (j/m3/k)
  real (kind=kind_phys), parameter :: cice   = 2.094e06  !< specific heat capacity of ice (j/m3/k)
  real (kind=kind_phys), parameter :: cpair  = 1004.64   !< heat capacity dry air at const pres (j/kg/k)
  real (kind=kind_phys), parameter :: tkwat  = 0.6       !< thermal conductivity of water (w/m/k)
  real (kind=kind_phys), parameter :: tkice  = 2.2       !< thermal conductivity of ice (w/m/k)
  real (kind=kind_phys), parameter :: tkair  = 0.023     !< thermal conductivity of air (w/m/k) (not used mb: 20140718)
  real (kind=kind_phys), parameter :: denh2o = 1000.     !< density of water (kg/m3)
  real (kind=kind_phys), parameter :: denice = 917.      !< density of ice (kg/m3)

contains

!== begin thermoprop ===============================================================================

!>\ingroup NoahMP_LSM
  subroutine thermoprop (parameters,nsoil   ,nsnow   ,isnow   ,ist     ,dzsnso  , & !in
                         dt      ,snowh   ,snice   ,snliq   , & !in
                         smc     ,sh2o    ,tg      ,stc     ,ur      , & !in
                         lat     ,z0m     ,zlvl    ,vegtyp  , & !in
                         df      ,hcpct   ,snicev  ,snliqv  ,epore   , & !out
                         fact    )                                       !out
! ------------------------------------------------------------------------------------------------- 
  implicit none
! --------------------------------------------------------------------------------------------------
! inputs
  type (noahmp_parameters), intent(in) :: parameters      !< 
  integer                        , intent(in)  :: nsoil   !< number of soil layers
  integer                        , intent(in)  :: nsnow   !< maximum no. of snow layers        
  integer                        , intent(in)  :: isnow   !< actual no. of snow layers
  integer                        , intent(in)  :: ist     !< surface type
  real (kind=kind_phys)                           , intent(in)  :: dt      !< time step [s]
  real (kind=kind_phys), dimension(-nsnow+1:    0), intent(in)  :: snice   !< snow ice mass (kg/m2)
  real (kind=kind_phys), dimension(-nsnow+1:    0), intent(in)  :: snliq   !< snow liq mass (kg/m2)
  real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(in)  :: dzsnso  !< thickness of snow/soil layers [m]
  real (kind=kind_phys), dimension(       1:nsoil), intent(in)  :: smc     !< soil moisture (ice + liq.) [m3/m3]
  real (kind=kind_phys), dimension(       1:nsoil), intent(in)  :: sh2o    !< liquid soil moisture [m3/m3]
  real (kind=kind_phys)                           , intent(in)  :: snowh   !< snow height [m]
  real (kind=kind_phys),                            intent(in)  :: tg      !< surface temperature (k)
  real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(in)  :: stc     !< snow/soil/lake temp. (k)
  real (kind=kind_phys),                            intent(in)  :: ur      !< wind speed at zlvl (m/s)
  real (kind=kind_phys),                            intent(in)  :: lat     !< latitude (radians)
  real (kind=kind_phys),                            intent(in)  :: z0m     !< roughness length (m)
  real (kind=kind_phys),                            intent(in)  :: zlvl    !< reference height (m)
  integer                                         , intent(in)  :: vegtyp  !< vegtyp type

! outputs
  real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(out) :: df      !< thermal conductivity [w/m/k]
  real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(out) :: hcpct   !< heat capacity [j/m3/k]
  real (kind=kind_phys), dimension(-nsnow+1:    0), intent(out) :: snicev  !< partial volume of ice [m3/m3]
  real (kind=kind_phys), dimension(-nsnow+1:    0), intent(out) :: snliqv  !< partial volume of liquid water [m3/m3]
  real (kind=kind_phys), dimension(-nsnow+1:    0), intent(out) :: epore   !< effective porosity [m3/m3]
  real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(out) :: fact    !< computing energy for phase change
! --------------------------------------------------------------------------------------------------
! locals

  integer :: iz
  real (kind=kind_phys), dimension(-nsnow+1:    0)              :: cvsno   !volumetric specific heat (j/m3/k)
  real (kind=kind_phys), dimension(-nsnow+1:    0)              :: tksno   !snow thermal conductivity (j/m3/k)
  real (kind=kind_phys), dimension(       1:nsoil)              :: sice    !soil ice content
  real (kind=kind_phys), parameter :: sbeta = -2.0
! --------------------------------------------------------------------------------------------------

! compute snow thermal conductivity and heat capacity

    call csnow (parameters,isnow   ,nsnow   ,nsoil   ,snice   ,snliq   ,dzsnso  , & !in
                tksno   ,cvsno   ,snicev  ,snliqv  ,epore   )   !out

    do iz = isnow+1, 0
      df   (iz) = tksno(iz)
      hcpct(iz) = cvsno(iz)
    end do

! compute soil thermal properties

    do  iz = 1, nsoil
       sice(iz)  = smc(iz) - sh2o(iz)
       hcpct(iz) = sh2o(iz)*cwat + (1.0-parameters%smcmax(iz))*parameters%csoil &
                + (parameters%smcmax(iz)-smc(iz))*cpair + sice(iz)*cice
       call tdfcnd (parameters,iz,df(iz), smc(iz), sh2o(iz))
    end do
       
    if ( parameters%urban_flag ) then
       do iz = 1,nsoil
         df(iz) = 3.24
       end do
    endif

! heat flux reduction effect from the overlying green canopy, adapted from 
! section 2.1.2 of peters-lidard et al. (1997, jgr, vol 102(d4)).
! not in use because of the separation of the canopy layer from the ground.
! but this may represent the effects of leaf litter (niu comments)
!       df1 = df1 * exp (sbeta * shdfac)

! compute lake thermal properties 
! (no consideration of turbulent mixing for this version)

    if(ist == 2) then
       do iz = 1, nsoil 
         if(stc(iz) > tfrz) then
            hcpct(iz) = cwat
            df(iz)    = tkwat  !+ keddy * cwat 
         else
            hcpct(iz) = cice
            df(iz)    = tkice 
         end if
       end do
    end if

! combine a temporary variable used for melting/freezing of snow and frozen soil

    do iz = isnow+1,nsoil
     fact(iz) = dt/(hcpct(iz)*dzsnso(iz))
    end do

! snow/soil interface

    if(isnow == 0) then
       df(1) = (df(1)*dzsnso(1)+0.35*snowh)      / (snowh    +dzsnso(1)) 
    else
       df(1) = (df(1)*dzsnso(1)+df(0)*dzsnso(0)) / (dzsnso(0)+dzsnso(1))
    end if


  end subroutine thermoprop

!== begin csnow ====================================================================================

!>\ingroup NoahMP_LSM
!! snow bulk density,volumetric capacity, and thermal conductivity
  subroutine csnow (parameters,isnow   ,nsnow   ,nsoil   ,snice   ,snliq   ,dzsnso  , & !in
                    tksno   ,cvsno   ,snicev  ,snliqv  ,epore   )   !out
! --------------------------------------------------------------------------------------------------
! snow bulk density,volumetric capacity, and thermal conductivity
!---------------------------------------------------------------------------------------------------
  implicit none
!---------------------------------------------------------------------------------------------------
! inputs

  type (noahmp_parameters), intent(in) :: parameters     !< 
  integer,                          intent(in) :: isnow  !< number of snow layers (-)            
  integer                        ,  intent(in) :: nsnow  !< maximum no. of snow layers        
  integer                        ,  intent(in) :: nsoil  !< number of soil layers
  real (kind=kind_phys), dimension(-nsnow+1:    0),  intent(in) :: snice  !< snow ice mass (kg/m2)
  real (kind=kind_phys), dimension(-nsnow+1:    0),  intent(in) :: snliq  !< snow liq mass (kg/m2) 
  real (kind=kind_phys), dimension(-nsnow+1:nsoil),  intent(in) :: dzsnso !< snow/soil layer thickness [m]

! outputs

  real (kind=kind_phys), dimension(-nsnow+1:    0), intent(out) :: cvsno  !< volumetric specific heat (j/m3/k)
  real (kind=kind_phys), dimension(-nsnow+1:    0), intent(out) :: tksno  !< thermal conductivity (w/m/k)
  real (kind=kind_phys), dimension(-nsnow+1:    0), intent(out) :: snicev !< partial volume of ice [m3/m3]
  real (kind=kind_phys), dimension(-nsnow+1:    0), intent(out) :: snliqv !< partial volume of liquid water [m3/m3]
  real (kind=kind_phys), dimension(-nsnow+1:    0), intent(out) :: epore  !< effective porosity [m3/m3]

! locals

  integer :: iz
  real (kind=kind_phys), dimension(-nsnow+1:    0) :: bdsnoi  !bulk density of snow(kg/m3)

!---------------------------------------------------------------------------------------------------
! thermal capacity of snow

  do iz = isnow+1, 0
      snicev(iz)   = min(1., snice(iz)/(dzsnso(iz)*denice) )
      epore(iz)    = 1. - snicev(iz)
      snliqv(iz)   = min(epore(iz),snliq(iz)/(dzsnso(iz)*denh2o))
  enddo

  do iz = isnow+1, 0
      bdsnoi(iz) = (snice(iz)+snliq(iz))/dzsnso(iz)
      cvsno(iz) = cice*snicev(iz)+cwat*snliqv(iz)
!      cvsno(iz) = 0.525e06                          ! constant
  enddo

! thermal conductivity of snow

  do iz = isnow+1, 0
!     tksno(iz) = 3.2217e-6*bdsnoi(iz)**2.           ! stieglitz(yen,1965)
!    tksno(iz) = 2e-2+2.5e-6*bdsnoi(iz)*bdsnoi(iz)   ! anderson, 1976
!    tksno(iz) = 0.35                                ! constant
    tksno(iz) = 2.576e-6*bdsnoi(iz)**2. + 0.074    ! verseghy (1991)
!    tksno(iz) = 2.22*(bdsnoi(iz)/1000.)**1.88      ! douvill(yen, 1981)
  enddo

  end subroutine csnow

!== begin tdfcnd ===================================================================================

!>\ingroup NoahMP_LSM
!! calculate thermal diffusivity and conductivity of the soil. peters-lidard
!! approach (peters-lidard et al., 1998)
  subroutine tdfcnd (parameters, isoil, df, smc, sh2o)
! --------------------------------------------------------------------------------------------------
! calculate thermal diffusivity and conductivity of the soil.
! peters-lidard approach (peters-lidard et al., 1998)
! --------------------------------------------------------------------------------------------------
! code history:
! june 2001 changes: frozen soil condition.
! --------------------------------------------------------------------------------------------------
    implicit none
  type (noahmp_parameters), intent(in) :: parameters
    integer, intent(in)                     :: isoil  !< soil layer
    real (kind=kind_phys), intent(in)       :: smc    !< total soil water
    real (kind=kind_phys), intent(in)       :: sh2o   !< liq. soil water
    real (kind=kind_phys), intent(out)      :: df     !< thermal diffusivity

! local variables
    real (kind=kind_phys)  :: ake
    real (kind=kind_phys)  :: gammd
    real (kind=kind_phys)  :: thkdry
    real (kind=kind_phys)  :: thko     ! thermal conductivity for other soil components         
    real (kind=kind_phys)  :: thkqtz   ! thermal conductivity for quartz
    real (kind=kind_phys)  :: thksat   ! 
    real (kind=kind_phys)  :: thks     ! thermal conductivity for the solids
    real (kind=kind_phys)  :: thkw     ! water thermal conductivity
    real (kind=kind_phys)  :: satratio
    real (kind=kind_phys)  :: xu
    real (kind=kind_phys)  :: xunfroz
! --------------------------------------------------------------------------------------------------
! we now get quartz as an input argument (set in routine redprm):
!      data quartz /0.82, 0.10, 0.25, 0.60, 0.52,
!     &             0.35, 0.60, 0.40, 0.82/
! --------------------------------------------------------------------------------------------------
! if the soil has any moisture content compute a partial sum/product
! otherwise use a constant value which works well with most soils
! --------------------------------------------------------------------------------------------------
!  quartz ....quartz content (soil type dependent)
! --------------------------------------------------------------------------------------------------
! use as in peters-lidard, 1998 (modif. from johansen, 1975).

!                                  pablo grunmann, 08/17/98
! refs.:
!      farouki, o.t.,1986: thermal properties of soils. series on rock
!              and soil mechanics, vol. 11, trans tech, 136 pp.
!      johansen, o., 1975: thermal conductivity of soils. ph.d. thesis,
!              university of trondheim,
!      peters-lidard, c. d., et al., 1998: the effect of soil thermal
!              conductivity parameterization on surface energy fluxes
!              and temperatures. journal of the atmospheric sciences,
!              vol. 55, pp. 1209-1224.
! --------------------------------------------------------------------------------------------------
! needs parameters
! porosity(soil type):
!      poros = smcmax
! saturation ratio:
! parameters  w/(m.k)
    satratio = smc / parameters%smcmax(isoil)
    thkw = 0.57
!      if (quartz .le. 0.2) thko = 3.0
    thko = 2.0
! solids' conductivity
! quartz' conductivity
    thkqtz = 7.7

! unfrozen fraction (from 1., i.e., 100%liquid, to 0. (100% frozen))
    thks = (thkqtz ** parameters%quartz(isoil))* (thko ** (1. - parameters%quartz(isoil)))

! unfrozen volume for saturation (porosity*xunfroz)
    xunfroz = 1.0                       ! prevent divide by zero (suggested by d. mocko)
    if(smc > 0.) xunfroz = sh2o / smc
! saturated thermal conductivity
    xu = xunfroz * parameters%smcmax(isoil)

! dry density in kg/m3
    thksat = thks ** (1. - parameters%smcmax(isoil))* tkice ** (parameters%smcmax(isoil) - xu)* thkw **   &
         (xu)

! dry thermal conductivity in w.m-1.k-1
    gammd = (1. - parameters%smcmax(isoil))*2700.

    thkdry = (0.135* gammd+ 64.7)/ (2700. - 0.947* gammd)
! frozen
    if ( (sh2o + 0.0005) <  smc ) then
       ake = satratio
! unfrozen
! range of validity for the kersten number (ake)
    else

! kersten number (using "fine" formula, valid for soils containing at
! least 5% of particles with diameter less than 2.e-6 meters.)
! (for "coarse" formula, see peters-lidard et al., 1998).

       if ( satratio >  0.1 ) then

          ake = log10 (satratio) + 1.0

! use k = kdry
       else

          ake = 0.0
       end if
!  thermal conductivity

    end if

    df = ake * (thksat - thkdry) + thkdry


  end subroutine tdfcnd

!== begin tsnosoi ==================================================================================

!>\ingroup NoahMP_LSM
!! compute snow (up to 3l) and soil (4l) temperature. note that snow
!! temperatures during melting season may exceed melting point (tfrz) but later
!! in phasechange subroutine the snow temperatures are reset to tfrz for melting
!! snow.
  subroutine tsnosoi (parameters,ice     ,nsoil   ,nsnow   ,isnow   ,ist     , & !in
                      tbot    ,zsnso   ,ssoil   ,df      ,hcpct   , & !in
                      sag     ,dt      ,snowh   ,dzsnso  , & !in
                      tg      ,iloc    ,jloc    ,                   & !in
                      stc     ,errmsg  ,errflg)                       !inout
! --------------------------------------------------------------------------------------------------
! compute snow (up to 3l) and soil (4l) temperature. note that snow temperatures
! during melting season may exceed melting point (tfrz) but later in phasechange
! subroutine the snow temperatures are reset to tfrz for melting snow.
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
!input

  type (noahmp_parameters), intent(in) :: parameters       !<
    integer,                         intent(in)  :: iloc   !<
    integer,                         intent(in)  :: jloc   !<
    integer,                         intent(in)  :: ice    !<
    integer,                         intent(in)  :: nsoil  !< no of soil layers (4)
    integer,                         intent(in)  :: nsnow  !< maximum no of snow layers (3)
    integer,                         intent(in)  :: isnow  !< actual no of snow layers
    integer,                         intent(in)  :: ist    !< surface type

    real (kind=kind_phys),                            intent(in)  :: dt     !< time step (s)
    real (kind=kind_phys),                            intent(in)  :: tbot   !< 
    real (kind=kind_phys),                            intent(in)  :: ssoil  !< ground heat flux (w/m2)
    real (kind=kind_phys),                            intent(in)  :: sag    !< solar rad. absorbed by ground (w/m2)
    real (kind=kind_phys),                            intent(in)  :: snowh  !< snow depth (m)
    real (kind=kind_phys),                            intent(in)  :: tg     !< ground temperature (k)
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(in)  :: zsnso  !< layer-bot. depth from snow surf.(m)
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(in)  :: dzsnso !< snow/soil layer thickness (m)
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(in)  :: df     !< thermal conductivity
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(in)  :: hcpct  !< heat capacity (j/m3/k)

!input and output

    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(inout) :: stc  !<
    character(len=*)               , intent(inout) :: errmsg
    integer                        , intent(inout) :: errflg

!local

    integer                                      :: iz
    real (kind=kind_phys)                                         :: zbotsno   !zbot from snow surface
    real (kind=kind_phys), dimension(-nsnow+1:nsoil)              :: ai, bi, ci, rhsts
    real (kind=kind_phys)                                         :: eflxb !energy influx from soil bottom (w/m2)
    real (kind=kind_phys), dimension(-nsnow+1:nsoil)              :: phi   !light through water (w/m2)

    real (kind=kind_phys), dimension(-nsnow+1:nsoil) :: tbeg
    real (kind=kind_phys)                            :: err_est !heat storage error  (w/m2)
    real (kind=kind_phys)                            :: ssoil2  !ground heat flux (w/m2) (for energy check)
    real (kind=kind_phys)                            :: eflxb2  !heat flux from the bottom (w/m2) (for energy check)
    character(len=256)              :: message
! ----------------------------------------------------------------------
! compute solar penetration through water, needs more work

    phi(isnow+1:nsoil) = 0.

! adjust zbot from soil surface to zbotsno from snow surface

    zbotsno = parameters%zbot - snowh    !from snow surface

! snow/soil heat storage for energy balance check

    do iz = isnow+1, nsoil
       tbeg(iz) = stc(iz)
    enddo

! compute soil temperatures

      call hrt   (parameters,nsnow     ,nsoil     ,isnow     ,zsnso     , &
                  stc       ,tbot      ,zbotsno   ,dt        , &
                  df        ,hcpct     ,ssoil     ,phi       , &
                  ai        ,bi        ,ci        ,rhsts     , &
                  eflxb     )

      call hstep (parameters,nsnow     ,nsoil     ,isnow     ,dt        , &
                  ai        ,bi        ,ci        ,rhsts     , &
                  stc       ) 

! update ground heat flux just for energy check, but not for final output
! otherwise, it would break the surface energy balance

    if(opt_tbot == 1) then
       eflxb2  = 0.
    else if(opt_tbot == 2) then
       eflxb2  = df(nsoil)*(tbot-stc(nsoil)) / &
            (0.5*(zsnso(nsoil-1)+zsnso(nsoil)) - zbotsno)
    end if

    ! skip the energy balance check for now, until we can make it work
    ! right for small time steps.
    return

! energy balance check

    err_est = 0.0
    do iz = isnow+1, nsoil
       err_est = err_est + (stc(iz)-tbeg(iz)) * dzsnso(iz) * hcpct(iz) / dt
    enddo

    if (opt_stc == 1 .or. opt_stc == 3) then   ! semi-implicit
       err_est = err_est - (ssoil +eflxb)
    else                     ! full-implicit
       ssoil2 = df(isnow+1)*(tg-stc(isnow+1))/(0.5*dzsnso(isnow+1))   !m. barlage
       err_est = err_est - (ssoil2+eflxb2)
    endif

    if (abs(err_est) > 1.) then    ! w/m2
       write(message,*) 'tsnosoi is losing(-)/gaining(+) false energy',err_est,' w/m2'
       errmsg = trim(message)
       write(message,'(i6,1x,i6,1x,i3,f18.13,5f20.12)') &
            iloc, jloc, ist,err_est,ssoil,snowh,tg,stc(isnow+1),eflxb
       errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
       !niu      stop
    end if

  end subroutine tsnosoi

!== begin hrt ======================================================================================

!>\ingroup NoahMP_LSM
!! calculate the right hand side of the time tendency term of the soil
!! thermal diffusion equation. also to compute (prepare) the matrix 
!! coefficients for the tri-diagonal matrix of the implicit time scheme.
  subroutine hrt (parameters,nsnow     ,nsoil     ,isnow     ,zsnso     , &
                  stc       ,tbot      ,zbot      ,dt        , &
                  df        ,hcpct     ,ssoil     ,phi       , &
                  ai        ,bi        ,ci        ,rhsts     , &
                  botflx    )
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! calculate the right hand side of the time tendency term of the soil
! thermal diffusion equation.  also to compute ( prepare ) the matrix
! coefficients for the tri-diagonal matrix of the implicit time scheme.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters       !<
    integer,                         intent(in)  :: nsoil  !< no of soil layers (4)
    integer,                         intent(in)  :: nsnow  !< maximum no of snow layers (3)
    integer,                         intent(in)  :: isnow  !, actual no of snow layers
    real (kind=kind_phys),                            intent(in)  :: tbot   !< bottom soil temp. at zbot (k)
    real (kind=kind_phys),                            intent(in)  :: zbot   !< depth of lower boundary condition (m)
                                                                            !! from soil surface not snow surface
    real (kind=kind_phys),                            intent(in)  :: dt     !< time step (s)
    real (kind=kind_phys),                            intent(in)  :: ssoil  !< ground heat flux (w/m2)
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(in)  :: zsnso  !< depth of layer-bottom of snow/soil (m)
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(in)  :: stc    !< snow/soil temperature (k)
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(in)  :: df     !< thermal conductivity [w/m/k]
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(in)  :: hcpct  !< heat capacity [j/m3/k]
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(in)  :: phi    !< light through water (w/m2)

! output

    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(out) :: rhsts  !< right-hand side of the matrix
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(out) :: ai     !< left-hand side coefficient
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(out) :: bi     !< left-hand side coefficient
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(out) :: ci     !< left-hand side coefficient
    real (kind=kind_phys),                            intent(out) :: botflx !< energy influx from soil bottom (w/m2)

! local

    integer                                      :: k
    real (kind=kind_phys), dimension(-nsnow+1:nsoil)              :: ddz
    real (kind=kind_phys), dimension(-nsnow+1:nsoil)              :: dz
    real (kind=kind_phys), dimension(-nsnow+1:nsoil)              :: denom
    real (kind=kind_phys), dimension(-nsnow+1:nsoil)              :: dtsdz
    real (kind=kind_phys), dimension(-nsnow+1:nsoil)              :: eflux
    real (kind=kind_phys)                                         :: temp1
! ----------------------------------------------------------------------

    do k = isnow+1, nsoil
        if (k == isnow+1) then
           denom(k)  = - zsnso(k) * hcpct(k)
           temp1     = - zsnso(k+1)
           ddz(k)    = 2.0 / temp1
           dtsdz(k)  = 2.0 * (stc(k) - stc(k+1)) / temp1
           eflux(k)  = df(k) * dtsdz(k) - ssoil - phi(k)
        else if (k < nsoil) then
           denom(k)  = (zsnso(k-1) - zsnso(k)) * hcpct(k)
           temp1     = zsnso(k-1) - zsnso(k+1)
           ddz(k)    = 2.0 / temp1
           dtsdz(k)  = 2.0 * (stc(k) - stc(k+1)) / temp1
           eflux(k)  = (df(k)*dtsdz(k) - df(k-1)*dtsdz(k-1)) - phi(k)
        else if (k == nsoil) then
           denom(k)  = (zsnso(k-1) - zsnso(k)) * hcpct(k)
           temp1     =  zsnso(k-1) - zsnso(k)
           if(opt_tbot == 1) then
               botflx     = 0. 
           end if
           if(opt_tbot == 2) then
               dtsdz(k)  = (stc(k) - tbot) / ( 0.5*(zsnso(k-1)+zsnso(k)) - zbot)
               botflx    = -df(k) * dtsdz(k)
           end if
           eflux(k)  = (-botflx - df(k-1)*dtsdz(k-1) ) - phi(k)
        end if
    end do

    do k = isnow+1, nsoil
        if (k == isnow+1) then
           ai(k)    =   0.0
           ci(k)    = - df(k)   * ddz(k) / denom(k)
           if (opt_stc == 1 .or. opt_stc == 3 ) then
              bi(k) = - ci(k)
           end if                                        
           if (opt_stc == 2) then
              bi(k) = - ci(k) + df(k)/(0.5*zsnso(k)*zsnso(k)*hcpct(k))
           end if
        else if (k < nsoil) then
           ai(k)    = - df(k-1) * ddz(k-1) / denom(k) 
           ci(k)    = - df(k  ) * ddz(k  ) / denom(k) 
           bi(k)    = - (ai(k) + ci (k))
        else if (k == nsoil) then
           ai(k)    = - df(k-1) * ddz(k-1) / denom(k) 
           ci(k)    = 0.0
           bi(k)    = - (ai(k) + ci(k))
        end if
           rhsts(k)  = eflux(k)/ (-denom(k))
    end do

  end subroutine hrt

!== begin hstep ====================================================================================

!>\ingroup NoahMP_LSM
!! calculate/update the soil temperature fields.
  subroutine hstep (parameters,nsnow     ,nsoil     ,isnow     ,dt        ,  &
                    ai        ,bi        ,ci        ,rhsts     ,  &
                    stc       )  
! ----------------------------------------------------------------------
! calculate/update the soil temperature field.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
    integer,                         intent(in)    :: nsoil
    integer,                         intent(in)    :: nsnow
    integer,                         intent(in)    :: isnow
    real (kind=kind_phys),                            intent(in)    :: dt

! output & input
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(inout) :: rhsts
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(inout) :: ai
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(inout) :: bi
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(inout) :: ci
    real (kind=kind_phys), dimension(-nsnow+1:nsoil), intent(inout) :: stc

! local
    integer                                        :: k
    real (kind=kind_phys), dimension(-nsnow+1:nsoil)                :: rhstsin
    real (kind=kind_phys), dimension(-nsnow+1:nsoil)                :: ciin
! ----------------------------------------------------------------------

    do k = isnow+1,nsoil
       rhsts(k) =   rhsts(k) * dt
       ai(k)    =      ai(k) * dt
       bi(k)    = 1. + bi(k) * dt
       ci(k)    =      ci(k) * dt
    end do

! copy values for input variables before call to rosr12

    do k = isnow+1,nsoil
       rhstsin(k) = rhsts(k)
       ciin(k)    = ci(k)
    end do

! solve the tri-diagonal matrix equation

    call rosr12 (ci,ai,bi,ciin,rhstsin,rhsts,isnow+1,nsoil,nsnow)

! update snow & soil temperature

    do k = isnow+1,nsoil
       stc (k) = stc (k) + ci (k)
    end do

  end subroutine hstep

!== begin rosr12 ===================================================================================

  SUBROUTINE ROSR12 (P,A,B,C,D,DELTA,NTOP,NSOIL,NSNOW)
! ----------------------------------------------------------------------
! SUBROUTINE ROSR12
! ----------------------------------------------------------------------
! INVERT (SOLVE) THE TRI-DIAGONAL MATRIX PROBLEM SHOWN BELOW:
! ###                                            ### ###  ###   ###  ###
! #B(1), C(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
! #A(2), B(2), C(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
! # 0  , A(3), B(3), C(3),  0  ,   . . .  ,    0   # #      #   # D(3) #
! # 0  ,  0  , A(4), B(4), C(4),   . . .  ,    0   # # P(4) #   # D(4) #
! # 0  ,  0  ,  0  , A(5), B(5),   . . .  ,    0   # # P(5) #   # D(5) #
! # .                                          .   # #  .   # = #   .  #
! # .                                          .   # #  .   #   #   .  #
! # .                                          .   # #  .   #   #   .  #
! # 0  , . . . , 0 , A(M-2), B(M-2), C(M-2),   0   # #P(M-2)#   #D(M-2)#
! # 0  , . . . , 0 ,   0   , A(M-1), B(M-1), C(M-1)# #P(M-1)#   #D(M-1)#
! # 0  , . . . , 0 ,   0   ,   0   ,  A(M) ,  B(M) # # P(M) #   # D(M) #
! ###                                            ### ###  ###   ###  ###
! ----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: NTOP           
    INTEGER, INTENT(IN)   :: NSOIL,NSNOW
    INTEGER               :: K, KK

    REAL, DIMENSION(-NSNOW+1:NSOIL),INTENT(IN):: A, B, D
    REAL, DIMENSION(-NSNOW+1:NSOIL),INTENT(INOUT):: C,P,DELTA

! ----------------------------------------------------------------------
! INITIALIZE EQN COEF C FOR THE LOWEST SOIL LAYER
! ----------------------------------------------------------------------
    C (NSOIL) = 0.0
    P (NTOP) = - C (NTOP) / B (NTOP)
! ----------------------------------------------------------------------
! SOLVE THE COEFS FOR THE 1ST SOIL LAYER
! ----------------------------------------------------------------------
    DELTA (NTOP) = D (NTOP) / B (NTOP)
! ----------------------------------------------------------------------
! SOLVE THE COEFS FOR SOIL LAYERS 2 THRU NSOIL
! ----------------------------------------------------------------------
    DO K = NTOP+1,NSOIL
       P (K) = - C (K) * ( 1.0 / (B (K) + A (K) * P (K -1)) )
       DELTA (K) = (D (K) - A (K)* DELTA (K -1))* (1.0/ (B (K) + A (K)&
            * P (K -1)))
    END DO
! ----------------------------------------------------------------------
! SET P TO DELTA FOR LOWEST SOIL LAYER
! ----------------------------------------------------------------------
    P (NSOIL) = DELTA (NSOIL)
! ----------------------------------------------------------------------
! ADJUST P FOR SOIL LAYERS 2 THRU NSOIL
! ----------------------------------------------------------------------
    DO K = NTOP+1,NSOIL
       KK = NSOIL - K + (NTOP-1) + 1
       P (KK) = P (KK) * P (KK +1) + DELTA (KK)
    END DO
! ----------------------------------------------------------------------
  END SUBROUTINE ROSR12

  subroutine noahmp_options(iopt_tbot, iopt_stc )

  implicit none

  integer,  intent(in) :: iopt_tbot !< lower boundary of soil temperature (1->zero-flux; 2->noah)

  integer,  intent(in) :: iopt_stc  !< snow/soil temperature time scheme (only layer 1)
                                    !! 1 -> semi-implicit; 2 -> full implicit (original noah)

! -------------------------------------------------------------------------------------------------

  opt_tbot = iopt_tbot 
  opt_stc  = iopt_stc
  
  end subroutine noahmp_options

end module routines
