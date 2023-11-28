module output

  use netcdf

  integer           :: ncid
  integer           :: iret
  integer           :: time_dim
  integer           :: soil_dim
  integer           :: varid

contains

  subroutine initialize_output(output_filename, ntime_in, output_freq, nsoil)
 
    implicit none
    integer       :: ntime_in
    integer       :: ntime
    integer       :: output_freq
    integer       :: nsoil
    character*256 :: output_filename
    real          :: fillvalue = 1.d30
    
    ntime = (ntime_in - 1) / output_freq + 1
 
    iret = nf90_create(trim(output_filename), NF90_CLOBBER, ncid)
     if (iret /= nf90_noerr) call handle_err(iret,"nf90_create")

    iret = nf90_def_dim(ncid, "time", ntime, time_dim)
     if (iret /= nf90_noerr) call handle_err(iret,"define time dimension")

    iret = nf90_def_dim(ncid, "soil_levels", nsoil, soil_dim)
     if (iret /= nf90_noerr) call handle_err(iret,"define soil dimension")

    iret = nf90_def_var(ncid, "time",                  NF90_INT, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define time variable")

    iret = nf90_def_var(ncid, "temperature_ground",    NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define temperature_ground variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define temperature_ground attribute")

    iret = nf90_def_var(ncid, "soil_temperature", NF90_DOUBLE, (/time_dim,soil_dim/), varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define soil_temperature variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define soil_temperature attribute")

    iret = nf90_def_var(ncid, "thermal_cond",     NF90_DOUBLE, (/time_dim,soil_dim/), varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define thermal_cond variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define thermal_cond attribute")

    iret = nf90_def_var(ncid, "heat_capacity",    NF90_DOUBLE, (/time_dim,soil_dim/), varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define heat_capacity variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define heat_capacity attribute")

    iret = nf90_def_var(ncid, "theory_soiltemp",  NF90_DOUBLE, (/time_dim,soil_dim/), varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define theory_soiltemp variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define theory_soiltemp attribute")

    iret = nf90_def_var(ncid, "error_soiltemp",   NF90_DOUBLE, (/time_dim,soil_dim/), varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define error_soiltemp variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define error_soiltemp attribute")
      iret = nf90_put_att(ncid, varid, "long_name", "calculated - theoretical temperature")
       if (iret /= nf90_noerr) call handle_err(iret,"define error_soiltemp long name")

    iret = nf90_def_var(ncid, "energy_balance",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define  variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define  attribute")
      iret = nf90_put_att(ncid, varid, "long_name", "soil column energy balance")
       if (iret /= nf90_noerr) call handle_err(iret,"define energy_balance long name")

    iret = nf90_def_var(ncid, "top_flux",           NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define ground_heat variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define ground_heat attribute")

    iret = nf90_def_var(ncid, "theoretical_top_flux",           NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define theory ground_heat variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define theory ground_heat attribute")

    iret = nf90_def_var(ncid, "error_top_flux",           NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define error ground_heat variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define error ground_heat attribute")

    iret = nf90_def_var(ncid, "bottom_flux",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define bottom_flux variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define bottom_flux attribute")
      iret = nf90_put_att(ncid, varid, "long_name", "soil bottom energy flux (+down)")
       if (iret /= nf90_noerr) call handle_err(iret,"define bottom_flux long name")

    iret = nf90_enddef(ncid)
     if (iret /= nf90_noerr) call handle_err(iret,"ending define mode")
  
   end subroutine initialize_output

   subroutine add_to_output(itime_in, output_freq,&
                            nsoil,                &
                            temperature_ground,   &
                            temperature_soil,     &
                            thermal_conductivity, &
                            heat_capacity,        &
                            ground_heat,          &
                            theoretical_temperature,&
                            energy_balance,       &
                            bottom_flux,          &
                            theoretical_top_flux  )

    implicit none

     integer                :: itime
     integer                :: itime_in
     integer                :: output_freq
     integer                :: nsoil
     real                   :: temperature_ground
     real                   :: ground_heat
     real                   :: energy_balance
     real                   :: bottom_flux
     real                   :: theoretical_top_flux
     real, dimension(nsoil) :: temperature_soil
     real, dimension(nsoil) :: thermal_conductivity
     real, dimension(nsoil) :: heat_capacity
     real, dimension(nsoil) :: theoretical_temperature
     
     itime = itime_in/output_freq
     
     iret = nf90_inq_varid(ncid, "time", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire time variable")
     iret = nf90_put_var(ncid, varid, itime,                 start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put time variable")

     iret = nf90_inq_varid(ncid, "temperature_ground", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire temperature_ground variable")
     iret = nf90_put_var(ncid, varid, temperature_ground,    start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put temperature_ground variable")

     iret = nf90_inq_varid(ncid, "energy_balance", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire energy_balance variable")
     iret = nf90_put_var(ncid, varid, energy_balance ,     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put energy_balance variable")

     iret = nf90_inq_varid(ncid, "top_flux", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire ground_heat variable")
     iret = nf90_put_var(ncid, varid, ground_heat,           start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put ground_heat variable")

     iret = nf90_inq_varid(ncid, "theoretical_top_flux", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire theoretical_ground_heat variable")
     iret = nf90_put_var(ncid, varid, theoretical_top_flux,           start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put theoretical_ground_heat variable")

     iret = nf90_inq_varid(ncid, "bottom_flux", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire bottom_flux variable")
     iret = nf90_put_var(ncid, varid, bottom_flux ,     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put bottom_flux variable")
  
     iret = nf90_inq_varid(ncid, "soil_temperature", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire soil_temperature variable")
     iret = nf90_put_var(ncid, varid,  temperature_soil,     start=(/itime+1,1/), count=(/1,nsoil/))
      if (iret /= nf90_noerr) call handle_err(iret,"put soil_temperature variable")

     iret = nf90_inq_varid(ncid, "thermal_cond", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire thermal_cond variable")
     iret = nf90_put_var(ncid, varid,  thermal_conductivity,     start=(/itime+1,1/), count=(/1,nsoil/))
      if (iret /= nf90_noerr) call handle_err(iret,"put thermal_cond variable")

     iret = nf90_inq_varid(ncid, "heat_capacity", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire heat_capacity variable")
     iret = nf90_put_var(ncid, varid,  heat_capacity,     start=(/itime+1,1/), count=(/1,nsoil/))
      if (iret /= nf90_noerr) call handle_err(iret,"put heat_capacity variable")

     iret = nf90_inq_varid(ncid, "theory_soiltemp", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire theory_soiltemp variable")
     iret = nf90_put_var(ncid, varid,  theoretical_temperature,     start=(/itime+1,1/), count=(/1,nsoil/))
      if (iret /= nf90_noerr) call handle_err(iret,"put theory_soiltemp variable")

  if(itime>0) then
     iret = nf90_inq_varid(ncid, "error_soiltemp", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire error_soiltemp variable")
     iret = nf90_put_var(ncid, varid,  temperature_soil-theoretical_temperature,     start=(/itime+1,1/), count=(/1,nsoil/))
      if (iret /= nf90_noerr) call handle_err(iret,"put error_soiltemp variable")

     iret = nf90_inq_varid(ncid, "error_top_flux", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire error_top_flux variable")
     iret = nf90_put_var(ncid, varid,  ground_heat-theoretical_top_flux,     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put error_top_flux variable")
  end if
  
   end subroutine add_to_output

   subroutine finalize_output()

     implicit none

     iret = nf90_close(ncid)
      if (iret /= nf90_noerr) call handle_err(iret,"closing file")

   end subroutine finalize_output
   
  subroutine handle_err(status,sometext)
    integer, intent ( in) :: status
    character(len=*), optional :: sometext
 
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      if(present(sometext)) print *, "Working on: ",trim(sometext)
      stop "Stopped"
    end if
  end subroutine handle_err

end module output

