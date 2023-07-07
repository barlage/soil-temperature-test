module output

  use netcdf

  integer           :: ncid
  integer           :: iret
  integer           :: time_dim
  integer           :: varid

contains

  subroutine initialize_output(output_filename, ntime, nsoil)
 
    implicit none
    integer       :: ntime
    integer       :: nsoil
    character*256 :: output_filename
    real          :: fillvalue = -999.0
 
    iret = nf90_create(trim(output_filename), NF90_CLOBBER, ncid)
     if (iret /= nf90_noerr) call handle_err(iret,"nf90_create")

    iret = nf90_def_dim(ncid, "time", ntime, time_dim)
     if (iret /= nf90_noerr) call handle_err(iret,"define time dimension")

    iret = nf90_def_var(ncid, "time",                  NF90_INT, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define time variable")

    iret = nf90_def_var(ncid, "temperature_ground",    NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define temperature_ground variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define temperature_ground attribute")

    iret = nf90_def_var(ncid, "soil_temperature_lev1", NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define soil_temperature_lev1 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define soil_temperature_lev1 attribute")

    iret = nf90_def_var(ncid, "soil_temperature_lev2", NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define soil_temperature_lev2 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define soil_temperature_lev2 attribute")

    iret = nf90_def_var(ncid, "soil_temperature_lev3", NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define soil_temperature_lev3 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define soil_temperature_lev3 attribute")

    iret = nf90_def_var(ncid, "soil_temperature_lev4", NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define soil_temperature_lev4 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define soil_temperature_lev4 attribute")

    iret = nf90_def_var(ncid, "ground_heat",           NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define ground_heat variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define ground_heat attribute")

    iret = nf90_def_var(ncid, "soil_heat_lev1",        NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define soil_heat_lev1 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define soil_heat_lev1 attribute")

    iret = nf90_def_var(ncid, "soil_heat_lev2",        NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define soil_heat_lev2 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define soil_heat_lev2 attribute")

    iret = nf90_def_var(ncid, "soil_heat_lev3",        NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define soil_heat_lev3 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define soil_heat_lev3 attribute")

    iret = nf90_def_var(ncid, "soil_heat_lev4",        NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define soil_heat_lev4 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define soil_heat_lev4 attribute")

    iret = nf90_def_var(ncid, "df_lev1",               NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define df_lev1 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define df_lev1 attribute")

    iret = nf90_def_var(ncid, "df_lev2",               NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define df_lev2 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define df_lev2 attribute")

    iret = nf90_def_var(ncid, "df_lev3",               NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define df_lev3 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define df_lev3 attribute")

    iret = nf90_def_var(ncid, "df_lev4",               NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define df_lev4 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define df_lev4 attribute")

    iret = nf90_def_var(ncid, "hcpct_lev1",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define hcpct_lev1 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define hcpct_lev1 attribute")

    iret = nf90_def_var(ncid, "hcpct_lev2",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define hcpct_lev2 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define hcpct_lev2 attribute")

    iret = nf90_def_var(ncid, "hcpct_lev3",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define hcpct_lev3 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define hcpct_lev3 attribute")

    iret = nf90_def_var(ncid, "hcpct_lev4",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define hcpct_lev4 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define hcpct_lev4 attribute")

    iret = nf90_def_var(ncid, "theory_lev1",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define theory_lev1 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define theory_lev1 attribute")

    iret = nf90_def_var(ncid, "theory_lev2",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define theory_lev2 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define theory_lev2 attribute")

    iret = nf90_def_var(ncid, "theory_lev3",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define theory_lev3 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define theory_lev3 attribute")

    iret = nf90_def_var(ncid, "theory_lev4",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define theory_lev4 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define theory_lev4 attribute")

    iret = nf90_def_var(ncid, "diff_lev1",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev1 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev1 attribute")
      iret = nf90_put_att(ncid, varid, "long_name", "calculated - theoretical temperature")
       if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev1 long name")

    iret = nf90_def_var(ncid, "diff_lev2",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev2 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev2 attribute")
      iret = nf90_put_att(ncid, varid, "long_name", "calculated - theoretical temperature")
       if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev2 long name")

    iret = nf90_def_var(ncid, "diff_lev3",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev3 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev3 attribute")
      iret = nf90_put_att(ncid, varid, "long_name", "calculated - theoretical temperature")
       if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev3 long name")

    iret = nf90_def_var(ncid, "diff_lev4",            NF90_DOUBLE, time_dim, varid)
     if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev4 variable")
      iret = nf90_put_att(ncid, varid, "_FillValue", fillvalue)
       if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev4 attribute")
      iret = nf90_put_att(ncid, varid, "long_name", "calculated - theoretical temperature")
       if (iret /= nf90_noerr) call handle_err(iret,"define diff_lev4 long name")

    iret = nf90_enddef(ncid)
     if (iret /= nf90_noerr) call handle_err(iret,"ending define mode")
  
   end subroutine initialize_output

   subroutine add_to_output(itime, nsoil,         &
                            temperature_ground,   &
                            temperature_soil,     &
                            thermal_conductivity, &
                            heat_capacity,        &
                            ground_heat,          &
                            theoretical_temperature)!,          &
                           ! soil_interface_flux)

    implicit none
     integer                :: itime
     integer                :: nsoil
     real                   :: temperature_ground
     real                   :: ground_heat
     real, dimension(nsoil) :: temperature_soil
     real, dimension(nsoil) :: thermal_conductivity
     real, dimension(nsoil) :: heat_capacity
     real, dimension(nsoil) :: soil_interface_flux
     real, dimension(nsoil) :: theoretical_temperature
     
     iret = nf90_inq_varid(ncid, "time", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire time variable")
     iret = nf90_put_var(ncid, varid, itime,                 start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put time variable")

     iret = nf90_inq_varid(ncid, "temperature_ground", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire temperature_ground variable")
     iret = nf90_put_var(ncid, varid, temperature_ground,    start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put temperature_ground variable")

     iret = nf90_inq_varid(ncid, "ground_heat", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire ground_heat variable")
     iret = nf90_put_var(ncid, varid, ground_heat,           start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put ground_heat variable")

     iret = nf90_inq_varid(ncid, "soil_temperature_lev1", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire soil_temperature_lev1 variable")
     iret = nf90_put_var(ncid, varid,  temperature_soil(1),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put soil_temperature_lev1 variable")

     iret = nf90_inq_varid(ncid, "soil_temperature_lev2", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire soil_temperature_lev2 variable")
     iret = nf90_put_var(ncid, varid,  temperature_soil(2),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put soil_temperature_lev2 variable")

     iret = nf90_inq_varid(ncid, "soil_temperature_lev3", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire soil_temperature_lev3 variable")
     iret = nf90_put_var(ncid, varid,  temperature_soil(3),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put soil_temperature_lev3 variable")

     iret = nf90_inq_varid(ncid, "soil_temperature_lev4", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire soil_temperature_lev4 variable")
     iret = nf90_put_var(ncid, varid,  temperature_soil(4),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put soil_temperature_lev4 variable")

     iret = nf90_inq_varid(ncid, "df_lev1", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire df_lev1 variable")
     iret = nf90_put_var(ncid, varid,  thermal_conductivity(1),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put df_lev1 variable")

     iret = nf90_inq_varid(ncid, "df_lev2", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire df_lev2 variable")
     iret = nf90_put_var(ncid, varid,  thermal_conductivity(2),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put df_lev2 variable")

     iret = nf90_inq_varid(ncid, "df_lev3", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire df_lev3 variable")
     iret = nf90_put_var(ncid, varid,  thermal_conductivity(3),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put df_lev3 variable")

     iret = nf90_inq_varid(ncid, "df_lev4", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire df_lev4 variable")
     iret = nf90_put_var(ncid, varid,  thermal_conductivity(4),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put df_lev4 variable")

     iret = nf90_inq_varid(ncid, "hcpct_lev1", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire hcpct_lev1 variable")
     iret = nf90_put_var(ncid, varid,  heat_capacity(1),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put hcpct_lev1 variable")

     iret = nf90_inq_varid(ncid, "hcpct_lev2", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire hcpct_lev2 variable")
     iret = nf90_put_var(ncid, varid,  heat_capacity(2),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put hcpct_lev2 variable")

     iret = nf90_inq_varid(ncid, "hcpct_lev3", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire hcpct_lev3 variable")
     iret = nf90_put_var(ncid, varid,  heat_capacity(3),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put hcpct_lev3 variable")

     iret = nf90_inq_varid(ncid, "hcpct_lev4", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire hcpct_lev4 variable")
     iret = nf90_put_var(ncid, varid,  heat_capacity(4),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put hcpct_lev4 variable")

     iret = nf90_inq_varid(ncid, "theory_lev1", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire theory_lev1 variable")
     iret = nf90_put_var(ncid, varid,  theoretical_temperature(1),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put theory_lev1 variable")

     iret = nf90_inq_varid(ncid, "theory_lev2", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire theory_lev2 variable")
     iret = nf90_put_var(ncid, varid,  theoretical_temperature(2),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put theory_lev2 variable")

     iret = nf90_inq_varid(ncid, "theory_lev3", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire theory_lev3 variable")
     iret = nf90_put_var(ncid, varid,  theoretical_temperature(3),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put theory_lev3 variable")

     iret = nf90_inq_varid(ncid, "theory_lev4", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire theory_lev4 variable")
     iret = nf90_put_var(ncid, varid,  theoretical_temperature(4),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put theory_lev4 variable")
  
  if(itime>0) then
     iret = nf90_inq_varid(ncid, "diff_lev1", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire diff_lev1 variable")
     iret = nf90_put_var(ncid, varid,  temperature_soil(1)-theoretical_temperature(1),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put diff_lev1 variable")

     iret = nf90_inq_varid(ncid, "diff_lev2", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire diff_lev2 variable")
     iret = nf90_put_var(ncid, varid,  temperature_soil(2)-theoretical_temperature(2),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put diff_lev2 variable")

     iret = nf90_inq_varid(ncid, "diff_lev3", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire diff_lev3 variable")
     iret = nf90_put_var(ncid, varid,  temperature_soil(3)-theoretical_temperature(3),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put diff_lev3 variable")

     iret = nf90_inq_varid(ncid, "diff_lev4", varid)
      if (iret /= nf90_noerr) call handle_err(iret,"inquire diff_lev4 variable")
     iret = nf90_put_var(ncid, varid,  temperature_soil(4)-theoretical_temperature(4),     start=(/itime+1/))
      if (iret /= nf90_noerr) call handle_err(iret,"put diff_lev4 variable")
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

