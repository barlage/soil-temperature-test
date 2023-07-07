# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f90

include ../user_build_options

OBJS = \
	../src/soil_temperature_routines.o \
	../driver/soil_temperature_output.o \
	../driver/soil_temperature_driver.o

CMD = soil_temperature.exe
all:	$(CMD)

soil_temperature.exe: $(OBJS)
	$(COMPILERF90) -o $(@) -I../driver -I../src $(OBJS) $(NETCDFLIB)

# This command cleans up
clean:
	rm -f $(CMD)

