FC      = gfortran
FCFLAGS = --coverage -Dstandalone -fprofile-arcs -ftest-coverage -fdefault-real-8 -fdefault-double-8 -g -fbacktrace

SRC = machine.F90 ccpp_kind_types.F90 module_state_description.F90 module_model_constants.F90 ../src/WRF/module_sf_mynnsfc_common.F90 ../src/module_sf_mynnsfc_land.F90 ../src/module_sf_mynnsfc_water.F90 ../src/module_sf_mynnsfc_ice.F90 ../src/module_sf_mynnsfc_driver.F90

SRC_CCPP = $(SRC)  module_sf_mynnsfc_ccpp_tests.F90 driver.F90
SRC_WRF  = $(SRC)  module_sf_mynnsfc_wrf_tests.F90 driver.F90

all: driver_ccpp driver_wrf

driver_ccpp:
	mkdir -p build_ccpp
	cd build_ccpp && $(FC) $(FCFLAGS) -c $(addprefix ../, $(SRC_CCPP)) -lgcov
	$(FC) $(FCFLAGS) -o $@ build_ccpp/*.o -lgcov

driver_wrf:
	mkdir -p build_wrf
	cd build_wrf && $(FC) $(FCFLAGS) -DEM_CORE=1 -c $(addprefix ../, $(SRC_WRF)) -lgcov
	$(FC) $(FCFLAGS) -o $@ build_wrf/*.o -lgcov


# run wrf and ccpp ci tests
test: all
	./driver_ccpp
	./driver_wrf

clean:
	@rm -rf build_ccpp build_wrf *.mod *.o driver_ccpp driver_wrf *.gcda *.gcno *.gcov
