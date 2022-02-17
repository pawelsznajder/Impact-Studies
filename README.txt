To compile:

cd build
cmake ..
make

If HepMC3 is not installed in your system use:

-DHepMC3_DIR=ABSOLUTE_PATH

where ABSOLUTE_PATH is the absolute path to directory containing HepMC3Config.cmake file (typically in share/HepMC3/cmake subfolder of HepMC3 instalation)

To pick up specific HepMC3 library (against that default one in your system) use this in cmake:

-DHEPMC3_INCLUDE_DIR=ABSOLUTE_PATH_INCLUDE -DHEPMC3_LIBRARIES=ABSOLUTE_PATH_LIB

where ABSOLUTE_PATH_INCLUDE is the absolute path to include directory and where ABSOLUTE_PATH_LIB is the path to the HepMC3 libraty (typically libHepMC3.dylib or libHepMC3.so)

To run:

cd bin
./eic_impact_studies directory(ies)_with_EpIC_output_files
