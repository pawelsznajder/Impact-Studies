To compile:

cd build
cmake ..
make

To pick up specific HepMC3 library use this in cmake:

-DHepMC3_DIR=ABSOLUTE_PATH

To run:

cd bin
./eic_impact_studies directory(ies)_with_EpIC_output_files
