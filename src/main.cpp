#include <cmath>
#include <iostream>
#include <cstdlib>
#include <HepMC3/ReaderAscii.h>

#ifdef __APPLE__
    #include <filesystem>
    namespace fs = std::filesystem;
#endif

#ifdef __linux__
    #if __GNUC__ >= 9
            #include <filesystem>
            namespace fs = std::filesystem;
    #else
            #include <experimental/filesystem>
            namespace fs = std::experimental::filesystem;
    #endif
#endif

#include "../include/analysis/AnalysisGeneral.h"
#include "../include/analysis/AnalysisALU.h"

int main(int argc, char* argv[]){

	//check arguments
	if(argc == 1){
		std::cout << __func__ << " error: usage: " << argv[0] << " dir1 dir2 ..." << std::endl;
		exit(0);
	}

	//analysis objects
	AnalysisGeneral analysisGeneral;
	AnalysisALU analysisALU;

	//loop over directories
	for(size_t i = 1; i < argc; i++){

		//check if exists
		if(! fs::exists(fs::path(argv[i]))){
			
			std::cout << __func__ << " warning: directory: " << argv[i] << " does not exist" << std::endl;
			continue;
		}

		//loop over files
		for (const auto& dirEntry : fs::recursive_directory_iterator(argv[i])){

			//skip directories
			if(! fs::is_regular_file(dirEntry)) continue;

			//txt
			if(dirEntry.path().extension() == ".txt"){
				
				//print
				std::cout << __func__ << 
					" info: reading: " << dirEntry.path() << std::endl;

				//variables
				std::pair<double, double> crossSection;
				size_t nEvents;
				int beamPolarisation;

				//read to collect atributes ===============
				{
					HepMC3::ReaderAscii inputFile(dirEntry.path());

					//loop over events
					for(;;){

						//event
	                	HepMC3::GenEvent evt(Units::GEV,Units::MM);
	               
	               		//read
	          			inputFile.read_event(evt);

	                	//if reading failed - exit loop
	                	if(inputFile.failed() ) break;
					}

					//run info
					std::shared_ptr<HepMC3::GenRunInfo> runInfo = inputFile.run_info();

					crossSection.first = 
						std::stod((runInfo->attributes().find("integrated_cross_section_value")->second)->unparsed_string());
					crossSection.second = 
						std::stod((runInfo->attributes().find("integrated_cross_section_uncertainty")->second)->unparsed_string());
					nEvents = 
						std::stoul((runInfo->attributes().find("generated_events_number")->second)->unparsed_string());
					beamPolarisation = 
						std::stoi((runInfo->attributes().find("beam_polarisation")->second)->unparsed_string());

					//close
    				inputFile.close();
				}

    			//read to process events ===============
    			{
					HepMC3::ReaderAscii inputFile(dirEntry.path());

					//loop over events
					for(;;){

						//event
	                	HepMC3::GenEvent evt(Units::GEV,Units::MM);
	               
	               		//read
	          			inputFile.read_event(evt);

	                	//if reading failed - exit loop
	                	if(inputFile.failed() ) break;

	                	//DVCS event 
	                	//TODO add beam charge and target polarisation
	               	 	DVCSEvent dvcsEvent(evt, beamPolarisation, -1, TVector3(0., 0., 0.));

	               	 	//fill
	               	 	//TODO add weight 
	               	 	analysisGeneral.fill(dvcsEvent, 1.);
	               	 	analysisALU.fill(dvcsEvent, 1.);
					}

					//close file
	    			inputFile.close();
    			}
			}
		}
	}
   
	//analyse
	analysisGeneral.analyse();
	analysisALU.analyse();

	//print
	analysisGeneral.plot("analysisGeneral.pdf");
	analysisALU.plot("analysisALU.pdf");

	return 0;
}