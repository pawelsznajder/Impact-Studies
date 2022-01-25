#include <cmath>
#include <filesystem>
#include <iostream>
#include <cstdlib>
#include "HepMC3/ReaderAscii.h"

#include "../include/analysis/AnalysisALU.h"

int main(int argc, char* argv[]){

	//check arguments
	if(argc == 1){
		std::cout << __func__ << " error: usage: " << argv[0] << " dir1 dir2 ..." << std::endl;
		exit(0);
	}

	//analysis objects
	AnalysisALU analysisALU;

	//loop over directories
	for(size_t i = 1; i < argc; i++){

		//check if exists
		if(! std::filesystem::exists(std::filesystem::path(argv[i]))){
			
			std::cout << __func__ << " warning: directory: " << argv[i] << " does not exist" << std::endl;
			continue;
		}

		//loop over files
		for (const auto& dirEntry : std::filesystem::recursive_directory_iterator(argv[i])){

			//skip directories
			if(! dirEntry.is_regular_file()) continue;

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
	               	 	analysisALU.fill(dvcsEvent, 1.);
					}

					//close file
	    			inputFile.close();
    			}
			}
		}
	}
   
	//analyse
	analysisALU.analyse();

	//print
	analysisALU.plot("test.pdf");

	return 0;
}