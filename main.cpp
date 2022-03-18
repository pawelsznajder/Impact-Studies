#include <cmath>
#include <iostream>
#include <cstdlib>
#include <HepMC3/ReaderAscii.h>

#include "../include/analysis/AnalysisGeneral.h"
#include "../include/analysis/AnalysisEpIC.h"
#include "../include/analysis/AnalysisGeneralRC.h"
#include "../include/analysis/AnalysisALU.h"
#include "../include/analysis/AnalysisTSlope.h"
#include "../include/other/SubProcessType.h"

#ifdef __APPLE__
    #include <filesystem>
    namespace fs = std::filesystem;
#endif

#ifdef __linux__
    #if __GNUC__ < 9
    	    #if __GNUC__ < 6

    			#define __USE_BOOST__

		 		#include <boost/filesystem.hpp>
           	 	namespace fs = boost::filesystem; 
	    #else
            	 #include <experimental/filesystem>
            	 namespace fs = std::experimental::filesystem;
    	    #endif
    #else
            #include <filesystem>
            namespace fs = std::filesystem;
    #endif
#endif

int main(int argc, char* argv[]){

	//check arguments
	if(argc == 1){
		std::cout << __func__ << " error: usage: " << argv[0] << " dir1 dir2 ..." << std::endl;
		exit(0);
	}

	//target luminosity (in fb-1)
	const double targetIntegratedLuminosityFb = 10;

	//target luminosity (in nb-1)
	const double targetIntegratedLuminosityNb = targetIntegratedLuminosityFb * 1E6;

	std::cout << __func__ << " info: target integrated luminosity: " << 
		targetIntegratedLuminosityFb << " [fb-1]" << std::endl;

	//analysis objects
	AnalysisGeneral analysisGeneral;
	AnalysisEpIC analysisEpIC;
	AnalysisGeneralRC analysisGeneralRC;
	AnalysisALU analysisALU;
	AnalysisTSlope analysisTSlope;

	//vectors to store info from files
	std::vector<std::pair<double, double> > crossSection;
	std::vector<size_t> nEvents;
	std::vector<int> beamPolarisation;
	std::vector<bool> isRCSample;
	std::vector<int> subProcessTypeMask;

	double integratedLumiWithDVCS = 0.;
	double integratedLumiWithoutDVCS = 0.;

	double totalEventsWithDVCS = 0.;
	double totalEventsWithoutDVCS = 0.;

	//loop (read each files two times)
	for(size_t loop = 0; loop < 2; loop++){

		//file index
		size_t iFile = 0;

		//loop over directories
		for(size_t i = 1; i < argc; i++){

			//check if exists
			if(! fs::exists(fs::path(argv[i]))){
				
				std::cout << __func__ << " warning: directory: " << argv[i] << " does not exist" << std::endl;
				continue;
			}

			//loop over files
			#ifdef __USE_BOOST__
				const fs::recursive_directory_iterator end;
				for(fs::recursive_directory_iterator dirEntry(fs::path(argv[i])); dirEntry != end; dirEntry++){
			#else
				for(const auto& dirEntry : fs::recursive_directory_iterator(fs::path(argv[i]))){
			#endif
			
				//skip directories
				#ifdef __USE_BOOST__
					 if(! fs::is_regular_file(dirEntry->status())) continue;
				#else
					 if(! fs::is_regular_file(dirEntry.status())) continue;
				#endif

				//txt
				#ifdef __USE_BOOST__
					if(dirEntry->path().extension() == ".txt"){
				#else
					if(dirEntry.path().extension() == ".txt"){
				#endif

					//get path
					#ifdef __USE_BOOST__
						std::string path = dirEntry->path().string();
					#else
						std::string path = dirEntry.path().string();
					#endif
					
					//print
					std::cout << __func__ << 
						" info: reading: " << path << ((loop == 0)?(" (init) "):( " (analyse) ")) << 
						"index: " << iFile << std::endl;

					//open
					HepMC3::ReaderAscii inputFile(path);

					//read to collect atributes ===============
					if(loop == 0){

						//to check if RC sample
						size_t lastParticleSize = 0;
						bool thisIsRCSample = false;

						//loop over events
						for(;;){

							//event
		                	HepMC3::GenEvent evt(Units::GEV,Units::MM);
		               
		               		//read
		          			inputFile.read_event(evt);

		          			//if the number of particles is not fixed, we have RC sample
		          			if(evt.particles().size() != 0 && evt.particles().size() != lastParticleSize){

		          				if(lastParticleSize == 0){
		          					lastParticleSize = evt.particles().size();
		          				}else{
		          					thisIsRCSample = true;
		          				}
		          			}

		                	//if reading failed - exit loop
		                	if(inputFile.failed() ) break;
						}

						//run info
						std::shared_ptr<HepMC3::GenRunInfo> runInfo = inputFile.run_info();

						crossSection.push_back(
							std::make_pair( 
								std::stod((runInfo->attributes().find("integrated_cross_section_value")->second)->unparsed_string()),
								std::stod((runInfo->attributes().find("integrated_cross_section_uncertainty")->second)->unparsed_string())
							)
						);

						nEvents.push_back(
							std::stoul((runInfo->attributes().find("generated_events_number")->second)->unparsed_string()));

						beamPolarisation.push_back(
							std::stoi((runInfo->attributes().find("beam_polarisation")->second)->unparsed_string()));

						subProcessTypeMask.push_back(
							SubProcessType::getSubProcessTypeMaskFromStdString(
								(runInfo->attributes().find("suprocesses_type")->second)->unparsed_string())
						);

						//rc
						isRCSample.push_back(thisIsRCSample);

						//luminosity
						if(subProcessTypeMask.back() & SubProcessType::DVCS){
							totalEventsWithDVCS += nEvents.back();
							//integratedLumiWithDVCS += nEvents.back() / crossSection.back().first;
						}else{
							totalEventsWithoutDVCS += nEvents.back();
							//integratedLumiWithoutDVCS += nEvents.back() / crossSection.back().first;
						}


						//print status
						std::cout << __func__ << " info: atribute: cross-section: " << crossSection.back().first 
							<< " +/- " << crossSection.back().second << std::endl;
						std::cout << __func__ << " info: atribute: number of events: " << nEvents.back() << std::endl;
						std::cout << __func__ << " info: atribute: beam polarisation: " << beamPolarisation.back() << std::endl;
						std::cout << __func__ << " info: atribute: sub process mask: " << 
							((subProcessTypeMask.back() & SubProcessType::BH)?("BH "):("")) <<
							((subProcessTypeMask.back() & SubProcessType::DVCS)?("DVCS "):("")) <<
							((subProcessTypeMask.back() & SubProcessType::INT)?("INT "):("")) << std::endl;
						std::cout << __func__ << " info: atribute: RC simulation: " << ((isRCSample.back())?("yes"):("no")) << std::endl;
					}

	    			//read to process events ===============
	    			if(loop == 1){

	    				//weight
	    				double thisWeight;

	    				if(subProcessTypeMask.at(iFile) & SubProcessType::DVCS){
	    					//thisWeight = targetIntegratedLuminosityNb / integratedLumiWithDVCS;
						thisWeight = targetIntegratedLuminosityNb * crossSection.back().first / totalEventsWithDVCS;
						}else{
							//thisWeight = targetIntegratedLuminosityNb / integratedLumiWithoutDVCS;
							thisWeight = targetIntegratedLuminosityNb * crossSection.back().first / totalEventsWithoutDVCS;
						}

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
		               	 	DVCSEvent dvcsEvent(evt, beamPolarisation.at(iFile), -1, TVector3(0., 0., 0.), isRCSample.at(iFile), 
		               	 		subProcessTypeMask.at(iFile));

		               	 	//fill
		               	 	analysisGeneral.fill(dvcsEvent, 1.);
					analysisEpIC.fill(dvcsEvent, 1.);
		               	 	analysisGeneralRC.fill(dvcsEvent, 1.);
		               	 	analysisALU.fill(dvcsEvent, 1.);
							analysisTSlope.fill(dvcsEvent, thisWeight);
						}
	    			}

	    			//file index
	    			iFile++;

					//close
    				inputFile.close();
				}
			}
		}
	}
   
	//analyse
	analysisGeneral.analyse();
	analysisEpIC.analyse();
	analysisGeneralRC.analyse();
	analysisALU.analyse();
	analysisTSlope.analyse();

	//print
	analysisGeneral.plot("analysisGeneral.pdf");
	analysisEpIC.plot("analysisEpIC.pdf");
	analysisGeneralRC.plot("analysisGeneralRC.pdf");
	analysisALU.plot("analysisALU.pdf");
	analysisTSlope.plot("analysisTSlope.pdf");

	return 0;
}
