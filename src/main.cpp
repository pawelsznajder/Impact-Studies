#include <cmath>
#include <iostream>
#include <cstdlib>
#include <regex>
#include <HepMC3/ReaderAscii.h>

#include "../include/analysis/AnalysisGeneral.h"
//#include "../include/analysis/AnalysisEpIC.h"
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
	const double targetIntegratedLuminosityFb = 10.;

	std::cout << __func__ << " info: target integrated luminosity: " << 
		targetIntegratedLuminosityFb << " [fb-1]" << std::endl;

	//accumulated luminosity for two beam polarisation states
	double integratedLuminosityFbALL[] = {0., 0.};
	double integratedLuminosityFbBH[] = {0., 0.};

	//analysis objects
	AnalysisGeneral analysisGeneral;
//	AnalysisEpIC analysisEpIC;
	// AnalysisGeneralRC analysisGeneralRC;
	AnalysisALU analysisALU(targetIntegratedLuminosityFb);
	AnalysisTSlope analysisTSlope(targetIntegratedLuminosityFb);

	//vectors to store info from files
	std::vector<std::pair<double, double> > crossSection;
	std::vector<size_t> nEvents;
	std::vector<int> beamPolarisation;
	std::vector<bool> isRCSample;
	std::vector<int> subProcessTypeMask;

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
					if(dirEntry->path().extension() != ".txt") continue;
				#else
					if(dirEntry.path().extension() != ".txt") continue;
				#endif

				//get path
				#ifdef __USE_BOOST__
					std::string path = dirEntry->path().string();
				#else
					std::string path = dirEntry.path().string();
				#endif

				//skip those with extension ".reconstructed.txt"
				if(std::regex_search(path, std::regex("\\.reconstructed\\.txt$"))) continue;

				//check if ".reconstructed.txt" exists
				std::string pathRec = std::regex_replace(path, std::regex("\\.txt$"), ".reconstructed.txt");
				
				if(! fs::exists(pathRec)){

					std::cout << __func__ << " error: file " << pathRec << " does not exist" << std::endl;
					exit(0);
				}

				//print
				std::cout << __func__ << 
					" info: reading: " << path << " and " << pathRec << ((loop == 0)?(" (init) "):( " (analyse) ")) << 
					"index: " << iFile << std::endl;

				//open
				HepMC3::ReaderAscii inputFile(path);
				HepMC3::ReaderAscii inputFileRec(pathRec);

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
						std::stoi((runInfo->attributes().find("lepton_polarisation")->second)->unparsed_string()));

					subProcessTypeMask.push_back(
						SubProcessType::getSubProcessTypeMaskFromStdString(
							(runInfo->attributes().find("suprocesses_type")->second)->unparsed_string())
					);

					//rc
					isRCSample.push_back(thisIsRCSample);

					//print status
					std::cout << __func__ << " info: atribute: cross-section: " << crossSection.back().first 
						<< " +/- " << crossSection.back().second << " [nb]"<< std::endl;
					std::cout << __func__ << " info: atribute: number of events: " << nEvents.back() << std::endl;
					std::cout << __func__ << " info: integrated luminosity: " << nEvents.back()/crossSection.back().first/1.E6 << " [fb-1]" << std::endl;
					std::cout << __func__ << " info: atribute: beam polarisation: " << beamPolarisation.back() << std::endl;
					std::cout << __func__ << " info: atribute: sub process mask: " << 
						((subProcessTypeMask.back() & SubProcessType::BH)?("BH "):("")) <<
						((subProcessTypeMask.back() & SubProcessType::DVCS)?("DVCS "):("")) <<
						((subProcessTypeMask.back() & SubProcessType::INT)?("INT "):("")) << std::endl;
					std::cout << __func__ << " info: atribute: RC simulation: " << ((isRCSample.back())?("yes"):("no")) << std::endl;

					//only ALL and BH
					if(!(
						(subProcessTypeMask.back() == SubProcessType::BH) || 
						(subProcessTypeMask.back() == (SubProcessType::BH | SubProcessType::INT | SubProcessType::DVCS))
						)
					){
						std::cout << __func__ << " error: only accepted files with subprocesses BH or ALL" << std::endl;
						exit(0);
					}
				}

    			//read to process events ===============
    			if(loop == 1){

    				//counter
					size_t iEvent = 0;

					//loop over events
					for(;;){

						//info
						if(iEvent%10000 == 0){

							std::cout << __func__ << " info: process event number: " << iEvent << std::endl;
							std::cout << __func__ << " info: integrated luminosity: BH+INT+DVCS (-1/+1): " << integratedLuminosityFbALL[0] << '/' << integratedLuminosityFbALL[1] << 
								" [fb-1], BH (-1/+1): " << integratedLuminosityFbBH[0] << '/' << integratedLuminosityFbBH[1] << std::endl; 
						}

						//event
	                	HepMC3::GenEvent evt(Units::GEV,Units::MM);
	                	HepMC3::GenEvent evtRec(Units::GEV,Units::MM);
	               
	               		//read
	          			inputFile.read_event(evt);
	          			inputFileRec.read_event(evtRec);

	                	//if reading failed - exit loop
	                	if(inputFile.failed()) break;
	                	if(inputFileRec.failed()) break;

	                	//DVCS event 
	                	//TODO add beam charge and target polarisation
	               	 	DVCSEvent dvcsEvent(evt, evtRec, beamPolarisation.at(iFile), -1, TVector3(0., 0., 0.), isRCSample.at(iFile), 
	               	 		subProcessTypeMask.at(iFile));
	       
	     				//fill
	               	 	analysisGeneral.fill(dvcsEvent, 1.);
		               	 	
	               	 	//luminosity
	               	 	size_t beamPolarisationState;

	               	 	switch(beamPolarisation.at(iFile)){

		               	 	case -1:{
		               	 		beamPolarisationState = 0;
		               	 		break;
		               	 	}

			               	case 1:{
			               		beamPolarisationState = 1;
			               	 	break;
			               	}

			               default:{
			               		std::cout << __func__ << " error: wrong beam polarisation state, " << beamPolarisation.at(iFile) << std::endl;
								exit(0);
			               }
	               	 	}

						if(subProcessTypeMask.at(iFile) == (SubProcessType::BH | SubProcessType::INT | SubProcessType::DVCS)){
							integratedLuminosityFbALL[beamPolarisationState] += 1/crossSection.at(iFile).first/1.E6;
						}

						if(subProcessTypeMask.at(iFile) == SubProcessType::BH){ 
							integratedLuminosityFbBH[beamPolarisationState] += 1/crossSection.at(iFile).first/1.E6;
						}

						//weight
						double weight = 1/crossSection.at(iFile).first/1.E6;

						//fill
						analysisALU.fill(dvcsEvent, weight);
						analysisTSlope.fill(dvcsEvent, weight);

						//counter
	               	 	iEvent++;
					}
    			}

    			//file index
    			iFile++;

				//close
				inputFile.close();
				inputFileRec.close();
			}
		}
	}

	// return 0;
   
	//analyse
	analysisGeneral.analyse();
	//analysisEpIC.analyse();
	// analysisGeneralRC.analyse();
	analysisALU.analyse();
	analysisTSlope.analyse();

	//print
	analysisGeneral.plot("analysisGeneral.pdf");
	//analysisEpIC.plot("analysisEpIC.pdf");
	// analysisGeneralRC.plot("analysisGeneralRC.pdf");
	analysisALU.plot("analysisALU.pdf");
	analysisTSlope.plot("analysisTSlope.pdf");

	std::cout << __func__ << ": info: target luminosity: " << 
		((integratedLuminosityFbALL[0] >= targetIntegratedLuminosityFb && integratedLuminosityFbALL[1] >= targetIntegratedLuminosityFb &&
			integratedLuminosityFbBH[0] >= targetIntegratedLuminosityFb && integratedLuminosityFbBH[1] >= targetIntegratedLuminosityFb)?
		(" reached"):(" NOT REACHED!!!")) << std::endl;

	return 0;
}
