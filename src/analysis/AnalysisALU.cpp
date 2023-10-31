#include "../../include/analysis/AnalysisALU.h"

#include <cmath>
#include <sstream>
#include <TCanvas.h>

#include "../../include/other/HashManager.h"

AnalysisALU::AnalysisALU(double targetLuminosity) : Analysis("AnalysisALU", targetLuminosity), 
	m_lumiM(0.), m_lumiP(0.){

	//set bin boundaries
	setBinBoundaries();

	//initialise bins
	initialiseBins();
}

AnalysisALU::~AnalysisALU(){
}

void AnalysisALU::fill(DVCSEvent& event, double weight){

	//only comming from DVCS+INT+BH sample
	if(! (
		event.checkSubProcessType(SubProcessType::BH) && 
		event.checkSubProcessType(SubProcessType::INT) && 
		event.checkSubProcessType(SubProcessType::DVCS)
		)
	) return;

	//only to reach target luminosity
	switch(event.getBeamPolarisation()){

   	 	case -1:{

   	 		m_lumiM += weight;

			if(m_lumiM >= m_targetLuminosity/2.){
				weight *= -1;
			}

   	 		break;
   	 	}

       	case 1:{

   	 		m_lumiP += weight;

			if(m_lumiP >= m_targetLuminosity/2.){
				weight *= -1;
			}

       	 	break;
       	}

       default:{
       		std::cout << __func__ << " error: wrong beam polarisation state, " << event.getBeamPolarisation() << std::endl;
			exit(0);
       }
	}

	//fill
	for(std::vector<BinALU>::iterator it = m_bins.begin(); 
		it != m_bins.end(); it++){

		if(event.getXB() < it->getRangeXB().first || 
			event.getXB() >= it->getRangeXB().second) continue; 

		if(event.getQ2() < it->getRangeQ2().first || 
			event.getQ2() >= it->getRangeQ2().second) continue; 

		if(fabs(event.getT()) < it->getRangeT().first || 
			fabs(event.getT()) >= it->getRangeT().second) continue; 

		it->fill(event, weight);
	}
}

void AnalysisALU::analyse(){

	//loop
	for(std::vector<BinALU>::iterator it = m_bins.begin(); 
		it != m_bins.end(); it++){

		//make analysis
		it->analyse();

		//get results
		FitResult* fitResult = it->getFitResult();

		//check if not empty
		if(fitResult == nullptr) continue;

		//print
		it->print();
		fitResult->print();
	}
}

void AnalysisALU::plot(const std::string& path){

	//canvases
	std::vector<TCanvas*> cans;

	//===============================================
	// ONE CANVAS ALL PLOTS
	//===============================================

	//loop
	for(std::vector<std::pair<double, double> >::const_iterator itT = m_binRangesT.begin(); 
		itT != m_binRangesT.end(); itT++){

		//labes
		std::stringstream ss;

		ss << itT->first << " #leq |t| < " << itT->second;

		//loop over canvases for this t bin
		for(size_t i = 0; i < 8; i++){

			//new
			cans.push_back(
				new TCanvas(
					(HashManager::getInstance()->getHash()).c_str(), ss.str().c_str())
			);

			//divide
			cans.back()->Divide(m_binRangesXB.size(), m_binRangesQ2.size());

			//plot
			for(std::vector<std::pair<double, double> >::const_iterator itXB = m_binRangesXB.begin(); 
				itXB != m_binRangesXB.end(); itXB++){
				for(std::vector<std::pair<double, double> >::const_iterator itQ2 = m_binRangesQ2.begin(); 
					itQ2 != m_binRangesQ2.end(); itQ2++){
				
					//iterator
					std::vector<BinALU>::const_iterator itBin;

					//look for bin
					for(itBin = m_bins.begin(); itBin != m_bins.end(); itBin++){
						if(
							itBin->getRangeT().first == itT->first && 
							itBin->getRangeT().second == itT->second && 
							itBin->getRangeXB().first == itXB->first && 
							itBin->getRangeXB().second == itXB->second && 
							itBin->getRangeQ2().first == itQ2->first && 
							itBin->getRangeQ2().second == itQ2->second 
						) break;
					}

					//check if found
					if(itBin == m_bins.end()){
						std::cout << getClassName() << "::" <<__func__ << " error: " << 
							"not able to find bin" << std::endl;
						exit(0);
					}

					//check if not empty
					if(itBin->getNEvents() == 0) continue;

					//set pad
					cans.back()->cd(1 + 
						size_t(itXB - m_binRangesXB.begin()) + 
						(size_t(m_binRangesQ2.end() - itQ2) - 1) * m_binRangesXB.size() 
					);

					if(i == 0){

						//histograms
						TH1* h1 = itBin->getHDistributions().first;
						TH1* h2 = itBin->getHDistributions().second;

						//set minima
						h1->SetMinimum(0.);
						h2->SetMinimum(0.);

						//colors
						h1->SetLineColor(2);
						h2->SetLineColor(4);

						//no stats
						h1->SetStats(0);
						h2->SetStats(0);

						//draw
						h1->Draw();
						h2->Draw("same");
					}

					if(i == 1){

						//histograms
						TH1* h1 = itBin->getHDistributionsBorn().first;
						TH1* h2 = itBin->getHDistributionsBorn().second;

						//set minima
						h1->SetMinimum(0.);
						h2->SetMinimum(0.);

						//colors
						h1->SetLineColor(2);
						h2->SetLineColor(4);

						//no stats
						h1->SetStats(0);
						h2->SetStats(0);

						//draw
						h1->Draw();
						h2->Draw("same");
					}

					if(i == 2){

						//histograms
						TH1* h1 = itBin->getHDistributionsRC().first;
						TH1* h2 = itBin->getHDistributionsRC().second;

						if(h1 == nullptr || h2 == nullptr) continue;

						//set minima
						h1->SetMinimum(0.);
						h2->SetMinimum(0.);

						//colors
						h1->SetLineColor(2);
						h2->SetLineColor(4);

						//no stats
						h1->SetStats(0);
						h2->SetStats(0);

						//draw
						h1->Draw();
						h2->Draw("same");
					}

					if(i == 3){

						//histograms
						TH1* h1 = itBin->getHDistributionsAcceptance().first;
						TH1* h2 = itBin->getHDistributionsAcceptance().second;

						if(h1 == nullptr || h2 == nullptr) continue;

						//set minima
						h1->SetMinimum(0.);
						h2->SetMinimum(0.);

						//colors
						h1->SetLineColor(2);
						h2->SetLineColor(4);

						//no stats
						h1->SetStats(0);
						h2->SetStats(0);

						//draw
						h1->Draw();
						h2->Draw("same");
					}

					if(i == 4){

						//histograms
						TH1* h1 = itBin->getHDistributionsCorrected().first;
						TH1* h2 = itBin->getHDistributionsCorrected().second;

						if(h1 == nullptr || h2 == nullptr) continue;

						//set minima
						h1->SetMinimum(0.);
						h2->SetMinimum(0.);

						//colors
						h1->SetLineColor(2);
						h2->SetLineColor(4);

						//no stats
						h1->SetStats(0);
						h2->SetStats(0);

						//draw
						h1->Draw();
						h2->Draw("same");
					}

					if(i == 5){

					 	//histogram
					 	TH1* h = itBin->getHSum();

					 	//check if not empty
					 	if(h != nullptr){

						 	//set minimum and maximum
						 	h->SetMinimum(0);

						 	//no stats
						 	h->SetStats(0);

						 	//draw
						 	h->Draw();
					 	}
					}

					if(i == 6){

					 	//histogram
					 	TH1* h = itBin->getHDifference();

					 	//check if not empty
					 	if(h != nullptr){

						 	//no stats
						 	h->SetStats(0);

						 	//draw
						 	h->Draw();
					 	}
					}

					if(i == 7){

					 	//histogram
					 	TH1* h = itBin->getHAsymmetry();

					 	//check if not empty
					 	if(h != nullptr){

						 	//set minimum and maximum
						 	h->SetMinimum(-1.);
						 	h->SetMaximum(1.);

						 	//no stats
						 	h->SetStats(0);

						 	//draw
						 	h->Draw();

						 	//draw associated function
						 	TObject* o = h->GetListOfFunctions()->First();

						 	if(o != 0){
						 		static_cast<TF1*>(o)->Draw("same");
						 	}
					 	}
					}
				}
			}
		}
	}

	//===============================================
	// ONE CANVAS ONE PLOT
	//===============================================

	for(std::vector<std::pair<double, double> >::const_iterator itT = m_binRangesT.begin(); 
		itT != m_binRangesT.end(); itT++){
		for(std::vector<std::pair<double, double> >::const_iterator itXB = m_binRangesXB.begin(); 
				itXB != m_binRangesXB.end(); itXB++){
				for(std::vector<std::pair<double, double> >::const_iterator itQ2 = m_binRangesQ2.begin(); 
					itQ2 != m_binRangesQ2.end(); itQ2++){

				//iterator
				std::vector<BinALU>::const_iterator itBin;

				//look for bin
				for(itBin = m_bins.begin(); itBin != m_bins.end(); itBin++){
					if(
						itBin->getRangeT().first == itT->first && 
						itBin->getRangeT().second == itT->second && 
						itBin->getRangeXB().first == itXB->first && 
						itBin->getRangeXB().second == itXB->second && 
						itBin->getRangeQ2().first == itQ2->first && 
						itBin->getRangeQ2().second == itQ2->second 
					) break;
				}

				//check if found
				if(itBin == m_bins.end()){
					std::cout << getClassName() << "::" <<__func__ << " error: " << 
						"not able to find bin" << std::endl;
					exit(0);
				}

				//check if not empty
				if(itBin->getNEvents() == 0) continue;

				//add canvas
				cans.push_back(new TCanvas(
					(HashManager::getInstance()->getHash()).c_str(), ""));

				//histograms
				TH1* h1 = itBin->getHDistributions().first;
				TH1* h2 = itBin->getHDistributions().second;

				//set minima
				h1->SetMinimum(0.);
				h2->SetMinimum(0.);

				//colors
				h1->SetLineColor(2);
				h2->SetLineColor(4);

				//no stats
				h1->SetStats(0);
				h2->SetStats(0);

				//draw
				h1->Draw();
				h2->Draw("same");

				//add canvas
				cans.push_back(new TCanvas(
					(HashManager::getInstance()->getHash()).c_str(), ""));

				//histogram
			 	TH1* h = itBin->getHAsymmetry();

			 	//check if not empty
			 	if(h != nullptr){

				 	//set minimum and maximum
				 	h->SetMinimum(-1.);
				 	h->SetMaximum(1.);

				 	//no stats
				 	h->SetStats(0);

				 	//draw
				 	h->Draw();

				 	//draw associated function
				 	TObject* o = h->GetListOfFunctions()->First();

				 	if(o != 0){
				 		static_cast<TF1*>(o)->Draw("same");
				 	}
			 	}
			}
		}
	}

	//===============================================
	// PRINT
	//===============================================

	//print
	for(size_t i = 0; i < cans.size(); i++){

		if(i == 0){
			cans[i]->Print((path+"(").c_str(), "pdf");
		}
		else if(i == cans.size() - 1){
			cans[i]->Print((path+")").c_str(), "pdf");
		}
		else{
			cans[i]->Print(path.c_str(), "pdf");
		}
	}
}

void AnalysisALU::setBinBoundaries(){

    m_binBoundsXB = {
    	0.0001, 0.000158489, 0.000251189, 0.000398107, 
    	0.000630957, 0.001, 0.00158489, 0.00251189, 
    	0.00398107, 0.00630957, 0.01, 0.0158489, 
    	0.0251189, 0.0398107, 0.0630957, 0.1, 
    	0.158489, 0.251189, 0.398107, 0.630957
    };

    m_binBoundsQ2 = {
	   	1.0, 1.77828, 3.16228, 5.62341, 
    	10.0, 17.7828, 31.6228, 56.2341, 
    	100.0, 177.828, 316.228, 562.341, 
    	1000.0
    };

    m_binBoundsT = {
    	0.2, 0.3, 0.4, 0.5, 
    	0.6, 0.7, 0.8, 0.9, 
    	1.0, 1.1, 1.2, 1.3
    };

    m_nBinsPhi = 20;
}

void AnalysisALU::initialiseBins(){

	//make ranges
	for(std::vector<double>::const_iterator itXB = (m_binBoundsXB.begin() + 1); 
			itXB != m_binBoundsXB.end(); itXB++){
		m_binRangesXB.push_back(std::make_pair(*(itXB - 1), *itXB));
	}

	for(std::vector<double>::const_iterator itQ2 = (m_binBoundsQ2.begin() + 1); 
			itQ2 != m_binBoundsQ2.end(); itQ2++){
		m_binRangesQ2.push_back(std::make_pair(*(itQ2 - 1), *itQ2));
	}

	for(std::vector<double>::const_iterator itT = (m_binBoundsT.begin() + 1); 
			itT != m_binBoundsT.end(); itT++){	
		m_binRangesT.push_back(std::make_pair(*(itT - 1), *itT));
	}

	//make bins
	for(std::vector<std::pair<double, double> >::const_iterator itXB = m_binRangesXB.begin(); 
			itXB != m_binRangesXB.end(); itXB++){
		for(std::vector<std::pair<double, double> >::const_iterator itQ2 = m_binRangesQ2.begin();
			itQ2 != m_binRangesQ2.end(); itQ2++){			
			for(std::vector<std::pair<double, double> >::const_iterator itT = m_binRangesT.begin(); 
					itT != m_binRangesT.end(); itT++){
				m_bins.push_back(
					BinALU(
						*itXB, *itQ2, *itT,
						m_nBinsPhi, std::make_pair(0., 2 * M_PI)
					)
				);
			}
		}
	}
}
