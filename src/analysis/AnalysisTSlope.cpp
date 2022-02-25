#include "../../include/analysis/AnalysisTSlope.h"

#include <cmath>
#include <sstream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>


#include "../../include/other/HashManager.h"
#include "../../include/other/SubProcessType.h"

AnalysisTSlope::AnalysisTSlope() : Analysis("AnalysisTSlope"){

	//set bin boundaries
	setBinBoundaries();

	//initialise bins
	initialiseBins();
}

AnalysisTSlope::~AnalysisTSlope(){
}

void AnalysisTSlope::fill(DVCSEvent& event, double weight){

	if(event.checkSubProcessType(SubProcessType::BH) && (! event.checkSubProcessType(SubProcessType::DVCS))) {
		weight *= -1.;
	}

	for(std::vector<BinTSlope>::iterator it = m_bins.begin(); 
		it != m_bins.end(); it++){

		if(event.getXB() < it->getRangeXB().first || 
			event.getXB() >= it->getRangeXB().second) continue; 

		if(event.getQ2() < it->getRangeQ2().first || 
			event.getQ2() >= it->getRangeQ2().second) continue; 

		it->fill(event, weight);
	}
}

void AnalysisTSlope::analyse(){

	//loop
	for(std::vector<BinTSlope>::iterator it = m_bins.begin(); 
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

void AnalysisTSlope::plot(const std::string& path){

	//canvases
	std::vector<TCanvas*> cans;

	//labes
	std::stringstream ss;

	//===============================================
	// ONE CANVAS ALL PLOTS
	//===============================================

	//loop over canvases for plotting xB vs. Q2 grid
	for(size_t i = 0; i < 1; i++){

		//new
		cans.push_back(new TCanvas(
				(HashManager::getInstance()->getHash()).c_str(), ss.str().c_str()));

		//divide
		cans.back()->Divide(m_binRangesXB.size(), m_binRangesQ2.size());

		//plot
		for(std::vector<std::pair<double, double> >::const_iterator itXB = m_binRangesXB.begin(); 
			itXB != m_binRangesXB.end(); itXB++){
			for(std::vector<std::pair<double, double> >::const_iterator itQ2 = m_binRangesQ2.begin(); 
				itQ2 != m_binRangesQ2.end(); itQ2++){
				
				//iterator
				std::vector<BinTSlope>::const_iterator itBin;

				//look for bin
				for(itBin = m_bins.begin(); itBin != m_bins.end(); itBin++){
					if( 
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

				//set pad
				cans.back()->cd(1 + 
					size_t(itXB - m_binRangesXB.begin()) + 
					(size_t(m_binRangesQ2.end() - itQ2) - 1) * m_binRangesXB.size() 
				);

				//set log-scale
				cans.back()->cd(1 + 
					size_t(itXB - m_binRangesXB.begin()) + 
					(size_t(m_binRangesQ2.end() - itQ2) - 1) * m_binRangesXB.size())->SetLogy();

				if(i == 0){

					//histogram
				  	TH1* h = itBin->getHTSlope();

					 //check if not empty
					 if(h != nullptr){

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

	//===============================================
	// ONE CANVAS ONE PLOT
	//===============================================

	//loop over canvases for individual xB vs. Q2 bins
	for(std::vector<std::pair<double, double> >::const_iterator itXB = m_binRangesXB.begin(); 
		itXB != m_binRangesXB.end(); itXB++){
		for(std::vector<std::pair<double, double> >::const_iterator itQ2 = m_binRangesQ2.begin(); 
			itQ2 != m_binRangesQ2.end(); itQ2++){
			
			//iterator
			std::vector<BinTSlope>::const_iterator itBin;

			//look for bin
			for(itBin = m_bins.begin(); itBin != m_bins.end(); itBin++){
				if( 
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

			//histogram
		  	TH1* h = itBin->getHTSlope();

			 //check if not empty
			 if(h != nullptr){

			 	//add canvas
				cans.push_back(new TCanvas(
					(HashManager::getInstance()->getHash()).c_str(), ""));

				//set log-scale
				cans.back()->SetLogy();

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

	//===============================================
	// ONE CANVAS ALL SLOPES
	//===============================================

	std::vector<std::vector<double> > dataSlopeX;
	std::vector<std::vector<double> > dataSlopeY;
	std::vector<std::vector<double> > dataSlopeErrX;
	std::vector<std::vector<double> > dataSlopeErrY;
	std::vector<string> dataLabels;

	auto legend = new TLegend(0.1,0.7,0.48,0.9);

	//loop over canvases for individual xB vs. Q2 bins
	for(std::vector<std::pair<double, double> >::const_iterator itQ2 = m_binRangesQ2.begin(); 
			itQ2 != m_binRangesQ2.end(); itQ2++){

		//add
		dataSlopeX.push_back(std::vector<double>());
		dataSlopeY.push_back(std::vector<double>());
		dataSlopeErrX.push_back(std::vector<double>());
		dataSlopeErrY.push_back(std::vector<double>());

		ss.str(std::string());
		ss << itQ2->first << " < Q2 < " << itQ2->second;
		dataLabels.push_back(ss.str());

		for(std::vector<std::pair<double, double> >::const_iterator itXB = m_binRangesXB.begin(); 
			itXB != m_binRangesXB.end(); itXB++){
		
			//iterator
			std::vector<BinTSlope>::const_iterator itBin;

			//look for bin
			for(itBin = m_bins.begin(); itBin != m_bins.end(); itBin++){
				if( 
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

			//get result
			FitResult* fitResult = itBin->getFitResult();

			//check if not empty
			if(fitResult == nullptr) continue;

			//add result
			dataSlopeX.back().push_back(itBin->getMeanXB());
			dataSlopeY.back().push_back(fitResult->getParameter(1).first);	//slope value
			dataSlopeErrX.back().push_back(0.);
			dataSlopeErrY.back().push_back(fitResult->getParameter(1).second);	//slope unc.
		}
	}

	//new
	cans.push_back(new TCanvas(
			(HashManager::getInstance()->getHash()).c_str(), ss.str().c_str()));

	//set log-scale
	cans.back()->SetLogx();

	//plot empty histogram that sets ranges
	TH1* hForGraphs = new TH1D((HashManager::getInstance()->getHash()).c_str(), "", 10, 0.0001, 1.);

 	//no stats
 	hForGraphs->SetStats(0);

 	//range
	hForGraphs->SetMinimum(0.);
	hForGraphs->SetMaximum(10.);

	//draw
	hForGraphs->Draw();

	//loop
	for(size_t i = 0; i < dataSlopeX.size(); i++){

		//check if any data
		if(dataSlopeX.at(i).size() == 0) continue;

		//create
		TGraphErrors* gr = new TGraphErrors(dataSlopeX.at(i).size(), &(dataSlopeX.at(i)[0]), 
			&(dataSlopeY.at(i)[0]), &(dataSlopeErrX.at(i)[0]), &(dataSlopeErrY.at(i)[0]));

		//draw
		gr->SetMarkerStyle(20);
		gr->SetMarkerColor(1+i);
		gr->Draw("same LP");

		legend->AddEntry(gr, dataLabels.at(i).c_str(),"lep");
	}
	
	legend->Draw();

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

void AnalysisTSlope::setBinBoundaries(){

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

    m_nBinsT = 20;
}

void AnalysisTSlope::initialiseBins(){

	//make ranges
	for(std::vector<double>::const_iterator itXB = (m_binBoundsXB.begin() + 1); 
			itXB != m_binBoundsXB.end(); itXB++){
		m_binRangesXB.push_back(std::make_pair(*(itXB - 1), *itXB));
	}

	for(std::vector<double>::const_iterator itQ2 = (m_binBoundsQ2.begin() + 1); 
			itQ2 != m_binBoundsQ2.end(); itQ2++){
		m_binRangesQ2.push_back(std::make_pair(*(itQ2 - 1), *itQ2));
	}

	//make bins
	for(std::vector<std::pair<double, double> >::const_iterator itXB = m_binRangesXB.begin(); 
			itXB != m_binRangesXB.end(); itXB++){
		for(std::vector<std::pair<double, double> >::const_iterator itQ2 = m_binRangesQ2.begin();
			itQ2 != m_binRangesQ2.end(); itQ2++){			
			
				m_bins.push_back(
					BinTSlope(
						*itXB, *itQ2, m_nBinsT, std::make_pair(0., 1.5)
					)
				);
		}
	}
}
