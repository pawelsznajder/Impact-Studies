#include "../../include/analysis/AnalysisTSlope.h"

#include <cmath>
#include <sstream>
#include <TCanvas.h>


#include "../../include/other/HashManager.h"

AnalysisTSlope::AnalysisTSlope() : Analysis("AnalysisTSlope"){

	//set bin boundaries
	setBinBoundaries();

	//initialise bins
	initialiseBins();
}

AnalysisTSlope::~AnalysisTSlope(){
}

void AnalysisTSlope::fill(const DVCSEvent& event, double weight){

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
		FitResult fitResult = it->analyse();

		//print if only something happend 
		if(fitResult.getStatusCode() != -10){

			it->print();
			fitResult.print();
		}
	}
}

void AnalysisTSlope::plot(const std::string& path){

	//canvases
	std::vector<TCanvas*> cans[2];

	//labes
	std::stringstream ss;

	//loop over canvases for this t bin
	for(size_t i = 0; i < 2; i++){

		//new
		cans[i].push_back(
			new TCanvas(
				(HashManager::getInstance()->getHash()).c_str(), ss.str().c_str())
		);

		//divide
		cans[i].back()->Divide(m_binRangesXB.size(), m_binRangesQ2.size());

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

				//check if not empty
				if(itBin->getNEvents() == 0) continue;

				//set pad
				cans[i].back()->cd(1 + 
					size_t(itXB - m_binRangesXB.begin()) + 
					(size_t(m_binRangesQ2.end() - itQ2) - 1) * m_binRangesXB.size() 
				);

				//set log-scale
				cans[i].back()->cd(1 + 
					size_t(itXB - m_binRangesXB.begin()) + 
					(size_t(m_binRangesQ2.end() - itQ2) - 1) * m_binRangesXB.size())->SetLogy();

				if(i == 0){

					//histograms
					TH1* h = itBin->getHDistributions();

					//set minima
					h->SetMinimum(1);

					//colors
					h->SetLineColor(2);

					//no stats
					h->SetStats(0);

					//draw
					h->Draw();
				}

				if(i == 1){

					//histogram
				  	TH1* h = itBin->getHTSlope();

					 //check if not empty
					 if(h != nullptr){

						 //set minimum and maximum
						 h->SetMinimum(1.);
					 	 //h->SetMaximum(10000.);

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

	//print
	for(size_t j = 0; j < cans[0].size(); j++){
		for(size_t i = 0; i < 2; i++){

			if(cans[0].size() > 1 && i == 0 && j == 0){
				cans[i][j]->Print((path+"(").c_str(), "pdf");
			}
			else if(i == 1 && j == cans[0].size() - 1){
				cans[i][j]->Print((path+")").c_str(), "pdf");
			}
			else{
				cans[i][j]->Print(path.c_str(), "pdf");
			}
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

    m_nBinsT = 11;

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
						*itXB, *itQ2, m_nBinsT, std::make_pair(-1.3, -0.2)
					)
				);
			
		}
	}
}
