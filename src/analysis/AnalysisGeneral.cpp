#include "../../include/analysis/AnalysisGeneral.h"

#include <cmath>
#include <sstream>
#include <TCanvas.h>

#include "../../include/other/HashManager.h"

AnalysisGeneral::AnalysisGeneral() : Analysis("AnalysisGeneral"){

	//set histograms
	m_hXBvsQ2 = new TH2D((HashManager::getInstance()->getHash()).c_str(), "xB vs. Q2", 
		100, -4., 0., 100, 0., 3.);
}

AnalysisGeneral::~AnalysisGeneral(){
}

void AnalysisGeneral::fill(const DVCSEvent& event, double weight){

	//fill histograms
	m_hXBvsQ2->Fill(log10(event.getXB()), log10(event.getQ2()), weight);
}

void AnalysisGeneral::analyse(){
	//nothing to be done here
}

void AnalysisGeneral::plot(const std::string& path){

	//canvases
	std::vector<TCanvas*> cans;

	//xB vs. Q2
	cans.push_back(
		new TCanvas(
			(HashManager::getInstance()->getHash()).c_str(), "")
	);

	m_hXBvsQ2->Draw("colz");
	cans.back()->SetLogz();

	//print
	for(size_t i = 0; i < cans.size(); i++){

			if(cans.size() > 1 && i == 0 ){
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

