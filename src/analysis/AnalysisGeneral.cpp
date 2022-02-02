#include "../../include/analysis/AnalysisGeneral.h"

#include <cmath>
#include <sstream>
#include <TCanvas.h>

#include "../../include/other/HashManager.h"

AnalysisGeneral::AnalysisGeneral() : Analysis("AnalysisGeneral"){

	int nbin = 100;
	//set 2D histograms
	m_hXBvsQ2 = new TH2D((HashManager::getInstance()->getHash()).c_str(), "xB vs. Q2", 
		100, -4., 0., 100, 0., 3.);
	//set 1D histograms
	m_hXB = new TH1D((HashManager::getInstance()->getHash()).c_str(), "x_{B}",nbin, -4., 0.);
	m_hQ2 = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Q^{2}",nbin, 0., 3.);
	m_hT = new TH1D((HashManager::getInstance()->getHash()).c_str(), "t",nbin, -2., 0.);
	m_hPhi = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#varphi [rad]",nbin, 0., 2.*TMath::Pi());
	m_hPhiS = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#varphi_{S} [rad]",nbin, 0., 2.*TMath::Pi());
	m_hY = new TH1D((HashManager::getInstance()->getHash()).c_str(), "y",	nbin, -0.1, 1.1);
	m_hEtaEOut = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta e'",nbin, -4., 2.);
	m_hEtaPOut = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta p'",nbin, 0., 10.);
	m_hEtaGOut = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta #gamma'",nbin, -10., 4.);
}

AnalysisGeneral::~AnalysisGeneral(){
}

void AnalysisGeneral::fill(const DVCSEvent& event, double weight){

	//fill 2D histograms
	m_hXBvsQ2->Fill(log10(event.getXB()), log10(event.getQ2()), weight);
	//fill 1D histograms
	m_hXB->Fill(log10(event.getXB()), weight);
	m_hQ2->Fill(log10(event.getQ2()), weight);
	m_hT->Fill(event.getT(), weight);
	m_hPhi->Fill(event.getPhi(), weight);
	m_hPhiS->Fill(event.getPhiS(), weight);
	m_hY->Fill(event.getY(), weight);
	m_hEtaEOut->Fill(event.getEtaEOut(), weight);
	m_hEtaPOut->Fill(event.getEtaPOut(), weight);
	m_hEtaGOut->Fill(event.getEtaGOut(), weight);
}

void AnalysisGeneral::analyse(){
	//nothing to be done here
}

void AnalysisGeneral::plot(const std::string& path){

	//canvases
	std::vector<TCanvas*> cans[2];

	//xB vs. Q2
	for(size_t i = 0; i < 2; i++){
		cans[i].push_back(
			new TCanvas(
				(HashManager::getInstance()->getHash()).c_str(), "")
		);
		if (i == 0){
			cans[i].back()->cd();
			m_hXBvsQ2->Draw("colz");
			cans[i].back()->SetLogz();
		} if (i == 1){
			cans[i].back()->Divide(3,3);
			cans[i].back()->cd(1);
			m_hXB->Draw();
			cans[i].back()->cd(2);
			m_hQ2->Draw();
			cans[i].back()->cd(3);
			m_hT->Draw();
			cans[i].back()->cd(4);
			m_hPhi->Draw();
			cans[i].back()->cd(5);
			m_hPhiS->Draw();
			cans[i].back()->cd(6);
			m_hY->Draw();
			cans[i].back()->cd(7);
			m_hEtaEOut->Draw();
			cans[i].back()->cd(8);
			m_hEtaPOut->Draw();
			cans[i].back()->cd(9);
			m_hEtaGOut->Draw();
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

