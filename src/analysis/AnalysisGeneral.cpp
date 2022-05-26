#include "../../include/analysis/AnalysisGeneral.h"

#include <cmath>
#include <sstream>
#include <TCanvas.h>

#include "../../include/other/HashManager.h"

AnalysisGeneral::AnalysisGeneral() : Analysis("AnalysisGeneral"){

	//number of bins 
	int nbin = 100;

	//set 2D histograms
	m_hXBvsQ2 = new TH2D((HashManager::getInstance()->getHash()).c_str(), "xB vs. Q2", 
		100, -4., 0., 100, 0., 3.);

	//set 1D histograms
	for(size_t i = 0; i < 2; i++){

		m_hXB[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "x_{B}",nbin, -4., 0.);
		m_hQ2[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Q^{2}",nbin, 0., 3.);
		m_hT[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "|t|",nbin, 0., 2.);
		m_hPhi[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#varphi [rad]",nbin, 0., 2.*TMath::Pi());
		m_hPhiS[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#varphi_{S} [rad]",nbin, 0., 2.*TMath::Pi());
		m_hY[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "y",	nbin, -0.1, 1.1);
		m_hEtaEOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta e'",nbin, -4., 2.);
		m_hEtaPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta p'",nbin, 0., 10.);
		m_hEtaGOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta #gamma'",nbin, -10., 4.);
	}
}

AnalysisGeneral::~AnalysisGeneral(){
}

void AnalysisGeneral::fill(DVCSEvent& event, double weight){

	//fill 2D histograms
	m_hXBvsQ2->Fill(log10(event.getXB()), log10(event.getQ2()), weight);

	//fill 1D histograms
	for(size_t i = 0; i < 2; i++){

		//type
		KinematicsType::Type kinematicsType = (i == 0)?(KinematicsType::True):(KinematicsType::Observed);

		//check if reconstricted
		if(kinematicsType == KinematicsType::Observed && (! event.isReconstructed())) continue;

		m_hXB[i]->Fill(log10(event.getXB(kinematicsType)), weight);
		m_hQ2[i]->Fill(log10(event.getQ2(kinematicsType)), weight);
		m_hT[i]->Fill(fabs(event.getT(kinematicsType)), weight);
		m_hPhi[i]->Fill(event.getPhi(kinematicsType), weight);
		m_hPhiS[i]->Fill(event.getPhiS(kinematicsType), weight);
		m_hY[i]->Fill(event.getY(kinematicsType), weight);
		m_hEtaEOut[i]->Fill(event.getEtaEOut(kinematicsType), weight);
		m_hEtaPOut[i]->Fill(event.getEtaPOut(kinematicsType), weight);
		m_hEtaGOut[i]->Fill(event.getEtaGOut(kinematicsType), weight);
	}
}

void AnalysisGeneral::analyse(){
	//nothing to be done here
}

void AnalysisGeneral::plot(const std::string& path){

	//canvases
	std::vector<TCanvas*> cans;

	//loop over canvases
	for(size_t i = 0; i < 2; i++){

		//add canvas
		cans.push_back(
			new TCanvas(
				(HashManager::getInstance()->getHash()).c_str(), "")
		);

		if (i == 0){

			cans.back()->cd();
			m_hXBvsQ2->Draw("colz");
			cans.back()->SetLogz();

		} if (i == 1){

			cans.back()->Divide(3,3);

			std::string options;

			for(size_t j = 0; j < 2; j++){

				if(j == 0){
					options = "";
				}else{
					options = "same";
				}

				cans.back()->cd(1);
				m_hXB[j]->Draw(options.c_str());

				cans.back()->cd(2);
				cans.back()->cd(2)->SetLogy();
				m_hQ2[j]->Draw(options.c_str());

				cans.back()->cd(3);
				cans.back()->cd(3)->SetLogy();
				m_hT[j]->Draw(options.c_str());

				cans.back()->cd(4);
				m_hPhi[j]->Draw(options.c_str());

				cans.back()->cd(5);
				m_hPhiS[j]->Draw(options.c_str());

				cans.back()->cd(6);
				m_hY[j]->Draw(options.c_str());

				cans.back()->cd(7);
				m_hEtaEOut[j]->Draw(options.c_str());

				cans.back()->cd(8);
				m_hEtaPOut[j]->Draw(options.c_str());

				cans.back()->cd(9);
				m_hEtaGOut[j]->Draw(options.c_str());
			}
		}
	}

	//print
	for(size_t i = 0; i < cans.size(); i++){

		if(cans.size() > 1 && i == 0){
			cans[i]->Print((path+"(").c_str(), "pdf");
		}
		else if(cans.size() > 1 && i == cans.size() - 1){
			cans[i]->Print((path+")").c_str(), "pdf");
		}
		else{
			cans[i]->Print(path.c_str(), "pdf");
		}
	}
}

