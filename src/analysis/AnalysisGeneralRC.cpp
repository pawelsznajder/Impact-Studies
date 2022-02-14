#include "../../include/analysis/AnalysisGeneralRC.h"

#include <cmath>
#include <sstream>
#include <TCanvas.h>
#include <TGraph.h>

#include "../../include/other/HashManager.h"

AnalysisGeneralRC::AnalysisGeneralRC() : Analysis("AnalysisGeneralRC"){

	//number of bins 
	int nbin = 100;

	//set 2D histograms
	m_hXBvsXB = new TH2D((HashManager::getInstance()->getHash()).c_str(), 
		"xB_{Born} vs. (xB_{Born}-xB_{RC})/xB_{Born}", nbin, -4., 0., nbin, -8., 0.1);
	m_hQ2vsQ2 = new TH2D((HashManager::getInstance()->getHash()).c_str(), 
		"Q2_{Born} vs. (Q2_{Born}-Q2_{RC})/Q2_{Born}", nbin, 0., 4., nbin, -8., 2.);
	m_hTvsT = new TH2D((HashManager::getInstance()->getHash()).c_str(), 
		"|t_{Born}| vs. (|t_{Born}|-|t_{RC}|)/|t_{Born}|",nbin, 0., 2., nbin, -0.2, 0.2);
	m_hPhivsPhi = new TH2D((HashManager::getInstance()->getHash()).c_str(), 
		"#varphi_{Born} [rad] vs. (#varphi_{Born}-#varphi_{RC})/#varphi_{Born}",nbin, 0., 2.*TMath::Pi(), nbin, -10., 1.);
	m_hPhiSvsPhiS = new TH2D((HashManager::getInstance()->getHash()).c_str(), 
		"#varphi_{S}^{Born} [rad] vs. (#varphi_{s}^{Born}-#varphi_{S}^{RC})/#varphi_{S}^{Born}", nbin, 0., 2.*TMath::Pi(), nbin, -1., 1.);
	m_hYvsY = new TH2D((HashManager::getInstance()->getHash()).c_str(), 
		"y_{Born} vs. (y_{Born}-y_{RC})/y_{Born}", nbin, -0.1, 1.1, nbin, -15., 0.1);

	size_t maxNEvents = 1E4;

	m_analyserXB = new MeanSigmaAnalyser(m_hXBvsXB, maxNEvents);
	m_analyserQ2 = new MeanSigmaAnalyser(m_hQ2vsQ2, maxNEvents);
	m_analyserT = new MeanSigmaAnalyser(m_hTvsT, maxNEvents);
	m_analyserPhi = new MeanSigmaAnalyser(m_hPhivsPhi, maxNEvents);
	m_analyserPhiS = new MeanSigmaAnalyser(m_hPhiSvsPhiS, maxNEvents);
	m_analyserY = new MeanSigmaAnalyser(m_hYvsY, maxNEvents);
}

AnalysisGeneralRC::~AnalysisGeneralRC(){

	if(m_analyserXB){
		delete m_analyserXB; m_analyserXB = 0;
	}

	if(m_analyserQ2){
		delete m_analyserQ2; m_analyserQ2 = 0;
	}

	if(m_analyserT){
		delete m_analyserT; m_analyserT = 0;
	}

	if(m_analyserPhi){
		delete m_analyserPhi; m_analyserPhi = 0;
	}

	if(m_analyserPhiS){
		delete m_analyserPhiS; m_analyserPhiS = 0;
	}

	if(m_analyserY){
		delete m_analyserY; m_analyserY = 0;
	}
}

void AnalysisGeneralRC::fill(DVCSEvent& event, double weight){

	//define variables for 2D histograms
	double XBb = log10(event.getXB(RCType::Born));
	double XBrc = log10(event.getXB(RCType::ISR | RCType::FSR));
	double Q2b = log10(event.getQ2(RCType::Born));
	double Q2rc = log10(event.getQ2(RCType::ISR | RCType::FSR));
	double Tb = fabs(event.getT(RCType::Born));
	double Trc = fabs(event.getT(RCType::ISR | RCType::FSR));
	double Phib = event.getPhi(RCType::Born);
	double Phirc = event.getPhi(RCType::ISR | RCType::FSR);
	double PhiSb = event.getPhiS(RCType::Born);
	double PhiSrc = event.getPhiS(RCType::ISR | RCType::FSR);
	double Yb = event.getY(RCType::Born);
	double Yrc = event.getY(RCType::ISR | RCType::FSR);
	
	//fill 2D histograms
	m_hXBvsXB->Fill(XBb, (XBb-XBrc)/XBb, weight);
	m_hQ2vsQ2->Fill(Q2b, (Q2b-Q2rc)/Q2b, weight);
	m_hTvsT->Fill(Tb, (Tb-Trc)/Tb, weight);
	m_hPhivsPhi->Fill(Phib, (Phib-Phirc)/Phib, weight);
	m_hPhiSvsPhiS->Fill(PhiSb, (PhiSb-PhiSrc)/PhiSb, weight);
	m_hYvsY->Fill(Yb, (Yb-Yrc)/Yb, weight);

	//fill analysers
	m_analyserXB->fill(XBb, (XBb-XBrc)/XBb, weight);
	m_analyserQ2->fill(Q2b, (Q2b-Q2rc)/Q2b, weight);
	m_analyserT->fill(Tb, (Tb-Trc)/Tb, weight);
	m_analyserPhi->fill(Phib, (Phib-Phirc)/Phib, weight);
	m_analyserPhiS->fill(PhiSb, (PhiSb-PhiSrc)/PhiSb, weight);
	m_analyserY->fill(Yb, (Yb-Yrc)/Yb, weight);
}

void AnalysisGeneralRC::analyse(){
	//nothing to be done here
}

void AnalysisGeneralRC::plot(const std::string& path){

	//canvases
	std::vector<TCanvas*> cans;

	//loop over canvases
	for(size_t i = 0; i < 2; i++)	{

		//add canvas
		cans.push_back(
			new TCanvas(
				(HashManager::getInstance()->getHash()).c_str(), "")
		);

		cans.back()->Divide(3,2);

		if(i == 0){

			cans.back()->cd(1);
			cans.back()->cd(1)->SetLogz();
			m_hXBvsXB->SetStats(0);
			m_hXBvsXB->Draw("colz");

			cans.back()->cd(2);
			cans.back()->cd(2)->SetLogz();
			m_hQ2vsQ2->SetStats(0);
			m_hQ2vsQ2->Draw("colz");

			cans.back()->cd(3);
			cans.back()->cd(3)->SetLogz();
			m_hTvsT->SetStats(0);
			m_hTvsT->Draw("colz");

			cans.back()->cd(4);
			cans.back()->cd(4)->SetLogz();
			m_hPhivsPhi->SetStats(0);
			m_hPhivsPhi->Draw("colz");

			cans.back()->cd(5);
			cans.back()->cd(5)->SetLogz();
			m_hPhiSvsPhiS->SetStats(0);
			m_hPhiSvsPhiS->Draw("colz");

			cans.back()->cd(6);
			cans.back()->cd(6)->SetLogz();
			m_hYvsY->SetStats(0);
			m_hYvsY->Draw("colz");
		}

		if(i == 1){

			cans.back()->cd(1);
			drawMeanSigmaAnalyserHistogram(m_analyserXB->getH());

			cans.back()->cd(2);
			drawMeanSigmaAnalyserHistogram(m_analyserQ2->getH());

			cans.back()->cd(3);
			drawMeanSigmaAnalyserHistogram(m_analyserT->getH());

			cans.back()->cd(4);
			drawMeanSigmaAnalyserHistogram(m_analyserPhi->getH());

			cans.back()->cd(5);
			drawMeanSigmaAnalyserHistogram(m_analyserPhiS->getH());

			cans.back()->cd(6);
			drawMeanSigmaAnalyserHistogram(m_analyserY->getH());
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

void AnalysisGeneralRC::drawMeanSigmaAnalyserHistogram(TH1* h) const{

			TGraph* g = new TGraph(2);
			g->SetPoint(0, h->GetXaxis()->GetXmin(), 0.);
			g->SetPoint(1, h->GetXaxis()->GetXmax(), 0.);
			g->SetLineColor(2);

			h->SetStats(0);
			h->SetMinimum(-1.);
			h->SetMaximum(1.);
			h->Draw("e");
			g->Draw("same");
}

