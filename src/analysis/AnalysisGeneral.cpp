#include "../../include/analysis/AnalysisGeneral.h"

#include <cmath>
#include <sstream>
#include <TCanvas.h>
#include <TLegend.h>

#include "../../include/other/HashManager.h"

AnalysisGeneral::AnalysisGeneral(double targetLuminosity) : Analysis("AnalysisGeneral", targetLuminosity), 
	m_lumi(0.){

	//number of bins 
	int nbin = 100;

	//reconstruction probabilities
	for(size_t i = 0; i < 2; i++){

		m_resProbPOut[i] = 0;
		m_resProbEOut[i] = 0;
		m_resProbGOut[i] = 0;

		m_resProbExcl[i] = 0;
	}

	//set 2D histograms
	m_hXBvsQ2 = new TH2D((HashManager::getInstance()->getHash()).c_str(), "xB vs. Q2", 
		nbin, -4., 0., nbin, 0., 3.);

	//set 1D and 2D histograms
	for(size_t i = 0; i < 2; i++){
			
		std::string title;
	
		if(i == 0){
			title = " (Generated)";
		}else{
			title = " (Full reconstructed)";
		}
	
		m_hxBvsQ2[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, -3.5, -0.2, nbin, 0., 2.2);
		m_hXBvsT[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, -3.5, -0.2, nbin, 0., 1.6);
		m_hXBvsY[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, -3.5, -0.2, nbin, 0., 0.7);
		m_hQ2vsT[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, 0., 2.2, nbin, 0., 1.6);
		m_hYvsT[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, 0., 0.7, nbin, 0., 1.6);
		m_hQ2vsY[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, 0., 2.2, nbin, 0., 0.7);

		m_hXB[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "log(x_{B})",nbin, -4., 0.);
		m_hQ2[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "log(Q^{2})",nbin, 0., 3.);
		m_hT[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "|t|",nbin, 0., 2.);
		m_hPhi[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#varphi [rad]",nbin, 0., 2.*TMath::Pi());
		m_hPhiS[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#varphi_{S} [rad]",nbin, 0., 2.*TMath::Pi());
		m_hY[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "y",	nbin, -0.1, 1.1);
		m_hEtaEOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta e'",nbin, -4., 2.);
		m_hEtaPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta p'",nbin, 0., 10.);
		m_hEtaGOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta #gamma'",nbin, -10., 4.);
		m_hPPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "momentum p' [GeV/c]",nbin, 50., 110.);
		m_hPPtOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "transverse momentum p' [GeV/c]",nbin, 0., 1.4);
		m_hPThOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Polar angle #theta p' [mrad]",nbin, 0., 25.);
		m_hPPhOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Azimuthal angle #phi p' [rad]",nbin, -TMath::Pi()-0.2, TMath::Pi()+0.2);
		m_hEPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "momentum e' [GeV/c]",nbin, 0., 20.);
		m_hEPtOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "transverse momentum e' [GeV/c]",nbin, 0.,10.);
		m_hEThOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Polar angle #theta e' [#murad]",nbin, 0., 5.);
		m_hEPhOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Azimuthal angle #phi e' [rad]",nbin, -TMath::Pi()-0.2, TMath::Pi()+0.2);
		m_hGPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "momentum #gamma' [GeV/c]",nbin, 0., 40.);
		m_hGPtOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "transverse momentum #gamma' [GeV/c]",nbin, 0., 10.);
		m_hGThOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Polar angle #theta #gamma' [#murad]",nbin, 0., 5.);
		m_hGPhOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Azimuthal angle #phi #gamma' [rad]",nbin, -TMath::Pi()-0.2, TMath::Pi()+0.2);
	}
}

AnalysisGeneral::~AnalysisGeneral(){
}

TH1* AnalysisGeneral::evaluateAcceptance(TH1** h) const{

	TH1* acc = (TH1*)h[1]->Clone();
	acc->SetName(("acc_"+std::string(h[1]->GetName())).c_str());
	acc->Divide(h[0]);
	return acc;
}

void AnalysisGeneral::fill(DVCSEvent& event, double weight){

	//only comming from DVCS+INT+BH sample
	if(! (
		event.checkSubProcessType(SubProcessType::BH) && 
		event.checkSubProcessType(SubProcessType::INT) && 
		event.checkSubProcessType(SubProcessType::DVCS)
		)
	) return;

	//luminosity
	if(m_lumi >= m_targetLuminosity) return;
	m_lumi += weight;

	//fill 1D histograms
	for(size_t i = 0; i < 2; i++){

		//type
		KinematicsType::Type kinematicsType = (i == 0)?(KinematicsType::True):(KinematicsType::Observed);

		if (event.getPOut(kinematicsType).E() > 0.) {

			m_resProbPOut[i]++;

			m_hPPOut[i]->Fill(event.getPOut(kinematicsType).P());
			m_hPPtOut[i]->Fill(event.getPOut(kinematicsType).Pt());
			m_hPThOut[i]->Fill(event.getPOut(kinematicsType).Theta()*1000.);
			m_hPPhOut[i]->Fill(event.getPOut(kinematicsType).Phi());
			m_hEtaPOut[i]->Fill(event.getPOut(kinematicsType).Eta());
		}
		if (event.getEOut(kinematicsType).E() > 0.) { 

			m_resProbEOut[i]++;

			m_hEPOut[i]->Fill(event.getEOut(kinematicsType).P());
			m_hEPtOut[i]->Fill(event.getEOut(kinematicsType).Pt());
			m_hEThOut[i]->Fill(event.getEOut(kinematicsType).Theta());
			m_hEPhOut[i]->Fill(event.getEOut(kinematicsType).Phi());
			m_hEtaEOut[i]->Fill(event.getEOut(kinematicsType).Eta());
		}
		if (event.getGammaOut(kinematicsType).E() > 0.) { 

			m_resProbGOut[i]++;

			m_hGPOut[i]->Fill(event.getGammaOut(kinematicsType).P());
			m_hGPtOut[i]->Fill(event.getGammaOut(kinematicsType).Pt());
			m_hGThOut[i]->Fill(event.getGammaOut(kinematicsType).Theta());
			m_hGPhOut[i]->Fill(event.getGammaOut(kinematicsType).Phi());
			m_hEtaGOut[i]->Fill(event.getGammaOut(kinematicsType).Eta());
		}

		if (event.getPOut(kinematicsType).E() < 0. || event.getEOut(kinematicsType).E() < 0. || event.getGammaOut(kinematicsType).E() < 0.) continue; 

			m_resProbExcl[i]++;

			m_hXBvsQ2->Fill(log10(event.getXB()), log10(event.getQ2()));

			m_hXB[i]->Fill(log10(event.getXB(kinematicsType)));
			m_hQ2[i]->Fill(log10(event.getQ2(kinematicsType)));
			m_hT[i]->Fill(fabs(event.getT(kinematicsType)));
			m_hPhi[i]->Fill(event.getPhi(kinematicsType));
			m_hPhiS[i]->Fill(event.getPhiS(kinematicsType));
			m_hY[i]->Fill(event.getY(kinematicsType));

			m_hxBvsQ2[i]->Fill(log10(event.getXB(kinematicsType)), log10(event.getQ2(kinematicsType)));
			m_hXBvsT[i]->Fill(log10(event.getXB(kinematicsType)), fabs(event.getT(kinematicsType)));
			m_hXBvsY[i]->Fill(log10(event.getXB(kinematicsType)), event.getY(kinematicsType));
			m_hQ2vsT[i]->Fill(log10(event.getQ2(kinematicsType)), fabs(event.getT(kinematicsType)));
			m_hYvsT[i]->Fill(event.getY(kinematicsType), fabs(event.getT(kinematicsType)));
			m_hQ2vsY[i]->Fill(log10(event.getQ2(kinematicsType)), event.getY(kinematicsType));
	}
}

void AnalysisGeneral::analyse(){
	//nothing to be done here
}

void AnalysisGeneral::plot(const std::string& path){

	//print reconstruction probabilities
	std::cout << "info: " << __func__ << ":  p': generated: " << m_resProbPOut[0] << "\treconstructed: " << 
		m_resProbPOut[1] << "\tratio: " << m_resProbPOut[1]/double(m_resProbPOut[0]) << std::endl;
	std::cout << "info: " << __func__ << ":  e': generated: " << m_resProbEOut[0] << "\treconstructed: " << 
		m_resProbEOut[1] << "\tratio: " << m_resProbEOut[1]/double(m_resProbEOut[0]) << std::endl;
	std::cout << "info: " << __func__ << ": gam: generated: " << m_resProbGOut[0] << "\treconstructed: " << 
		m_resProbGOut[1] << "\tratio: " << m_resProbGOut[1]/double(m_resProbGOut[0]) << std::endl;
	std::cout << "info: " << __func__ << ": all: generated: " << m_resProbExcl[0] << "\treconstructed: " << 
		m_resProbExcl[1] << "\tratio: " << m_resProbExcl[1]/double(m_resProbExcl[0]) << std::endl;

	//Clone and make ratio of histograms
	m_hRatio[0] = evaluateAcceptance(m_hPPOut);
	m_hRatio[1] = evaluateAcceptance(m_hPPtOut);
	m_hRatio[2] = evaluateAcceptance(m_hPThOut);
	m_hRatio[3] = evaluateAcceptance(m_hPPhOut);
	m_hRatio[4] = evaluateAcceptance(m_hEtaPOut);

	m_hRatio[5] = evaluateAcceptance(m_hEPOut);
	m_hRatio[6] = evaluateAcceptance(m_hEPtOut);
	m_hRatio[7] = evaluateAcceptance(m_hEThOut);
	m_hRatio[8] = evaluateAcceptance(m_hEPhOut);
	m_hRatio[9] = evaluateAcceptance(m_hEtaEOut);
	
	m_hRatio[10] = evaluateAcceptance(m_hGPOut);
	m_hRatio[11] = evaluateAcceptance(m_hGPtOut);
	m_hRatio[12] = evaluateAcceptance(m_hGThOut);
	m_hRatio[13] = evaluateAcceptance(m_hGPhOut);
	m_hRatio[14] = evaluateAcceptance(m_hEtaGOut);
	
	m_hRatio[15] = evaluateAcceptance(m_hXB);
	m_hRatio[16] = evaluateAcceptance(m_hQ2);
	m_hRatio[17] = evaluateAcceptance(m_hT);
	m_hRatio[18] = evaluateAcceptance(m_hY);
	m_hRatio[19] = evaluateAcceptance(m_hPhi);
	m_hRatio[20] = evaluateAcceptance(m_hPhiS);

	//canvases
	std::vector<TCanvas*> cans;
	
	leg[0] = new TLegend(0.4,0.65,0.75,0.8);
	leg[1] = new TLegend(0.52,0.55,0.85,0.7);
	leg[2] = new TLegend(0.5,0.65,0.85,0.8);
	leg[3] = new TLegend(0.15,0.65,0.5,0.8);
	leg[4] = new TLegend(0.5,0.15,0.85,0.3);
	leg[5] = new TLegend(0.15,0.65,0.5,0.8);

	//write
	m_hXBvsQ2->Write();

	for(size_t i = 0; i < 2; i++){

		m_hXB[i]->Write();
		m_hQ2[i]->Write();
		m_hT[i]->Write();
		m_hY[i]->Write();

		m_hEtaPOut[i]->Write();
		m_hEtaEOut[i]->Write();
		m_hEtaGOut[i]->Write();
	}

	m_hRatio[15]->Write();
	m_hRatio[16]->Write();
	m_hRatio[17]->Write();
	m_hRatio[18]->Write();

	m_hRatio[4]->Write();
	m_hRatio[9]->Write();
	m_hRatio[14]->Write();

	//loop over canvases
	for(size_t i = 0; i < 6; i++){

		//add canvas
		cans.push_back(
			new TCanvas(
				(HashManager::getInstance()->getHash()).c_str(), "")
		);
		leg[i]->SetFillColor(10);
		leg[i]->SetLineColor(10);
		leg[i]->SetShadowColor(10);
		leg[i]->AddEntry(m_hXB[0],"Generated");
		leg[i]->AddEntry(m_hXB[1],"Full reconstruction");
			
		for(size_t k = 0; k <21; k++){
			m_hRatio[k]->SetStats(0);
			m_hRatio[k]->SetTitle("Reconstructed / Generated");
			m_hRatio[k]->SetMinimum(0.);
			m_hRatio[k]->SetMaximum(3.);
		}

		if (i == 0){

			cans.back()->cd();
			m_hXBvsQ2->Draw("colz");
			m_hXBvsQ2->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsQ2->GetYaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			cans.back()->SetLogz();

		} if (i == 1){

			cans.back()->Divide(4,3);

			std::string options;
			int color;

			for(size_t j = 0; j < 2; j++){

				if(j == 0){
					options = "";
					color = 1;
				}else{
					options = "same";
					color = 2;
				}

				cans.back()->cd(1);
				m_hXB[j]->SetLineColor(color);
				m_hXB[j]->Draw(options.c_str());
				m_hXB[j]->GetXaxis()->SetTitle("log(X_{B})");
				leg[1]->Draw();

				cans.back()->cd(2);
				m_hRatio[15]->GetXaxis()->SetTitle("log(X_{B})");
				m_hRatio[15]->Draw();

				cans.back()->cd(3);
				cans.back()->cd(3)->SetLogy();
				m_hQ2[j]->SetLineColor(color);
				m_hQ2[j]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
				m_hQ2[j]->Draw(options.c_str());
				
				cans.back()->cd(4);
				m_hRatio[16]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
				m_hRatio[16]->Draw();

				cans.back()->cd(5);
				cans.back()->cd(5)->SetLogy();
				m_hT[j]->SetLineColor(color);
				m_hT[j]->GetXaxis()->SetTitle("|t| [(GeV/c)^{2}]");
				m_hT[j]->Draw(options.c_str());

				cans.back()->cd(6);
				m_hRatio[17]->GetXaxis()->SetTitle("|t| [(GeV/c)^{2}]");
				m_hRatio[17]->Draw();

				cans.back()->cd(7);
				m_hY[j]->SetLineColor(color);
				m_hY[j]->GetXaxis()->SetTitle("y");
				m_hY[j]->Draw(options.c_str());
				
				cans.back()->cd(8);
				m_hRatio[18]->GetXaxis()->SetTitle("y");
				m_hRatio[18]->Draw();

				cans.back()->cd(9);
				m_hPhi[j]->SetLineColor(color);
				m_hPhi[j]->SetMinimum(0.);
				m_hPhi[j]->GetXaxis()->SetTitle("#phi [rad]");
				m_hPhi[j]->Draw(options.c_str());

				cans.back()->cd(10);
				m_hRatio[19]->GetXaxis()->SetTitle("#phi [rad]");
				m_hRatio[19]->Draw();

				cans.back()->cd(11);
				m_hPhiS[j]->SetLineColor(color);
				m_hPhiS[j]->SetMinimum(0.);
				m_hPhiS[j]->GetXaxis()->SetTitle("#phi_{S} [rad]");
				m_hPhiS[j]->Draw(options.c_str());

				cans.back()->cd(12);
				m_hRatio[20]->GetXaxis()->SetTitle("#phi_{S} [rad]");
				m_hRatio[20]->Draw();

			}
		}if (i == 2){
			
			cans.back()->Divide(4,4);
			
			cans.back()->cd(1);
			m_hxBvsQ2[0]->Draw("colz");
			m_hxBvsQ2[0]->GetXaxis()->SetTitle("log(X_{B})");
			m_hxBvsQ2[0]->GetYaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			cans.back()->SetLogz();
			
			cans.back()->cd(2);
			m_hxBvsQ2[1]->Draw("colz");
			m_hxBvsQ2[1]->GetXaxis()->SetTitle("log(X_{B})");
			m_hxBvsQ2[1]->GetYaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			cans.back()->SetLogz();

			cans.back()->cd(3);
			m_hXBvsT[0]->Draw("colz");
			m_hXBvsT[0]->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsT[0]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();
			
			cans.back()->cd(4);
			m_hXBvsT[1]->Draw("colz");
			m_hXBvsT[1]->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsT[1]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();
			
			cans.back()->cd(5);
			m_hQ2vsT[0]->Draw("colz");
			m_hQ2vsT[0]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			m_hQ2vsT[0]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();
			
			cans.back()->cd(6);
			m_hQ2vsT[1]->Draw("colz");
			m_hQ2vsT[1]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			m_hQ2vsT[1]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();

			cans.back()->cd(7);
			m_hYvsT[0]->Draw("colz");
			m_hYvsT[0]->GetXaxis()->SetTitle("y");
			m_hYvsT[0]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();
			
			cans.back()->cd(8);
			m_hYvsT[1]->Draw("colz");
			m_hYvsT[1]->GetXaxis()->SetTitle("y");
			m_hYvsT[1]->GetYaxis()->SetTitle("|t| [(GeV/c)^{2}]");
			cans.back()->SetLogz();
			
			cans.back()->cd(9);
			m_hQ2vsY[0]->Draw("colz");
			m_hQ2vsY[0]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			m_hQ2vsY[0]->GetYaxis()->SetTitle("y");
			cans.back()->SetLogz();
			
			cans.back()->cd(10);
			m_hQ2vsY[1]->Draw("colz");
			m_hQ2vsY[1]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
			m_hQ2vsY[1]->GetYaxis()->SetTitle("y");
			cans.back()->SetLogz();
			
			cans.back()->cd(11);
			m_hXBvsY[0]->Draw("colz");
			m_hXBvsY[0]->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsY[0]->GetYaxis()->SetTitle("y");
			cans.back()->SetLogz();
			
			cans.back()->cd(12);
			m_hXBvsY[1]->Draw("colz");
			m_hXBvsY[1]->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsY[1]->GetYaxis()->SetTitle("y");
			cans.back()->SetLogz();

		} if (i == 3){

			cans.back()->Divide(4,3);

			std::string options;
			int color;
	
			for(size_t j = 0; j < 2; j++){

				if(j == 0){
					options = "";
					color = 1;
				}else{
					options = "same";
					color = 2;
				}
				
				cans.back()->cd(1);
				cans.back()->cd(1)->SetLogy();
				m_hPPOut[j]->SetLineColor(color);
				m_hPPOut[j]->GetXaxis()->SetTitle("p p' [GeV/c]");
				m_hPPOut[j]->Draw(options.c_str());
				leg[3]->Draw();
			
				cans.back()->cd(2);
				m_hRatio[0]->GetXaxis()->SetTitle("p p' [GeV/c]");
				m_hRatio[0]->Draw();

				cans.back()->cd(3);
				cans.back()->cd(3)->SetLogy();
				m_hPPtOut[j]->SetLineColor(color);
				m_hPPtOut[j]->GetXaxis()->SetTitle("pt p' [GeV/c]");
				m_hPPtOut[j]->Draw(options.c_str());
				
				cans.back()->cd(4);
				m_hRatio[1]->GetXaxis()->SetTitle("pt p' [GeV/c]");
				m_hRatio[1]->Draw();

				cans.back()->cd(5);
				cans.back()->cd(5)->SetLogy();
				m_hEPOut[j]->SetLineColor(color);
				m_hEPOut[j]->GetXaxis()->SetTitle("p e' [GeV/c]");
				m_hEPOut[j]->Draw(options.c_str());

				cans.back()->cd(6);
				m_hRatio[5]->GetXaxis()->SetTitle("p e' [GeV/c]");
				m_hRatio[5]->Draw();

				cans.back()->cd(7);
				cans.back()->cd(7)->SetLogy();
				m_hEPtOut[j]->SetLineColor(color);
				m_hEPtOut[j]->GetXaxis()->SetTitle("pt e' [GeV/c]");
				m_hEPtOut[j]->Draw(options.c_str());

				cans.back()->cd(8);
				m_hRatio[6]->GetXaxis()->SetTitle("pt e' [GeV/c]");
				m_hRatio[6]->Draw();

				cans.back()->cd(9);
				cans.back()->cd(9)->SetLogy();
				m_hGPOut[j]->SetLineColor(color);
				m_hGPOut[j]->GetXaxis()->SetTitle("p #gamma' [GeV/c]");
				m_hGPOut[j]->Draw(options.c_str());

				cans.back()->cd(10);
				m_hRatio[10]->GetXaxis()->SetTitle("p #gamma' [GeV/c]");
				m_hRatio[10]->Draw();

				cans.back()->cd(11);
				cans.back()->cd(11)->SetLogy();
				m_hGPtOut[j]->SetLineColor(color);
				m_hGPtOut[j]->GetXaxis()->SetTitle("pt #gamma' [GeV/c]");
				m_hGPtOut[j]->Draw(options.c_str());
				
				cans.back()->cd(12);
				m_hRatio[11]->GetXaxis()->SetTitle("pt #gamma' [GeV/c]");
				m_hRatio[11]->Draw();

			}
		} if (i == 4){

			cans.back()->Divide(4,3);

			std::string options;
			int color;
	
			for(size_t j = 0; j < 2; j++){

				if(j == 0){
					options = "";
					color = 1;
				}else{
					options = "same";
					color = 2;
				}
				

				cans.back()->cd(1);
				cans.back()->cd(1)->SetLogy();
				m_hPThOut[j]->SetLineColor(color);
				m_hPThOut[j]->GetXaxis()->SetTitle("#theta p' [mrad]");
				m_hPThOut[j]->Draw(options.c_str());
				leg[4]->Draw();
				
				cans.back()->cd(2);
				m_hRatio[2]->GetXaxis()->SetTitle("#theta p' [mrad]");
				m_hRatio[2]->Draw();

				cans.back()->cd(3);
				//cans.back()->cd(3)->SetLogy();
				m_hPPhOut[j]->SetLineColor(color);
				//m_hPPhOut[j]->SetMinimum(100.);
				m_hPPhOut[j]->GetXaxis()->SetTitle("#phi p' [rad]");
				m_hPPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(4);
				m_hRatio[3]->GetXaxis()->SetTitle("#phi p' [rad]");
				m_hRatio[3]->Draw();

				cans.back()->cd(5);
				cans.back()->cd(5)->SetLogy();
				m_hEThOut[j]->SetLineColor(color);
				m_hEThOut[j]->GetXaxis()->SetTitle("#theta e' [#murad]");
				m_hEThOut[j]->Draw(options.c_str());

				cans.back()->cd(6);
				m_hRatio[7]->GetXaxis()->SetTitle("#theta e' [#murad]");
				m_hRatio[7]->Draw();

				cans.back()->cd(7);
				//cans.back()->cd(7)->SetLogy();
				m_hEPhOut[j]->SetLineColor(color);
				//m_hEPhOut[j]->SetMinimum(100.);
				m_hEPhOut[j]->GetXaxis()->SetTitle("#phi e' [rad]");
				m_hEPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(8);
				m_hRatio[8]->GetXaxis()->SetTitle("#phi e' [rad]");
				m_hRatio[8]->Draw();

				cans.back()->cd(9);
				cans.back()->cd(9)->SetLogy();
				m_hGThOut[j]->SetLineColor(color);
				m_hGThOut[j]->GetXaxis()->SetTitle("#theta #gamma' [#murad]");
				m_hGThOut[j]->Draw(options.c_str());

				cans.back()->cd(10);
				m_hRatio[12]->GetXaxis()->SetTitle("#theta #gamma' [#murad]");
				m_hRatio[12]->Draw();

				cans.back()->cd(11);
				//cans.back()->cd(11)->SetLogy();
				m_hGPhOut[j]->SetLineColor(color);
				//m_hGPhOut[j]->SetMinimum(100.);
				m_hGPhOut[j]->GetXaxis()->SetTitle("#phi #gamma' [rad]");
				m_hGPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(12);
				m_hRatio[13]->GetXaxis()->SetTitle("#phi #gamma' [rad]");
				m_hRatio[13]->Draw();

			}
		} if (i == 5){

			cans.back()->Divide(4,3);

			std::string options;
			int color;

			for(size_t j = 0; j < 2; j++){

				if(j == 0){
					options = "";
					color = 1;
				}else{
					options = "same";
					color = 2;
				}
				
				cans.back()->cd(1);
				m_hEtaPOut[j]->SetLineColor(color);
				m_hEtaPOut[j]->GetXaxis()->SetTitle("#eta p'");
				m_hEtaPOut[j]->Draw(options.c_str());
				leg[5]->Draw();
				
				cans.back()->cd(2);
				m_hRatio[4]->GetXaxis()->SetTitle("#eta p'");
				m_hRatio[4]->Draw();
				
				cans.back()->cd(5);
				m_hEtaEOut[j]->SetLineColor(color);
				m_hEtaEOut[j]->GetXaxis()->SetTitle("#eta e'");
				m_hEtaEOut[j]->Draw(options.c_str());
				
				cans.back()->cd(6);
				m_hRatio[9]->GetXaxis()->SetTitle("#eta e'");
				m_hRatio[9]->Draw();
				
				cans.back()->cd(9);
				m_hEtaGOut[j]->SetLineColor(color);
				m_hEtaGOut[j]->GetXaxis()->SetTitle("#eta #gamma'");
				m_hEtaGOut[j]->Draw(options.c_str());
			
				cans.back()->cd(10);
				m_hRatio[14]->GetXaxis()->SetTitle("#eta #gamma'");
				m_hRatio[14]->Draw();
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

