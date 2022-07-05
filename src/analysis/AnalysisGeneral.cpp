#include "../../include/analysis/AnalysisGeneral.h"

#include <cmath>
#include <sstream>
#include <TCanvas.h>
#include <TLegend.h>

#include "../../include/other/HashManager.h"

AnalysisGeneral::AnalysisGeneral() : Analysis("AnalysisGeneral"){

	//number of bins 
	int nbin = 100;

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

void AnalysisGeneral::fill(DVCSEvent& event, double weight){

	//fill 2D histograms
	m_hXBvsQ2->Fill(log10(event.getXB()), log10(event.getQ2()), weight);

	//fill 1D histograms
	for(size_t i = 0; i < 2; i++){

		//type
		KinematicsType::Type kinematicsType = (i == 0)?(KinematicsType::True):(KinematicsType::Observed);

		//check if reconstricted
	//	if(kinematicsType == KinematicsType::Observed && (! event.isReconstructed())) continue;
		if (event.getY(kinematicsType) > 0.05) {

		if (event.getPOut(kinematicsType).E() > 0.) {//continue;
		m_hPPOut[i]->Fill(event.getPOut(kinematicsType).P(), weight);
		m_hPPtOut[i]->Fill(event.getPOut(kinematicsType).Pt(), weight);
		m_hPThOut[i]->Fill(event.getPOut(kinematicsType).Theta()*1000., weight);
		m_hPPhOut[i]->Fill(event.getPOut(kinematicsType).Phi(), weight);
		m_hEtaPOut[i]->Fill(event.getEtaPOut(kinematicsType), weight);
		}
		if (event.getEOut(kinematicsType).E() > 0.) { //continue;
		m_hEPOut[i]->Fill(event.getEOut(kinematicsType).P(), weight);
		m_hEPtOut[i]->Fill(event.getEOut(kinematicsType).Pt(), weight);
		m_hEThOut[i]->Fill(event.getEOut(kinematicsType).Theta(), weight);
		m_hEPhOut[i]->Fill(event.getEOut(kinematicsType).Phi(), weight);
		m_hEtaEOut[i]->Fill(event.getEtaEOut(kinematicsType), weight);
		}
		if (event.getGammaOut(kinematicsType).E() > 0.) { //continue;
		m_hGPOut[i]->Fill(event.getGammaOut(kinematicsType).P(), weight);
		m_hGPtOut[i]->Fill(event.getGammaOut(kinematicsType).Pt(), weight);
		m_hGThOut[i]->Fill(event.getGammaOut(kinematicsType).Theta(), weight);
		m_hGPhOut[i]->Fill(event.getGammaOut(kinematicsType).Phi(), weight);
		m_hEtaGOut[i]->Fill(event.getEtaGOut(kinematicsType), weight);
		}

		if (event.getPOut(kinematicsType).E() < 0. || event.getEOut(kinematicsType).E() < 0. || event.getGammaOut(kinematicsType).E() < 0.) continue; 

		m_hXB[i]->Fill(log10(event.getXB(kinematicsType)), weight);
		m_hQ2[i]->Fill(log10(event.getQ2(kinematicsType)), weight);
		m_hT[i]->Fill(fabs(event.getT(kinematicsType)), weight);
		m_hPhi[i]->Fill(event.getPhi(kinematicsType), weight);
		m_hPhiS[i]->Fill(event.getPhiS(kinematicsType), weight);
		m_hY[i]->Fill(event.getY(kinematicsType), weight);

		m_hxBvsQ2[i]->Fill(log10(event.getXB(kinematicsType)), log10(event.getQ2(kinematicsType)), weight);
		m_hXBvsT[i]->Fill(log10(event.getXB(kinematicsType)), fabs(event.getT(kinematicsType)), weight);
		m_hXBvsY[i]->Fill(log10(event.getXB(kinematicsType)), event.getY(kinematicsType), weight);
		m_hQ2vsT[i]->Fill(log10(event.getQ2(kinematicsType)), fabs(event.getT(kinematicsType)), weight);
		m_hYvsT[i]->Fill(event.getY(kinematicsType), fabs(event.getT(kinematicsType)), weight);
		m_hQ2vsY[i]->Fill(log10(event.getQ2(kinematicsType)), event.getY(kinematicsType), weight);
		}
	}
}

void AnalysisGeneral::analyse(){
	//nothing to be done here
}

void AnalysisGeneral::plot(const std::string& path){

	//canvases
	std::vector<TCanvas*> cans;
	
	leg[0] = new TLegend(0.4,0.65,0.75,0.8);
	leg[1] = new TLegend(0.25,0.65,0.6,0.8);
	leg[2] = new TLegend(0.25,0.65,0.6,0.8);
	leg[3] = new TLegend(0.25,0.65,0.6,0.8);

	//loop over canvases
	for(size_t i = 0; i < 4; i++){

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

		if (i == 0){

			cans.back()->cd();
			m_hXBvsQ2->Draw("colz");
			cans.back()->SetLogz();

		} if (i == 1){

			cans.back()->Divide(3,3);

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
				leg[0]->Draw();

				cans.back()->cd(2);
				cans.back()->cd(2)->SetLogy();
				m_hQ2[j]->SetLineColor(color);
				m_hQ2[j]->Draw(options.c_str());

				cans.back()->cd(3);
				cans.back()->cd(3)->SetLogy();
				m_hT[j]->SetLineColor(color);
				m_hT[j]->Draw(options.c_str());

				cans.back()->cd(4);
				m_hPhi[j]->SetLineColor(color);
				m_hPhi[j]->Draw(options.c_str());

				cans.back()->cd(5);
				m_hPhiS[j]->SetLineColor(color);
				m_hPhiS[j]->SetMinimum(0.);
				m_hPhiS[j]->Draw(options.c_str());

				cans.back()->cd(6);
				m_hY[j]->SetLineColor(color);
				m_hY[j]->Draw(options.c_str());

			}
		}if (i == 2){
			
			cans.back()->Divide(4,4);
			
			cans.back()->cd(1);
			m_hxBvsQ2[0]->Draw("colz");
			m_hxBvsQ2[0]->GetXaxis()->SetTitle("log(X_{B})");
			m_hxBvsQ2[0]->GetYaxis()->SetTitle("log(Q^{2})");
			cans.back()->SetLogz();
			
			cans.back()->cd(2);
			m_hxBvsQ2[1]->Draw("colz");
			m_hxBvsQ2[1]->GetXaxis()->SetTitle("log(X_{B})");
			m_hxBvsQ2[1]->GetYaxis()->SetTitle("log(Q^{2})");
			cans.back()->SetLogz();

			cans.back()->cd(3);
			m_hXBvsT[0]->Draw("colz");
			m_hXBvsT[0]->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsT[0]->GetYaxis()->SetTitle("|t|");
			cans.back()->SetLogz();
			
			cans.back()->cd(4);
			m_hXBvsT[1]->Draw("colz");
			m_hXBvsT[1]->GetXaxis()->SetTitle("log(X_{B})");
			m_hXBvsT[1]->GetYaxis()->SetTitle("|t|");
			cans.back()->SetLogz();
			
			cans.back()->cd(5);
			m_hQ2vsT[0]->Draw("colz");
			m_hQ2vsT[0]->GetXaxis()->SetTitle("log(Q^{2})");
			m_hQ2vsT[0]->GetYaxis()->SetTitle("|t|");
			cans.back()->SetLogz();
			
			cans.back()->cd(6);
			m_hQ2vsT[1]->Draw("colz");
			m_hQ2vsT[1]->GetXaxis()->SetTitle("log(Q^{2})");
			m_hQ2vsT[1]->GetYaxis()->SetTitle("|t|");
			cans.back()->SetLogz();

			cans.back()->cd(7);
			m_hYvsT[0]->Draw("colz");
			m_hYvsT[0]->GetXaxis()->SetTitle("y");
			m_hYvsT[0]->GetYaxis()->SetTitle("|t|");
			cans.back()->SetLogz();
			
			cans.back()->cd(8);
			m_hYvsT[1]->Draw("colz");
			m_hYvsT[1]->GetXaxis()->SetTitle("y");
			m_hYvsT[1]->GetYaxis()->SetTitle("|t|");
			cans.back()->SetLogz();
			
			cans.back()->cd(9);
			m_hQ2vsY[0]->Draw("colz");
			m_hQ2vsY[0]->GetXaxis()->SetTitle("log(Q^{2})");
			m_hQ2vsY[0]->GetYaxis()->SetTitle("y");
			cans.back()->SetLogz();
			
			cans.back()->cd(10);
			m_hQ2vsY[1]->Draw("colz");
			m_hQ2vsY[1]->GetXaxis()->SetTitle("log(Q^{2})");
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

			cans.back()->Divide(5,3);

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
				m_hPPOut[j]->Draw(options.c_str());
				leg[1]->Draw();

				cans.back()->cd(2);
				cans.back()->cd(2)->SetLogy();
				m_hPPtOut[j]->SetLineColor(color);
				m_hPPtOut[j]->Draw(options.c_str());

				cans.back()->cd(3);
				cans.back()->cd(3)->SetLogy();
				m_hPThOut[j]->SetLineColor(color);
				m_hPThOut[j]->Draw(options.c_str());

				cans.back()->cd(4);
				cans.back()->cd(4)->SetLogy();
				m_hPPhOut[j]->SetLineColor(color);
				m_hPPhOut[j]->SetMinimum(100.);
				m_hPPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(5);
				m_hEtaPOut[j]->SetLineColor(color);
				m_hEtaPOut[j]->Draw(options.c_str());
				
				cans.back()->cd(6);
				cans.back()->cd(6)->SetLogy();
				m_hEPOut[j]->SetLineColor(color);
				m_hEPOut[j]->Draw(options.c_str());

				cans.back()->cd(7);
				cans.back()->cd(7)->SetLogy();
				m_hEPtOut[j]->SetLineColor(color);
				m_hEPtOut[j]->Draw(options.c_str());

				cans.back()->cd(8);
				cans.back()->cd(8)->SetLogy();
				m_hEThOut[j]->SetLineColor(color);
				m_hEThOut[j]->Draw(options.c_str());

				cans.back()->cd(9);
				cans.back()->cd(9)->SetLogy();
				m_hEPhOut[j]->SetLineColor(color);
				m_hEPhOut[j]->SetMinimum(100.);
				m_hEPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(10);
				m_hEtaEOut[j]->SetLineColor(color);
				m_hEtaEOut[j]->Draw(options.c_str());
				
				cans.back()->cd(11);
				cans.back()->cd(11)->SetLogy();
				m_hGPOut[j]->SetLineColor(color);
				m_hGPOut[j]->Draw(options.c_str());

				cans.back()->cd(12);
				cans.back()->cd(12)->SetLogy();
				m_hGPtOut[j]->SetLineColor(color);
				m_hGPtOut[j]->Draw(options.c_str());

				cans.back()->cd(13);
				cans.back()->cd(13)->SetLogy();
				m_hGThOut[j]->SetLineColor(color);
				m_hGThOut[j]->Draw(options.c_str());

				cans.back()->cd(14);
				cans.back()->cd(14)->SetLogy();
				m_hGPhOut[j]->SetLineColor(color);
				m_hGPhOut[j]->SetMinimum(100.);
				m_hGPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(15);
				m_hEtaGOut[j]->SetLineColor(color);
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

