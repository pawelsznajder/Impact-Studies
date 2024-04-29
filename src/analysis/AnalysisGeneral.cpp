#include "../../include/analysis/AnalysisGeneral.h"

#include <cmath>
#include <sstream>
#include <TCanvas.h>
#include <TLegend.h>

#include "../../include/other/HashManager.h"
#include "../../include/kinematic_cuts/KinematicCuts.h"

AnalysisGeneral::AnalysisGeneral(double targetLuminosity) : Analysis("AnalysisGeneral", targetLuminosity), 
	m_lumiM(0.), m_lumiP(0.){

	//number of bins 
	int nbin = 100;

	//cuts
	m_resCut[0] = 0;
	m_resCut[1] = 0;

	//reconstruction probabilities
	for(size_t i = 0; i < 3; i++){

		m_resProbPOut[i] = 0;
		m_resProbEOut[i] = 0;
		m_resProbGOut[i] = 0;

		m_resProbExcl[i] = 0;
	}

	//set 1D and 2D histograms
	for(size_t i = 0; i < 3; i++){
			
		std::string title;
	
		if(i == 0) title = " (observed)";
		if(i == 1) title = " (true)";
		if(i == 2) title = " (born)";
	
		m_hxBvsQ2[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, -5., 0., nbin, 0., 2.);
		m_hxBvsTOverQ2[i] = new TH2D((HashManager::getInstance()->getHash()).c_str(), title.c_str(), nbin, -5., 0., nbin, -3., 1.);

		m_hXB[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "log(x_{B})",nbin, -5., 0.);
		m_hQ2[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "log(Q^{2})",nbin, 0., 2.);
		m_hT[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "|t|",nbin, 0.05, 1.2);
		m_hPhi[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#varphi [rad]",nbin, 0., 2.*TMath::Pi());
		m_hPhiS[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#varphi_{S} [rad]",nbin, 0., 2.*TMath::Pi());
		m_hY[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "y",nbin, 0.05, 0.6);
		m_hEtaEOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta e'",nbin, -4., 0.);
		m_hEtaPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta p'",nbin, 3., 9.);
		m_hEtaGOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "#eta #gamma'",nbin, -10., 5.);
		m_hPPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Momentum p' [GeV/c]",nbin, 50., 110.);
		m_hPPtOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Transverse momentum p' [GeV/c]",nbin, 0., 1.4);
		m_hPThOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Polar angle #theta p' [mrad]",nbin, 0., 25.);
		m_hPPhOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Azimuthal angle #phi p' [rad]",nbin, -TMath::Pi()-0.2, TMath::Pi()+0.2);
		m_hEPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Momentum e' [GeV/c]",nbin, 0., 20.);
		m_hEPtOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Transverse momentum e' [GeV/c]",nbin, 0.,10.);
		m_hEThOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Polar angle #theta e' [#murad]",nbin, 0., 5.);
		m_hEPhOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Azimuthal angle #phi e' [rad]",nbin, -TMath::Pi()-0.2, TMath::Pi()+0.2);
		m_hGPOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Momentum #gamma' [GeV/c]",nbin, 0., 40.);
		m_hGPtOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Transverse momentum #gamma' [GeV/c]",nbin, 0., 10.);
		m_hGThOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Polar angle #theta #gamma' [#murad]",nbin, 0., 5.);
		m_hGPhOut[i] = new TH1D((HashManager::getInstance()->getHash()).c_str(), "Azimuthal angle #phi #gamma' [rad]",nbin, -TMath::Pi()-0.2, TMath::Pi()+0.2);
	}
}

AnalysisGeneral::~AnalysisGeneral(){
}

TH1* AnalysisGeneral::evaluateAcceptance(TH1** h, bool isRC) const{

	TH1* acc;

	if(! isRC){
		acc = (TH1*)h[0]->Clone();
		acc->SetName(("acc_"+std::string(h[0]->GetName())).c_str());
		acc->Divide(h[1]);
	}else{
		acc = (TH1*)h[1]->Clone();
		acc->SetName(("rc_"+std::string(h[1]->GetName())).c_str());
		acc->Divide(h[2]);
	}

	return acc;
}

TH2* AnalysisGeneral::evaluateAcceptance(TH2** h, bool isRC) const{

	TH2* acc;

	if(! isRC){
		acc = (TH2*)h[0]->Clone();
		acc->SetName(("acc_"+std::string(h[0]->GetName())).c_str());
		acc->Divide(h[1]);
	}else{
		acc = (TH2*)h[1]->Clone();
		acc->SetName(("rc_"+std::string(h[1]->GetName())).c_str());
		acc->Divide(h[2]);
	}

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

	//only to reach target luminosity
	switch(event.getBeamPolarisation()){

   	 	case -1:{

   	 		if(m_lumiM >= m_targetLuminosity/2.) return;
   	 		m_lumiM += weight;

   	 		break;
   	 	}

       	case 1:{

			if(m_lumiP >= m_targetLuminosity/2.) return;
   	 		m_lumiP += weight;

       	 	break;
       	}

       default:{
       		std::cout << __func__ << " error: wrong beam polarisation state, " << event.getBeamPolarisation() << std::endl;
			exit(0);
       }
	}

	//number of events after restricting the kinematic donain
	m_resCut[1]++;
	if(KinematicCuts::checkKinematicCuts(event, KinematicsType::True)) m_resCut[0]++;

	//fill 1D histograms
	//only for interesting (physics-wise motivated) kinematics (i.e. restricted kinematics)
	if(KinematicCuts::checkKinematicCuts(event, KinematicsType::True)){

		for(size_t i = 0; i < 3; i++){

			//skip Born
			if(i == 2) continue;

			//type (kinematicsTypePlotted - only geometric acceptance here)
			KinematicsType::Type kinematicsType, kinematicsTypePlotted;

			if(i == 0) {kinematicsType = KinematicsType::Observed; kinematicsTypePlotted = KinematicsType::True;}
			if(i == 1) {kinematicsType = KinematicsType::True; kinematicsTypePlotted = kinematicsType;}
			if(i == 2) {kinematicsType = KinematicsType::Born; kinematicsTypePlotted = kinematicsType;}

			if (event.getPOut(kinematicsType).E() > 0.) {

				m_resProbPOut[i]++;

				m_hPPOut[i]->Fill(event.getPOut(kinematicsTypePlotted).P());
				m_hPPtOut[i]->Fill(event.getPOut(kinematicsTypePlotted).Pt());
				m_hPThOut[i]->Fill(event.getPOut(kinematicsTypePlotted).Theta()*1000.);
				m_hPPhOut[i]->Fill(event.getPOut(kinematicsTypePlotted).Phi());
				m_hEtaPOut[i]->Fill(event.getPOut(kinematicsTypePlotted).Eta());
			}

			if (event.getEOut(kinematicsType).E() > 0.) { 

				m_resProbEOut[i]++;

				m_hEPOut[i]->Fill(event.getEOut(kinematicsTypePlotted).P());
				m_hEPtOut[i]->Fill(event.getEOut(kinematicsTypePlotted).Pt());
				m_hEThOut[i]->Fill(event.getEOut(kinematicsTypePlotted).Theta());
				m_hEPhOut[i]->Fill(event.getEOut(kinematicsTypePlotted).Phi());
				m_hEtaEOut[i]->Fill(event.getEOut(kinematicsTypePlotted).Eta());
			}

			if (event.getGammaOut(kinematicsType).E() > 0.) { 

				m_resProbGOut[i]++;

				m_hGPOut[i]->Fill(event.getGammaOut(kinematicsTypePlotted).P());
				m_hGPtOut[i]->Fill(event.getGammaOut(kinematicsTypePlotted).Pt());
				m_hGThOut[i]->Fill(event.getGammaOut(kinematicsTypePlotted).Theta());
				m_hGPhOut[i]->Fill(event.getGammaOut(kinematicsTypePlotted).Phi());
				m_hEtaGOut[i]->Fill(event.getGammaOut(kinematicsTypePlotted).Eta());
			}

			if(event.getPOut(kinematicsType).E() > 0. && event.getEOut(kinematicsType).E() > 0. && event.getGammaOut(kinematicsType).E() > 0.){
				m_resProbExcl[i]++;
			}
		}
	}

	for(size_t i = 0; i < 3; i++){

		//skip Born
		if(i == 2) continue;

		//type
		KinematicsType::Type kinematicsType;

		if(i == 0) kinematicsType = KinematicsType::Observed;
		if(i == 1) kinematicsType = KinematicsType::True;
		if(i == 2) kinematicsType = KinematicsType::Born;

		//the rest acts as backround
		if (event.getPOut(kinematicsType).E() > 0. && event.getEOut(kinematicsType).E() > 0. && event.getGammaOut(kinematicsType).E() > 0.){

			if(KinematicCuts::checkKinematicCuts(event, kinematicsType), 3)
				m_hXB[i]->Fill(log10(event.getXB(kinematicsType)));
			if(KinematicCuts::checkKinematicCuts(event, kinematicsType), 2)
				m_hQ2[i]->Fill(log10(event.getQ2(kinematicsType)));

			m_hPhi[i]->Fill(event.getPhi(kinematicsType));
			m_hPhiS[i]->Fill(event.getPhiS(kinematicsType));

			if(KinematicCuts::checkKinematicCuts(event, kinematicsType, 0))
				m_hT[i]->Fill(fabs(event.getT(kinematicsType)));

			if(KinematicCuts::checkKinematicCuts(event, kinematicsType, 1))
				m_hY[i]->Fill(event.getY(kinematicsType));

			if(KinematicCuts::checkKinematicCuts(event, kinematicsType))
				m_hxBvsQ2[i]->Fill(log10(event.getXB(kinematicsType)), log10(event.getQ2(kinematicsType)));
			if(KinematicCuts::checkKinematicCuts(event, kinematicsType))
				m_hxBvsTOverQ2[i]->Fill(log10(event.getXB(kinematicsType)), log10(fabs(event.getT(kinematicsType))/event.getQ2(kinematicsType)));
		}
	}
}

void AnalysisGeneral::analyse(){
	//nothing to be done here
}

void AnalysisGeneral::plot(const std::string& path){

	//print reconstruction probabilities
	std::cout << "info: " << __func__ << ":  resticted kinematic domain: all " << m_resCut[1] << "\tin: " << 
		m_resCut[0] << "\tratio: " << m_resCut[0]/double(m_resCut[1]) << std::endl;
	std::cout << "info: " << __func__ << ":  p': generated: " << m_resProbPOut[1] << "\treconstructed: " << 
		m_resProbPOut[0] << "\tratio: " << m_resProbPOut[0]/double(m_resProbPOut[1]) << std::endl;
	std::cout << "info: " << __func__ << ":  e': generated: " << m_resProbEOut[1] << "\treconstructed: " << 
		m_resProbEOut[0] << "\tratio: " << m_resProbEOut[0]/double(m_resProbEOut[1]) << std::endl;
	std::cout << "info: " << __func__ << ": gam: generated: " << m_resProbGOut[1] << "\treconstructed: " << 
		m_resProbGOut[0] << "\tratio: " << m_resProbGOut[0]/double(m_resProbGOut[1]) << std::endl;
	std::cout << "info: " << __func__ << ": exc: generated: " << m_resProbExcl[1] << "\treconstructed: " << 
		m_resProbExcl[0] << "\tratio: " << m_resProbExcl[0]/double(m_resProbExcl[1]) << std::endl;

	//Clone and make ratio of histograms
	for(int i = 0; i < 2; i++){

		m_hxBvsQ2[3+i] = evaluateAcceptance(m_hxBvsQ2, i);
		m_hxBvsTOverQ2[3+i] = evaluateAcceptance(m_hxBvsTOverQ2, i);

		m_hPPOut[3+i] = evaluateAcceptance(m_hPPOut, i);
		m_hPPtOut[3+i] = evaluateAcceptance(m_hPPtOut, i);
		m_hPThOut[3+i] = evaluateAcceptance(m_hPThOut, i);
		m_hPPhOut[3+i] = evaluateAcceptance(m_hPPhOut, i);
		m_hEtaPOut[3+i] = evaluateAcceptance(m_hEtaPOut, i);

		m_hEPOut[3+i] = evaluateAcceptance(m_hEPOut, i);
		m_hEPtOut[3+i] = evaluateAcceptance(m_hEPtOut, i);
		m_hEThOut[3+i] = evaluateAcceptance(m_hEThOut, i);
		m_hEPhOut[3+i] = evaluateAcceptance(m_hEPhOut, i);
		m_hEtaEOut[3+i] = evaluateAcceptance(m_hEtaEOut, i);
		
		m_hGPOut[3+i] = evaluateAcceptance(m_hGPOut, i);
		m_hGPtOut[3+i] = evaluateAcceptance(m_hGPtOut, i);
		m_hGThOut[3+i] = evaluateAcceptance(m_hGThOut, i);
		m_hGPhOut[3+i] = evaluateAcceptance(m_hGPhOut, i);
		m_hEtaGOut[3+i] = evaluateAcceptance(m_hEtaGOut, i);
		
		m_hXB[3+i] = evaluateAcceptance(m_hXB, i);
		m_hQ2[3+i] = evaluateAcceptance(m_hQ2, i);
		m_hT[3+i] = evaluateAcceptance(m_hT, i);
		m_hY[3+i] = evaluateAcceptance(m_hY, i);
		m_hPhi[3+i] = evaluateAcceptance(m_hPhi, i);
		m_hPhiS[3+i] = evaluateAcceptance(m_hPhiS, i);
	}

	//canvases
	std::vector<TCanvas*> cans;
	
	leg[0] = new TLegend(0.4,0.65,0.75,0.8);
	leg[1] = new TLegend(0.52,0.55,0.85,0.7);
	leg[2] = new TLegend(0.5,0.65,0.85,0.8);
	leg[3] = new TLegend(0.15,0.65,0.5,0.8);
	leg[4] = new TLegend(0.5,0.15,0.85,0.3);
	leg[5] = new TLegend(0.15,0.65,0.5,0.8);

	//write
	for(size_t i = 0; i < 5; i++){

		m_hxBvsQ2[i]->Write();
		m_hxBvsTOverQ2[i]->Write();

		m_hXB[i]->Write();
		m_hQ2[i]->Write();
		m_hT[i]->Write();
		m_hY[i]->Write();

		m_hEtaPOut[i]->Write();
		m_hEtaEOut[i]->Write();
		m_hEtaGOut[i]->Write();
	}

	//loop over canvases
	for(size_t i = 0; i < 5; i++){

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

			cans.back()->Divide(2,2);

			cans.back()->cd(1);
			m_hxBvsQ2[0]->Draw("colz");
			m_hxBvsQ2[0]->GetXaxis()->SetTitle("log(x_{B})");
			m_hxBvsQ2[0]->GetYaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");

			cans.back()->cd(2);
			m_hxBvsQ2[1]->Draw("colz");
			m_hxBvsQ2[1]->GetXaxis()->SetTitle("log(x_{B})");
			m_hxBvsQ2[1]->GetYaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");

			cans.back()->cd(3);
			m_hxBvsTOverQ2[0]->Draw("colz");
			m_hxBvsTOverQ2[0]->GetXaxis()->SetTitle("log(x_{B})");
			m_hxBvsTOverQ2[0]->GetYaxis()->SetTitle("log(|t|/Q^{2})");

			cans.back()->cd(4);
			m_hxBvsTOverQ2[1]->Draw("colz");
			m_hxBvsTOverQ2[1]->GetXaxis()->SetTitle("log(x_{B})");
			m_hxBvsTOverQ2[1]->GetYaxis()->SetTitle("log(|t|/Q^{2})");
		} 

		if (i == 1){

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
				m_hXB[j]->GetXaxis()->SetTitle("log(x_{B})");
				leg[1]->Draw();

				cans.back()->cd(2);
				m_hXB[3]->GetXaxis()->SetTitle("log(x_{B})");
				m_hXB[3]->Draw();

				cans.back()->cd(3);
				cans.back()->cd(3)->SetLogy();
				m_hQ2[j]->SetLineColor(color);
				m_hQ2[j]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
				m_hQ2[j]->Draw(options.c_str());
				
				cans.back()->cd(4);
				m_hQ2[3]->GetXaxis()->SetTitle("log(Q^{2} [(GeV/c)^{2}])");
				m_hQ2[3]->Draw();

				cans.back()->cd(5);
				cans.back()->cd(5)->SetLogy();
				m_hT[j]->SetLineColor(color);
				m_hT[j]->GetXaxis()->SetTitle("|t| [(GeV/c)^{2}]");
				m_hT[j]->Draw(options.c_str());

				cans.back()->cd(6);
				m_hT[3]->GetXaxis()->SetTitle("|t| [(GeV/c)^{2}]");
				m_hT[3]->Draw();

				cans.back()->cd(7);
				m_hY[j]->SetLineColor(color);
				m_hY[j]->GetXaxis()->SetTitle("y");
				m_hY[j]->Draw(options.c_str());
				
				cans.back()->cd(8);
				m_hY[3]->GetXaxis()->SetTitle("y");
				m_hY[3]->Draw();

				cans.back()->cd(9);
				m_hPhi[j]->SetLineColor(color);
				m_hPhi[j]->SetMinimum(0.);
				m_hPhi[j]->GetXaxis()->SetTitle("#phi [rad]");
				m_hPhi[j]->Draw(options.c_str());

				cans.back()->cd(10);
				m_hPhi[3]->GetXaxis()->SetTitle("#phi [rad]");
				m_hPhi[3]->Draw();

				cans.back()->cd(11);
				m_hPhiS[j]->SetLineColor(color);
				m_hPhiS[j]->SetMinimum(0.);
				m_hPhiS[j]->GetXaxis()->SetTitle("#phi_{S} [rad]");
				m_hPhiS[j]->Draw(options.c_str());

				cans.back()->cd(12);
				m_hPhiS[3]->GetXaxis()->SetTitle("#phi_{S} [rad]");
				m_hPhiS[3]->Draw();

			}
		}

		if (i == 2){

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
				m_hPPOut[3]->GetXaxis()->SetTitle("p p' [GeV/c]");
				m_hPPOut[3]->Draw();

				cans.back()->cd(3);
				cans.back()->cd(3)->SetLogy();
				m_hPPtOut[j]->SetLineColor(color);
				m_hPPtOut[j]->GetXaxis()->SetTitle("pt p' [GeV/c]");
				m_hPPtOut[j]->Draw(options.c_str());
				
				cans.back()->cd(4);
				m_hPPtOut[3]->GetXaxis()->SetTitle("pt p' [GeV/c]");
				m_hPPtOut[3]->Draw();

				cans.back()->cd(5);
				cans.back()->cd(5)->SetLogy();
				m_hEPOut[j]->SetLineColor(color);
				m_hEPOut[j]->GetXaxis()->SetTitle("p e' [GeV/c]");
				m_hEPOut[j]->Draw(options.c_str());

				cans.back()->cd(6);
				m_hEPOut[3]->GetXaxis()->SetTitle("p e' [GeV/c]");
				m_hEPOut[3]->Draw();

				cans.back()->cd(7);
				cans.back()->cd(7)->SetLogy();
				m_hEPtOut[j]->SetLineColor(color);
				m_hEPtOut[j]->GetXaxis()->SetTitle("pt e' [GeV/c]");
				m_hEPtOut[j]->Draw(options.c_str());

				cans.back()->cd(8);
				m_hEPtOut[3]->GetXaxis()->SetTitle("pt e' [GeV/c]");
				m_hEPtOut[3]->Draw();

				cans.back()->cd(9);
				cans.back()->cd(9)->SetLogy();
				m_hGPOut[j]->SetLineColor(color);
				m_hGPOut[j]->GetXaxis()->SetTitle("p #gamma' [GeV/c]");
				m_hGPOut[j]->Draw(options.c_str());

				cans.back()->cd(10);
				m_hGPOut[3]->GetXaxis()->SetTitle("p #gamma' [GeV/c]");
				m_hGPOut[3]->Draw();

				cans.back()->cd(11);
				cans.back()->cd(11)->SetLogy();
				m_hGPtOut[j]->SetLineColor(color);
				m_hGPtOut[j]->GetXaxis()->SetTitle("pt #gamma' [GeV/c]");
				m_hGPtOut[j]->Draw(options.c_str());
				
				cans.back()->cd(12);
				m_hGPtOut[3]->GetXaxis()->SetTitle("pt #gamma' [GeV/c]");
				m_hGPtOut[3]->Draw();

			}
		} 

		if (i == 3){

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
				m_hPThOut[3]->GetXaxis()->SetTitle("#theta p' [mrad]");
				m_hPThOut[3]->Draw();

				cans.back()->cd(3);
				//cans.back()->cd(3)->SetLogy();
				m_hPPhOut[j]->SetLineColor(color);
				//m_hPPhOut[j]->SetMinimum(100.);
				m_hPPhOut[j]->GetXaxis()->SetTitle("#phi p' [rad]");
				m_hPPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(4);
				m_hPPhOut[3]->GetXaxis()->SetTitle("#phi p' [rad]");
				m_hPPhOut[3]->Draw();

				cans.back()->cd(5);
				cans.back()->cd(5)->SetLogy();
				m_hEThOut[j]->SetLineColor(color);
				m_hEThOut[j]->GetXaxis()->SetTitle("#theta e' [#murad]");
				m_hEThOut[j]->Draw(options.c_str());

				cans.back()->cd(6);
				m_hEThOut[3]->GetXaxis()->SetTitle("#theta e' [#murad]");
				m_hEThOut[3]->Draw();

				cans.back()->cd(7);
				//cans.back()->cd(7)->SetLogy();
				m_hEPhOut[j]->SetLineColor(color);
				//m_hEPhOut[j]->SetMinimum(100.);
				m_hEPhOut[j]->GetXaxis()->SetTitle("#phi e' [rad]");
				m_hEPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(8);
				m_hEPhOut[3]->GetXaxis()->SetTitle("#phi e' [rad]");
				m_hEPhOut[3]->Draw();

				cans.back()->cd(9);
				cans.back()->cd(9)->SetLogy();
				m_hGThOut[j]->SetLineColor(color);
				m_hGThOut[j]->GetXaxis()->SetTitle("#theta #gamma' [#murad]");
				m_hGThOut[j]->Draw(options.c_str());

				cans.back()->cd(10);
				m_hGThOut[3]->GetXaxis()->SetTitle("#theta #gamma' [#murad]");
				m_hGThOut[3]->Draw();

				cans.back()->cd(11);
				//cans.back()->cd(11)->SetLogy();
				m_hGPhOut[j]->SetLineColor(color);
				//m_hGPhOut[j]->SetMinimum(100.);
				m_hGPhOut[j]->GetXaxis()->SetTitle("#phi #gamma' [rad]");
				m_hGPhOut[j]->Draw(options.c_str());
				
				cans.back()->cd(12);
				m_hGPhOut[3]->GetXaxis()->SetTitle("#phi #gamma' [rad]");
				m_hGPhOut[3]->Draw();

			}
		} 

		if (i == 4){

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
				m_hEtaPOut[3]->GetXaxis()->SetTitle("#eta p'");
				m_hEtaPOut[3]->Draw();
				
				cans.back()->cd(5);
				m_hEtaEOut[j]->SetLineColor(color);
				m_hEtaEOut[j]->GetXaxis()->SetTitle("#eta e'");
				m_hEtaEOut[j]->Draw(options.c_str());
				
				cans.back()->cd(6);
				m_hEtaEOut[3]->GetXaxis()->SetTitle("#eta e'");
				m_hEtaEOut[3]->Draw();
				
				cans.back()->cd(9);
				m_hEtaGOut[j]->SetLineColor(color);
				m_hEtaGOut[j]->GetXaxis()->SetTitle("#eta #gamma'");
				m_hEtaGOut[j]->Draw(options.c_str());
			
				cans.back()->cd(10);
				m_hEtaGOut[3]->GetXaxis()->SetTitle("#eta #gamma'");
				m_hEtaGOut[3]->Draw();
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

