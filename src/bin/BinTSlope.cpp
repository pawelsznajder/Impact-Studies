#include "../../include/bin/BinTSlope.h"

#include <iostream>
#include <sstream>

#include "../../include/other/HashManager.h"

BinTSlope::BinTSlope(
	const std::pair<double, double>& rangeXB, 
	const std::pair<double, double>& rangeQ2, 
	size_t nTBins,
	const std::pair<double, double>& rangeT
) : Bin("BinTSlope", rangeXB, rangeQ2, rangeT, std::make_pair(0., 2 * M_PI)), m_lumiALL(0.), m_lumiBH(0.) {

	//reset
	reset();
 
	//make histograms
 	for(size_t i = 0; i < 6; i++){

 		//labels
		std::stringstream ss;

		if(i == 0) ss << "ALL rec lumi: ";
		if(i == 1) ss << "BH true: ";
		if(i == 2) ss << "ALL true: ";
		if(i == 3) ss << "ALL rec: ";

		if(i == 4) ss << "BH born: ";
		if(i == 5) ss << "ALL born: ";

		ss << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
			m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second;
		
		//make histograms
		m_hDistributions.push_back(new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
						nTBins, rangeT.first, rangeT.second));

		//set sumw2
		m_hDistributions.back()->Sumw2();
	}

	//function for fitting
	m_fFit = new TF1((HashManager::getInstance()->getHash()).c_str(), "[0]*exp(-1*[1]*x)", 0., 2.);

	m_fFit->SetParameter(0, 100.);
	m_fFit->SetParameter(1, 5.);

}

BinTSlope::~BinTSlope(){
}

void BinTSlope::reset(){

	//run for parent class
	Bin::reset();

	//reset histograms
	m_hDistributions.clear();
	m_hTSlope = nullptr;
	m_hTAcceptance = nullptr;
	m_hTRC = nullptr;
}

void BinTSlope::fill(DVCSEvent& event, double weight){

	//if BH+INT+DVCS
	if( event.checkSubProcessType(SubProcessType::BH) && 
		event.checkSubProcessType(SubProcessType::INT) && 
		event.checkSubProcessType(SubProcessType::DVCS)
	){

		if(weight > 0.){

			//kinematics
			Bin::fill(event);

			//fill
			if(event.isReconstructed()) m_hDistributions.at(0)->Fill(-1 * event.getT());
		}

		m_hDistributions.at(2)->Fill(-1 * event.getT(KinematicsType::True));
		if(event.isReconstructed()) m_hDistributions.at(3)->Fill(-1 * event.getT());
		m_hDistributions.at(5)->Fill(-1 * event.getT(KinematicsType::Born));

		m_lumiALL += fabs(weight);
	}
	
	//if BH
	if( event.checkSubProcessType(SubProcessType::BH) && 
		(! event.checkSubProcessType(SubProcessType::INT)) && 
		(! event.checkSubProcessType(SubProcessType::DVCS))
	){

		m_hDistributions.at(1)->Fill(-1 * event.getT(KinematicsType::True));
		m_hDistributions.at(4)->Fill(-1 * event.getT(KinematicsType::Born));
		m_lumiBH += weight;
	}
}

void BinTSlope::analyse(){

	std::cout << "error: " << __func__ << ": do not use this function" << std::endl;
	exit(0);
}

void BinTSlope::analyse(double totalLumiALL, double totalLumiBH){
	
	//run for parent class
	Bin::analyse();

	//skip bins with low entry count
	if(getNEvents(KinematicsType::Observed) < 2000){
		return;
	}

	//labels
	std::stringstream ss;

	//make t-slope histogram
	m_hTSlope = static_cast<TH1*>(m_hDistributions.at(0)->Clone());

	ss.clear();
	ss << "DVCS corrected " << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
			m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second;
	m_hTSlope->SetTitle(ss.str().c_str());

	//calculate acceptance
	m_hTAcceptance = static_cast<TH1*>(m_hDistributions.at(3)->Clone());

	ss.clear();
	ss << "ALL acceptance " << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
			m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second;
	m_hTAcceptance->SetTitle(ss.str().c_str());

	m_hTAcceptance->Divide(m_hDistributions.at(2));

	//apply acceptance
	m_hTSlope->Divide(m_hTAcceptance);

	//get cumulative cross-sections
	double ALLCS = m_hTSlope->Integral() / totalLumiALL;
	double BHCS = m_hDistributions.at(1)->Integral() / totalLumiBH;

	//subtract BH
	m_hTSlope->Scale(1. / totalLumiALL);
	m_hTSlope->Add(m_hDistributions.at(1), -1 / totalLumiBH);

	std::cout << "debug: " << __func__ << ":" <<
		" BH: event fraction/total_lumi: " << m_lumiBH/totalLumiBH << "/" << totalLumiBH << 
		" ALL: event fraction/total_lumi: " << m_lumiALL/totalLumiALL << "/" << totalLumiALL << 
		" BH fraction: " << BHCS/ALLCS << std::endl;

	//radiative correction
	m_hTRC = static_cast<TH1*>(m_hDistributions.at(2)->Clone());

	ss.clear();
	ss << "DVCS RC " << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
			m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second;
	m_hTRC->SetTitle(ss.str().c_str());

	m_hTRC->Scale(1. / totalLumiALL); //TODO DIFFERENT LUMIS FOR TRUE AND BORN?
	m_hTRC->Add(m_hDistributions.at(1), -1 / totalLumiBH);

	TH1* tmpH = static_cast<TH1*>(m_hDistributions.at(5)->Clone());
	tmpH->Scale(1. / totalLumiALL); //TODO DIFFERENT LUMIS FOR TRUE AND BORN?
	tmpH->Add(m_hDistributions.at(4), -1 / totalLumiBH);

	m_hTRC->Divide(tmpH);

	//fit
	//fit options: 
	//0: do not attempt to draw function

	int statusCode = int(m_hTSlope->Fit(m_fFit, "0IME", "", 0.1, 1.));
	//int statusCode = int(m_hTSlope->Fit(m_fFit, "I"));

	//store results
	if(m_fitResult){

		delete m_fitResult;
		m_fitResult = nullptr;
	}

	m_fitResult = new FitResult();

	m_fitResult->setStatusCode(statusCode);
	m_fitResult->setNPoints(m_fFit->GetNumberFitPoints());
	m_fitResult->setChi2(m_fFit->GetChisquare());

	for(size_t i = 0; i < m_fFit->GetNpar(); i++){
		m_fitResult->addParameter(std::make_pair(m_fFit->GetParameter(i), m_fFit->GetParError(i)));
	}

	//print
	// std::cout << "NUCLEON_TOMOGRAPHY{" << getMeanXB() << ", " << getMeanQ2() << ", " << m_hTSlope->GetNbinsX() << ", " << m_hTSlope->GetXaxis()->GetXmin() << ", " << m_hTSlope->GetXaxis()->GetXmax();
	for(size_t i = 1; i <= m_hTSlope->GetNbinsX(); i++){
		std::cout << ", " << m_hTSlope->GetBinContent(i) << ", " << m_hTSlope->GetBinError(i);
	}
	std::cout << "}," << std::endl;
}

void BinTSlope::print() const{

	//run for parent class
	Bin::print();

	//print
	std::cout << __func__ << " debug: " << 
		"range xB: min: " << m_rangeXB.first << " max: " << m_rangeXB.second << " mean (from events): " << 
		getMeanXB(KinematicsType::Observed) << ' ' << getMeanXB(KinematicsType::True) << ' ' << getMeanXB(KinematicsType::Born) << std::endl;
	std::cout << __func__ << " debug: " << 
		"range Q2: min: " << m_rangeQ2.first << " max: " << m_rangeQ2.second << " mean (from events): " << 
		getMeanQ2(KinematicsType::Observed) << ' ' << getMeanQ2(KinematicsType::True) << ' ' << getMeanQ2(KinematicsType::Born) << std::endl;
}

const std::vector<TH1*>& BinTSlope::getHDistributions() const{
	return m_hDistributions;
}

TH1* BinTSlope::getHTSlope() const{
	return m_hTSlope;
}

TH1* BinTSlope::getHTAcceptance() const{
	return m_hTAcceptance;
}

TH1* BinTSlope::getHTRC() const{
	return m_hTRC;
}
