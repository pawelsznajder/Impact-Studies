#include "../../include/bin/BinTSlope.h"

#include <iostream>
#include <sstream>

#include "../../include/other/HashManager.h"
#include "../../include/kinematic_cuts/KinematicCuts.h"

BinTSlope::BinTSlope(
	const std::pair<double, double>& rangeXB, 
	const std::pair<double, double>& rangeQ2, 
	size_t nTBins,
	const std::pair<double, double>& rangeT
) : Bin("BinTSlope", rangeXB, rangeQ2, rangeT, std::make_pair(0., 2 * M_PI)) {

	//reset
	reset();
 
	//make histograms
 	for(size_t i = 0; i < 3; i++){

 		//labels
		std::stringstream ss;

		if(i == 0) ss << "DVCS rec lumi: ";
		if(i == 1) ss << "DVCS true: ";
		if(i == 2) ss << "DVCS rec: ";

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
}

void BinTSlope::fill(DVCSEvent& event, double weight){

	for(size_t i = 0; i < 2; i++){

		KinematicsType::Type kinType = ((i == 0)?(KinematicsType::True):(KinematicsType::Observed));

		if(KinematicCuts::checkKinematicCuts(event, kinType)){

			if(event.getXB(kinType) < getRangeXB().first || 
				event.getXB(kinType) >= getRangeXB().second) continue; 

			if(event.getQ2(kinType) < getRangeQ2().first || 
				event.getQ2(kinType) >= getRangeQ2().second) continue; 

			if(i == 0){
				m_hDistributions.at(1)->Fill(-1 * event.getT(kinType));
			}

			if(i == 1){
				if(weight > 0.){

					Bin::fill(event);

					m_hDistributions.at(0)->Fill(-1 * event.getT(kinType));
				}else{
					m_hDistributions.at(2)->Fill(-1 * event.getT(kinType));
				}
			}
		}
	}
}

void BinTSlope::analyse(){
	
	//run for parent class
	Bin::analyse();

	//skip bins with low entry count
	if(getNEvents(KinematicsType::Observed) < 2000){
		return;
	}

	for(size_t i = 0; i < 3; i++){
		std::cout << "histogram size: " << i << ": " << m_hDistributions.at(i)->GetEntries() << std::endl;
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
	m_hTAcceptance = static_cast<TH1*>(m_hDistributions.at(2)->Clone());

	ss.clear();
	ss << "ALL acceptance " << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
			m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second;
	m_hTAcceptance->SetTitle(ss.str().c_str());

	m_hTAcceptance->Divide(m_hDistributions.at(1));

	//apply acceptance
	m_hTSlope->Divide(m_hTAcceptance);

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
	std::cout << "NUCLEON_TOMOGRAPHY{" << getMeanXB() << ", " << getMeanQ2() << ", " << m_hTSlope->GetNbinsX() << ", " << m_hTSlope->GetXaxis()->GetXmin() << ", " << m_hTSlope->GetXaxis()->GetXmax();
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
