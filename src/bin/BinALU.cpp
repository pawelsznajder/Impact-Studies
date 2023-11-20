#include "../../include/bin/BinALU.h"

#include <iostream>
#include <sstream>

#include "../../include/other/HashManager.h"
#include "../../include/other/SubProcessType.h"

BinALU::BinALU(
	const std::pair<double, double>& rangeXB, 
	const std::pair<double, double>& rangeQ2, 
	const std::pair<double, double>& rangeT,
	size_t nPhiBins, 
	const std::pair<double, double>& rangePhi
) : Bin("BinALU", rangeXB, rangeQ2, rangeT, rangePhi) {

	//reset
	reset();

	//labels
	std::stringstream ss;

	//make histograms
	ss.clear();
	ss << "observed (lumi): " << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
		m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second << " " <<
		m_rangeT.first << " #leq |t| < " << m_rangeT.second;
	
	m_hDistributions = std::make_pair(
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second), 
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second)
	);

	ss.clear();
	ss << "observed: " << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
		m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second << " " <<
		m_rangeT.first << " #leq |t| < " << m_rangeT.second;
	
	m_hDistributionsObserved = std::make_pair(
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second), 
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second)
	);

	ss.clear();
	ss << "true: " << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
		m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second << " " <<
		m_rangeT.first << " #leq |t| < " << m_rangeT.second;

	m_hDistributionsTrue = std::make_pair(
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second), 
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second)
	);

	ss.clear();
	ss << "born: " << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
		m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second << " " <<
		m_rangeT.first << " #leq |t| < " << m_rangeT.second;

	m_hDistributionsBorn = std::make_pair(
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second), 
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second)
	);

	//set sumw2
	m_hDistributions.first->Sumw2();
	m_hDistributions.second->Sumw2();

	m_hDistributionsObserved.first->Sumw2();
	m_hDistributionsObserved.second->Sumw2();

	m_hDistributionsTrue.first->Sumw2();
	m_hDistributionsTrue.second->Sumw2();

	m_hDistributionsBorn.first->Sumw2();
	m_hDistributionsBorn.second->Sumw2();

	//function for fitting
	m_fFit = new TF1((HashManager::getInstance()->getHash()).c_str(), "[0]*sin(x)", 0., 2 * M_PI);
}

BinALU::~BinALU(){
}

void BinALU::reset(){

	//run for parent class
	Bin::reset();

	//reset histograms
	m_hDistributions = std::make_pair(nullptr, nullptr);
	m_hDistributionsObserved = std::make_pair(nullptr, nullptr);
	m_hDistributionsTrue = std::make_pair(nullptr, nullptr);
	m_hDistributionsBorn = std::make_pair(nullptr, nullptr);
	m_hDistributionsAcceptance = std::make_pair(nullptr, nullptr);
	m_hDistributionsRC = std::make_pair(nullptr, nullptr);
	m_hDistributionsCorrected = std::make_pair(nullptr, nullptr);
	m_hAsymmetry = nullptr;
	m_hSum = nullptr;
	m_hDif = nullptr;
}

void BinALU::fill(DVCSEvent& event, double weight){

	//get beam polarisation state
	const int& beamPolarisation = event.getBeamPolarisation();

	//fill
	if(weight > 0.){

		//kinematics
		Bin::fill(event);

		//fill
		if(event.isReconstructed() && belongsToThisBin(event, KinematicsType::Observed)) 
			getH(m_hDistributions, beamPolarisation)->Fill(event.getPhi(KinematicsType::Observed));

	}

	//fill
	if(event.isReconstructed() && belongsToThisBin(event, KinematicsType::Observed)) 
		getH(m_hDistributionsObserved, beamPolarisation)->Fill(event.getPhi(KinematicsType::Observed));
	if(belongsToThisBin(event, KinematicsType::True)) getH(m_hDistributionsTrue, beamPolarisation)->Fill(event.getPhi(KinematicsType::True));
	if(belongsToThisBin(event, KinematicsType::Born)) getH(m_hDistributionsBorn, beamPolarisation)->Fill(event.getPhi(KinematicsType::Born));
}

void BinALU::analyse(){
	
	//run for parent class
	Bin::analyse();

	//skip bins with low entry count
	if(getNEvents(KinematicsType::Observed) < 500){
		return;
	}

	//calculate acceptance
	m_hDistributionsAcceptance.first = static_cast<TH1*>(m_hDistributionsObserved.first->Clone());
	m_hDistributionsAcceptance.second = static_cast<TH1*>(m_hDistributionsObserved.second->Clone());

	m_hDistributionsAcceptance.first->Divide(m_hDistributionsTrue.first);
	m_hDistributionsAcceptance.second->Divide(m_hDistributionsTrue.second);

	//calculate RC
	m_hDistributionsRC.first = static_cast<TH1*>(m_hDistributionsTrue.first->Clone());
	m_hDistributionsRC.second = static_cast<TH1*>(m_hDistributionsTrue.second->Clone());

	m_hDistributionsRC.first->Divide(m_hDistributionsBorn.first);
	m_hDistributionsRC.second->Divide(m_hDistributionsBorn.second);

	//get corrected distributions
	m_hDistributionsCorrected.first = static_cast<TH1*>(m_hDistributions.first->Clone());
	m_hDistributionsCorrected.second = static_cast<TH1*>(m_hDistributions.second->Clone());

	m_hDistributionsCorrected.first->Divide(m_hDistributionsAcceptance.first);
	m_hDistributionsCorrected.second->Divide(m_hDistributionsAcceptance.second);

	m_hDistributionsCorrected.first->Divide(m_hDistributionsRC.first);
	m_hDistributionsCorrected.second->Divide(m_hDistributionsRC.second);

	//sum and difference
	m_hSum = static_cast<TH1*>(m_hDistributionsCorrected.first->Clone());
	m_hSum->Add(m_hDistributionsCorrected.second);

	m_hDif = static_cast<TH1*>(m_hDistributionsCorrected.first->Clone());
	m_hDif->Add(m_hDistributionsCorrected.second, -1.);

	printHistogram(m_hSum, "ALUSUM");
	printHistogram(m_hDif, "ALUDIF");

	//make asymmetry histogram
	m_hAsymmetry = 
		m_hDistributions.first->GetAsymmetry(m_hDistributions.second);

	printHistogram(m_hAsymmetry, "ALUASS");

	//fit
	//fit options: 
	//				0: do not attempt to draw function
	int statusCode = int(m_hAsymmetry->Fit(m_fFit, "0"));

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
}

void BinALU::print() const{

	//run for parent class
	Bin::print();
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributions() const{
	return m_hDistributions;
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributionsObserved() const{
	return m_hDistributionsObserved;
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributionsTrue() const{
	return m_hDistributionsTrue;
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributionsBorn() const{
	return m_hDistributionsBorn;
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributionsAcceptance() const{
	return m_hDistributionsAcceptance;
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributionsRC() const{
	return m_hDistributionsRC;
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributionsCorrected() const{
	return m_hDistributionsCorrected;
}

TH1* BinALU::getHSum() const{
	return m_hSum;
}

TH1* BinALU::getHDifference() const{
	return m_hDif;
}

TH1* BinALU::getHAsymmetry() const{
	return m_hAsymmetry;
}

TH1* BinALU::getH(const std::pair<TH1*, TH1*>& histogramPair, int beamPolarisation) const{

	switch(beamPolarisation){

	 	case -1:{
	 		return histogramPair.first;
	 	}

	   	case 1:{
	 		return histogramPair.second;
	   	}

	   default:{
	   		std::cout << __func__ << " error: wrong beam polarisation state, " << beamPolarisation << std::endl;
			exit(0);
	   }
	}
}