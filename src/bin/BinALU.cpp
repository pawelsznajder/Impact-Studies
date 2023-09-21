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
) : Bin("BinALU") {

	//reset
	reset();

	//set ranges
	m_rangeXB = checkRange(rangeXB);
	m_rangeQ2 = checkRange(rangeQ2);
	m_rangeT = checkRange(rangeT);
	m_rangePhi = checkRange(rangePhi);

	//labels
	std::stringstream ss;

	//make histograms
	ss.clear();
	ss << "observed: " << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
		m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second;
	
	m_hDistributions = std::make_pair(
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second), 
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second)
	);

	ss.clear();
	ss << "true: " << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
		m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second;

	m_hDistributionsTrue = std::make_pair(
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second), 
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second)
	);

	ss.clear();
	ss << "born: " << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
		m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second;

	m_hDistributionsBorn = std::make_pair(
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second), 
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second)
	);

	//set sumw2
	m_hDistributions.first->Sumw2();
	m_hDistributions.second->Sumw2();

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

	//reset ranges
	m_rangeXB = std::make_pair(0., 0.);
	m_rangeQ2 = std::make_pair(0., 0.);
	m_rangeT = std::make_pair(0., 0.);

	//reset sums
	m_sumXB = 0.;
	m_sumQ2 = 0.;
	m_sumT = 0.;

	//reset number of events
	m_nEvents = 0;

	//reset histograms
	m_hDistributions = std::make_pair(nullptr, nullptr);
	m_hDistributionsTrue = std::make_pair(nullptr, nullptr);
	m_hDistributionsBorn = std::make_pair(nullptr, nullptr);
	m_hDistributionsRC = std::make_pair(nullptr, nullptr);
	m_hAsymmetry = nullptr;
}

void BinALU::fill(DVCSEvent& event, double weight){

	
	//only reconstructed
	if(! event.isReconstructed()) return;

	//fill
	switch(event.getBeamPolarisation()){

	 	case -1:{

	 		if(weight > 0.) m_hDistributions.first->Fill(event.getPhi());
	 		m_hDistributionsTrue.first->Fill(event.getPhi(KinematicsType::True));
	 		m_hDistributionsBorn.first->Fill(event.getPhi(KinematicsType::Born));

	 		break;
	 	}

	   	case 1:{

	   		if(weight > 0.) m_hDistributions.second->Fill(event.getPhi());
	   		m_hDistributionsTrue.second->Fill(event.getPhi(KinematicsType::True));
	   		m_hDistributionsBorn.second->Fill(event.getPhi(KinematicsType::Born));

	   	 	break;
	   	}

	   default:{
	   		std::cout << __func__ << " error: wrong beam polarisation state, " << event.getBeamPolarisation() << std::endl;
			exit(0);
	   }
	}

	//kinematics
	if(weight > 0.){

		Bin::fill(event, weight);

		m_sumXB += weight * event.getXB();
		m_sumQ2 += weight * event.getQ2();
		m_sumT += weight * event.getT();
		m_sumPhi += weight * event.getPhi();
	}
}

void BinALU::analyse(){

	//print
	std::cout << "=================================================================" << std::endl;
	
	//run for parent class
	Bin::analyse();

	//skip bins with low entry count
	if(m_nEvents < 100){
		return;
	}

	//calculate RC
	m_hDistributionsRC.first = static_cast<TH1*>(m_hDistributionsBorn.first->Clone());
	m_hDistributionsRC.second = static_cast<TH1*>(m_hDistributionsBorn.second->Clone());

	m_hDistributionsRC.first->Divide(m_hDistributionsTrue.first);
	m_hDistributionsRC.second->Divide(m_hDistributionsTrue.second);

	//make asymmetry histogram
	m_hAsymmetry = 
		m_hDistributions.first->GetAsymmetry(m_hDistributions.second);

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

	//print
	std::cout << __func__ << " debug: " << 
		"range xB: min: " << m_rangeXB.first << " max: " << m_rangeXB.second << " mean (from events): " << getMeanXB() << std::endl;
	std::cout << __func__ << " debug: " << 
		"range Q2: min: " << m_rangeQ2.first << " max: " << m_rangeQ2.second << " mean (from events): " << getMeanQ2() << std::endl;
	std::cout << __func__ << 
		"range |t|: min: " << m_rangeT.first << " max: " << m_rangeT.second << " mean (from events): " << getMeanT() << std::endl;
}

const std::pair<double, double>& BinALU::getRangeXB() const{
	return m_rangeXB;
}

const std::pair<double, double>& BinALU::getRangeQ2() const{
	return m_rangeQ2;
}

const std::pair<double, double>& BinALU::getRangeT() const{
	return m_rangeT;
}

const std::pair<double, double>& BinALU::getRangePhi() const{
	return m_rangePhi;
}

double BinALU::getMeanXB() const{
	return getMean(m_sumXB, m_sumWeights);
}

double BinALU::getMeanQ2() const{
	return getMean(m_sumQ2, m_sumWeights);
}

double BinALU::getMeanT() const{
	return getMean(m_sumT, m_sumWeights);
}

double BinALU::getMeanPhi() const{
	return getMean(m_sumPhi, m_sumWeights);
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributions() const{
	return m_hDistributions;
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributionsTrue() const{
	return m_hDistributionsTrue;
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributionsBorn() const{
	return m_hDistributionsBorn;
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributionsRC() const{
	return m_hDistributionsRC;
}

TH1* BinALU::getHAsymmetry() const{
	return m_hAsymmetry;
}
