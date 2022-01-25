#include "../../include/bin/BinALU.h"

#include <iostream>
#include <sstream>

#include "../../include/other/HashManager.h"

BinALU::BinALU(
	const std::pair<double, double>& rangeXB, 
	const std::pair<double, double>& rangeQ2, 
	const std::pair<double, double>& rangeT,
	size_t nPhiBins, 
	const std::pair<double, double>& rangePhi
) : Bin("BinALU") {

	//set ranges
	m_rangeXB = checkRange(rangeXB);
	m_rangeQ2 = checkRange(rangeQ2);
	m_rangeT = checkRange(rangeT);
	m_rangePhi = checkRange(rangePhi);

	//labes
	std::stringstream ss;

	ss << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
		m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second;
 
	//make histograms
	m_hDistributions = std::make_pair(
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second), 
		new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
			nPhiBins, rangePhi.first, rangePhi.second)
	);

	//set sumw2
	m_hDistributions.first->Sumw2();
	m_hDistributions.second->Sumw2();
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
	m_hAsymmetry = nullptr;
}

void BinALU::fill(const DVCSEvent& event, double weight){

	//run for parent class
	Bin::fill(event, weight);

	//pol -
	if(event.getBeamPolarisation() == -1){
		m_hDistributions.first->Fill(event.getPhi(), weight);
	}

	//pol +
	else if(event.getBeamPolarisation() == 1){
		m_hDistributions.second->Fill(event.getPhi(), weight);
	}

	//none
	else{

		std::cout << getClassName() << "::" << __func__ << " error: " << 
			"wrong polarisation, " << event.getBeamPolarisation() << std::endl;
		exit(0);
	}
}

std::pair<double, double> BinALU::analyse(){
	
	//run for parent class
	Bin::analyse();

	//make asymmetry histogram
	m_hAsymmetry = 
		m_hDistributions.first->GetAsymmetry(m_hDistributions.second);
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

double BinALU::getMeanXB() const{
	return getMean(m_sumXB, m_nEvents);
}

double BinALU::getMeanQ2() const{
	return getMean(m_sumQ2, m_nEvents);
}

double BinALU::getMeanT() const{
	return getMean(m_sumT, m_nEvents);
}

const std::pair<TH1*, TH1*>& BinALU::getHDistributions() const{
	return m_hDistributions;
}

TH1* BinALU::getHAsymmetry() const{

	if(! m_isLocked){

		std::cout << getClassName() << "::" << __func__ << " error: " << 
			"no analysis made so far for this bin" << std::endl;
		exit(0);
	}

	return m_hAsymmetry;
}
