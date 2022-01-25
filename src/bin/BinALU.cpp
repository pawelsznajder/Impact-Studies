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

	//reset
	reset();

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

	//add
	m_sumXB += weight * event.getXB();
	m_sumQ2 += weight * event.getQ2();
	m_sumT += weight * event.getT();
	m_sumPhi += weight * event.getPhi();
}

FitResult BinALU::analyse(){
	
	//run for parent class
	Bin::analyse();

	//skip bins with low entry
	if(m_nEvents < 100){
		return FitResult();
	}

	//make asymmetry histogram
	m_hAsymmetry = 
		m_hDistributions.first->GetAsymmetry(m_hDistributions.second);

	//fit
	//fit options: 
	//				0: do not attempt to draw function
	int statusCode = int(m_hAsymmetry->Fit(m_fFit, "0"));

	//store results
	FitResult result;

	result.setStatusCode(statusCode);
	result.setNPoints(m_fFit->GetNumberFitPoints());
	result.setChi2(m_fFit->GetChisquare());

	for(size_t i = 0; i < m_fFit->GetNpar(); i++){
		result.addParameter(std::make_pair(m_fFit->GetParameter(i), m_fFit->GetParError(i)));
	}

	//return
	return result;
}

void BinALU::print() const{

	//run for parent class
	Bin::print();

	//print
	std::cout << getClassName() << "::" << __func__ << " debug: " << 
		"range xB: min: " << m_rangeXB.first << " max: " << m_rangeXB.second << " mean (from events): " << getMeanXB() << std::endl;
	std::cout << getClassName() << "::" << __func__ << " debug: " << 
		"range Q2: min: " << m_rangeQ2.first << " max: " << m_rangeQ2.second << " mean (from events): " << getMeanQ2() << std::endl;
	std::cout << getClassName() << "::" << __func__ << " debug: " << 
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

TH1* BinALU::getHAsymmetry() const{
	return m_hAsymmetry;
}
