#include "../../include/bin/BinTSlope.h"

#include <iostream>
#include <sstream>

#include "../../include/other/HashManager.h"

BinTSlope::BinTSlope(
	const std::pair<double, double>& rangeXB, 
	const std::pair<double, double>& rangeQ2, 
	size_t nTBins,
	const std::pair<double, double>& rangeT
) : Bin("BinTSlope") {

	//reset
	reset();

	//set ranges
	m_rangeXB = checkRange(rangeXB);
	m_rangeQ2 = checkRange(rangeQ2);
	m_rangeT = checkRange(rangeT);

	//labels
	std::stringstream ss;

	ss << m_rangeXB.first << " #leq xB < " << m_rangeXB.second << " " <<
		m_rangeQ2.first << " #leq Q2 < " << m_rangeQ2.second;
 
	//make histograms
	m_hDistribution = new TH1D((HashManager::getInstance()->getHash()).c_str(), ss.str().c_str(), 
					nTBins, rangeT.first, rangeT.second);

	//set sumw2
	m_hDistribution->Sumw2();

	//function for fitting
	m_fFit = new TF1((HashManager::getInstance()->getHash()).c_str(), "[0]*exp(-1*[1]*x)", 0., 2.);

	m_fFit->SetParameter(0, m_sumWeights);
	m_fFit->SetParameter(1, 5.);

}

BinTSlope::~BinTSlope(){
}

void BinTSlope::reset(){

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
	m_hDistribution = nullptr;
	m_hTSlope = nullptr;
}

void BinTSlope::fill(DVCSEvent& event, double weight){

	//run for parent class
	Bin::fill(event, weight);

	//fill
	m_hDistribution->Fill(-1 * event.getT(), weight);
	
	//add
	m_sumXB += weight * event.getXB();
	m_sumQ2 += weight * event.getQ2();
	m_sumT += weight * event.getT();
}

void BinTSlope::analyse(){
	
	//run for parent class
	Bin::analyse();

	//skip bins with low entry
	if(m_nEvents < 500){
		return;
	}

	//skip bins with low summed weights
	if(m_sumWeights < 500.){
		return;
	}

	//make t-slope histogram
	m_hTSlope = static_cast<TH1*>(m_hDistribution->Clone());

	//fit
	//fit options: 
	//0: do not attempt to draw function

	int statusCode = int(m_hTSlope->Fit(m_fFit, "0B"));

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

void BinTSlope::print() const{

	//run for parent class
	Bin::print();

	//print
	std::cout << getClassName() << "::" << __func__ << " debug: " << 
		"range xB: min: " << m_rangeXB.first << " max: " << m_rangeXB.second << " mean (from events): " << getMeanXB() << std::endl;
	std::cout << getClassName() << "::" << __func__ << " debug: " << 
		"range Q2: min: " << m_rangeQ2.first << " max: " << m_rangeQ2.second << " mean (from events): " << getMeanQ2() << std::endl;
}

const std::pair<double, double>& BinTSlope::getRangeXB() const{
	return m_rangeXB;
}

const std::pair<double, double>& BinTSlope::getRangeQ2() const{
	return m_rangeQ2;
}

const std::pair<double, double>& BinTSlope::getRangeT() const{
	return m_rangeT;
}

double BinTSlope::getMeanXB() const{
	return getMean(m_sumXB, m_sumWeights);
}

double BinTSlope::getMeanQ2() const{
	return getMean(m_sumQ2, m_sumWeights);
}

double BinTSlope::getMeanT() const{
	return getMean(m_sumT, m_sumWeights);
}

TH1* BinTSlope::getHDistribution() const{
	return m_hDistribution;
}

TH1* BinTSlope::getHTSlope() const{
	return m_hTSlope;
}
