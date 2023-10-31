#include "../../include/bin/Bin.h"

#include <iostream>
#include <limits>

Bin::Bin(const std::string& className) : BaseObject(className), 
	m_nEvents(0), m_sumWeights(0.), m_fitResult(nullptr){
}

Bin::~Bin(){
	reset();
}

void Bin::reset(){

	m_nEvents = 0;
	m_sumWeights = 0.;

	if(m_fitResult){

		delete m_fitResult;
		m_fitResult = nullptr;
	}
}

void Bin::fill(DVCSEvent& event, double weight){

	//add
	m_nEvents++;
	m_sumWeights += weight;
}

void Bin::analyse(){
}

double Bin::getMean(double sum, double sumOfWeights) const{

	//check if zero
	if(sumOfWeights == 0.){

		std::cout << getClassName() << "::" << __func__ << " warning: " << "sum of weight is zero" << std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	}

	//return
	return sum / sumOfWeights;
}

std::pair<double, double> Bin::checkRange(const std::pair<double, double>& range) const{

	//check if equal
	if(range.first == range.second){

		std::cout << getClassName() << "::" << __func__ << " error: " << "invalid range, (" <<
			range.first << ", " << range.second << ")" << std::endl;
		exit(0);

		return std::make_pair(0., 0.);
	}

	//check if in order
	if(range.first > range.second){

		std::cout << getClassName() << "::" << __func__ << " warning: " << "inverted range, (" <<
			range.first << ", " << range.second << "), has been corrected" << std::endl;

		return std::make_pair(range.second, range.first);
	}

	return range;
}

void Bin::printHistogram(TH1* h, const std::string& token) const{

	std::cout << token;

	for(size_t i = 1; i <= h->GetNbinsX(); i++){
		std::cout << ' ' << h->GetBinCenter(i) << ' ' << h->GetBinContent(i) << ' ' << h->GetBinError(i);
	}

	std::cout << std::endl;
}

void Bin::print() const{
	std::cout << getClassName() << "::" << __func__ << " debug: " << 
		"number of events: " << m_nEvents << std::endl;
}

size_t Bin::getNEvents() const{
	return m_nEvents;
}

double Bin::getSumWeights() const{
	return m_sumWeights;
}

FitResult* Bin::getFitResult() const{
	return m_fitResult;
}