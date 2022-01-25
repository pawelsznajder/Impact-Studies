#include "../../include/bin/Bin.h"

#include <iostream>

Bin::Bin(const std::string& className) : BaseObject(className){
	reset();
}

Bin::~Bin(){
	reset();
}

void Bin::reset(){
	m_nEvents = 0;
}

void Bin::fill(const DVCSEvent& event, double weight){

	//check weight
	if(weight < 0.){
		std::cout << getClassName() << "::" << __func__ << " warning: " 
			<< "event weight is negative, " << weight << std::endl;
	}

	//add
	m_nEvents++;
	m_sumWeights += weight;
}

FitResult Bin::analyse(){
	
	//return 
	return FitResult();
}

double Bin::getMean(double sum, double sumOfWeights) const{

	//check if zero
	if(sumOfWeights == 0.){

		std::cout << getClassName() << "::" << __func__ << " error: " << "sum of weight is zero or negative, "
			<< sumOfWeights << std::endl;
		exit(0);
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