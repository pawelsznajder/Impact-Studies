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
	m_isLocked = false;
}

void Bin::fill(const DVCSEvent& event, double weight){

	//check weight
	if(weight < 0.){
		std::cout << getClassName() << "::" << __func__ << " warning: " 
			<< "event weight is negative, " << weight << std::endl;
	}

	//check if locked
	if(m_isLocked){

		std::cout << getClassName() << "::" << __func__ << " error: " 
			<< "bin was previously locked by fit() method" << std::endl;
		exit(0);
	}
}

std::pair<double, double> Bin::analyse(){

	//mark as locked	
	m_isLocked = true;
	
	//return 
	return std::make_pair(0., 0.);
}

double Bin::getMean(double sum, size_t nEvents) const{

	//check if zero
	if(nEvents == 0){

		std::cout << getClassName() << "::" << __func__ << " error: " << "number of events is zero" << std::endl;
		exit(0);
	}

	//return
	return sum / double(nEvents);
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