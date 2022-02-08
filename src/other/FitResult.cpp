#include "../../include/other/FitResult.h"

#include <iostream>

FitResult::FitResult() : BaseObject("FitResult"){

	m_statusCode = -10;
	m_nPoints = 0;
	m_chi2 = 0.;
}

FitResult::~FitResult(){
}

int FitResult::getStatusCode() const{
	return m_statusCode;
}		

size_t FitResult::getNPoints() const{
	return m_nPoints;
}		

double FitResult::getChi2() const{
	return m_chi2;
}	

double FitResult::getChi2Norm() const{
	return m_chi2 / (m_nPoints - m_parameters.size());
}		

void FitResult::setStatusCode(int statusCode){
	m_statusCode = statusCode;
}		 		

void FitResult::setNPoints(size_t nPoints){
	m_nPoints = nPoints;
}

void FitResult::setChi2(double chi2){
	m_chi2 = chi2;
}	

size_t FitResult::getNDF() const{
	return m_nPoints - m_parameters.size();
}

void FitResult::addParameter(const std::pair<double, double>& parameter){
	m_parameters.push_back(parameter);
}	

size_t FitResult::getNParameters() const{
	return m_parameters.size();
}	

const std::pair<double, double>& FitResult::getParameter(size_t i) const{

	if(! (i < m_parameters.size())){
		std::cout << getClassName() << "::" << __func__ << " error: " << 
			"wrong index, " << i << std::endl;
		exit(0);
	}

	return m_parameters.at(i);
}

void FitResult::print() const{

	std::cout << getClassName() << "::" << __func__ << " debug: " <<
		"status code: " << m_statusCode << " chi2/ndf: " << getChi2Norm() << std::endl;

	for(std::vector<std::pair<double, double> >::const_iterator it = m_parameters.begin(); it != m_parameters.end(); it++){
		std::cout << getClassName() << "::" << __func__ << " debug: " <<
			"par " << size_t(it - m_parameters.begin()) << ": value: " << it->first << " unc: " << it->second << std::endl;
	}
}
