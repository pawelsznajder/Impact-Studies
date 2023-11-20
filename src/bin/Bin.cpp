#include "../../include/bin/Bin.h"

#include <iostream>
#include <limits>

Bin::Bin(
		const std::string& className,
		const std::pair<double, double>& rangeXB, 
		const std::pair<double, double>& rangeQ2, 
		const std::pair<double, double>& rangeT,
		const std::pair<double, double>& rangePhi) : BaseObject(className), 
	m_fitResult(nullptr){

	//reset
	reset();

	//set ranges
	m_rangeXB = checkRange(rangeXB);
	m_rangeQ2 = checkRange(rangeQ2);
	m_rangeT = checkRange(rangeT);
	m_rangePhi = checkRange(rangePhi);
}

Bin::~Bin(){
	reset();
}

void Bin::reset(){

	if(m_fitResult){

		delete m_fitResult;
		m_fitResult = nullptr;
	}

	//reset sums
	m_sumXB.clear();
	m_sumQ2.clear();
	m_sumT.clear();
	m_sumPhi.clear();

	m_nEvents.clear();
	m_sumWeights.clear();

	for(size_t i = 0; i < 3; i++){

		KinematicsType::Type kinematicsType;

		switch(i){

		case 0:{

			kinematicsType = KinematicsType::Observed;
			break;
		}

		case 1:{
			kinematicsType = KinematicsType::True;
			break;
		}

		case 2:{
			kinematicsType = KinematicsType::Born;
			break;
		}
		}

		m_sumXB.insert(std::make_pair(kinematicsType, 0.));
		m_sumQ2.insert(std::make_pair(kinematicsType, 0.));
		m_sumT.insert(std::make_pair(kinematicsType, 0.));
		m_sumPhi.insert(std::make_pair(kinematicsType, 0.));

		m_nEvents.insert(std::make_pair(kinematicsType, 0));
		m_sumWeights.insert(std::make_pair(kinematicsType, 0.));
	}
}

void Bin::fill(DVCSEvent& event, double weight){

	for(size_t i = 0; i < 3; i++){

		KinematicsType::Type kinematicsType;

		switch(i){

		case 0:{

			if(! event.isReconstructed()) continue;

			kinematicsType = KinematicsType::Observed;
			break;
		}

		case 1:{
			kinematicsType = KinematicsType::True;
			break;
		}

		case 2:{
			kinematicsType = KinematicsType::Born;
			break;
		}
		}

		if(belongsToThisBin(event, kinematicsType)){

			m_sumXB.find(kinematicsType)->second += weight * event.getXB(kinematicsType);
			m_sumQ2.find(kinematicsType)->second += weight * event.getQ2(kinematicsType);
			m_sumT.find(kinematicsType)->second += weight * fabs(event.getT(kinematicsType));
			m_sumPhi.find(kinematicsType)->second += weight * event.getPhi(kinematicsType);

			m_nEvents.find(kinematicsType)->second ++;
			m_sumWeights.find(kinematicsType)->second += weight; 
		}
	}
}

void Bin::analyse(){
}

bool Bin::belongsToThisBin(DVCSEvent& event, KinematicsType::Type kinematicType) const{

	double xB = event.getXB(kinematicType);
	double Q2 = event.getQ2(kinematicType);
	double mT = fabs(event.getT(kinematicType));
	double phi = event.getPhi(kinematicType);

	if(xB < m_rangeXB.first || xB >= m_rangeXB.second) return false; 
	if(Q2 < m_rangeQ2.first || Q2 >= m_rangeQ2.second) return false; 
	if(mT < m_rangeT.first || mT >= m_rangeT.second) return false; 
	if(phi < m_rangePhi.first || phi >= m_rangePhi.second) return false; 

	return true;
}

double Bin::getMean(double sum, double sumOfWeights) const{

	//check if zero
	if(sumOfWeights == 0.){

		// std::cout << getClassName() << "::" << __func__ << " warning: " << "sum of weight is zero" << std::endl;
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

	std::cout << __func__ << " debug: " << 
		"range xB: min: " << m_rangeXB.first << " max: " << m_rangeXB.second << " mean (observed/true/born): " << 
		getMeanXB(KinematicsType::Observed) << '/' << getMeanXB(KinematicsType::True) << '/' << getMeanXB(KinematicsType::Born) << std::endl;
	std::cout << __func__ << " debug: " << 
		"range Q2: min: " << m_rangeQ2.first << " max: " << m_rangeQ2.second << " mean (observed/true/born): " << 
		getMeanQ2(KinematicsType::Observed) << '/' << getMeanQ2(KinematicsType::True) << '/' << getMeanQ2(KinematicsType::Born) << std::endl;
	std::cout << __func__ << " debug: " << 
		"range |t|: min: " << m_rangeT.first << " max: " << m_rangeT.second << " mean (observed/true/born): " << 
		getMeanT(KinematicsType::Observed) << '/' << getMeanT(KinematicsType::True) << '/' << getMeanT(KinematicsType::Born) << std::endl;
	std::cout << __func__ << " debug: " << 
		"range phi: min: " << m_rangePhi.first << " max: " << m_rangePhi.second << " mean (observed/true/born): " << 
		getMeanPhi(KinematicsType::Observed) << '/' << getMeanPhi(KinematicsType::True) << '/' << getMeanPhi(KinematicsType::Born) << std::endl;
	std::cout << __func__ << " debug: number of events (observed/true/born): " << 
		getNEvents(KinematicsType::Observed) << '/' << getNEvents(KinematicsType::True) << '/' << getNEvents(KinematicsType::Born) << std::endl;
	std::cout << __func__ << " debug: sum of weights (observed/true/born): " << 
		getSumWeights(KinematicsType::Observed) << '/' << getSumWeights(KinematicsType::True) << '/' << getSumWeights(KinematicsType::Born) << std::endl;
}

FitResult* Bin::getFitResult() const{
	return m_fitResult;
}

double Bin::getMeanXB(KinematicsType::Type kinematicType) const{
	return getMean(m_sumXB.find(kinematicType)->second, m_sumWeights.find(kinematicType)->second);
}

double Bin::getMeanQ2(KinematicsType::Type kinematicType) const{
	return getMean(m_sumQ2.find(kinematicType)->second, m_sumWeights.find(kinematicType)->second);
}

double Bin::getMeanT(KinematicsType::Type kinematicType) const{
	return getMean(m_sumT.find(kinematicType)->second, m_sumWeights.find(kinematicType)->second);
}

double Bin::getMeanPhi(KinematicsType::Type kinematicType) const{
	return getMean(m_sumPhi.find(kinematicType)->second, m_sumWeights.find(kinematicType)->second);
}

size_t Bin::getNEvents(KinematicsType::Type kinematicType) const{
	return m_nEvents.find(kinematicType)->second;
}

double Bin::getSumWeights(KinematicsType::Type kinematicType) const{
	return m_sumWeights.find(kinematicType)->second;
}

const std::pair<double, double>& Bin::getRangeXB() const{
	return m_rangeXB;
}

const std::pair<double, double>& Bin::getRangeQ2() const{
	return m_rangeQ2;
}

const std::pair<double, double>& Bin::getRangeT() const{
	return m_rangeT;
}

const std::pair<double, double>& Bin::getRangePhi() const{
	return m_rangePhi;
}