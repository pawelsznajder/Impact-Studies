#ifndef BIN_H
#define BIN_H

#include <cstring>
#include <TH1.h>

#include "../../include/other/BaseObject.h"
#include "../../include/other/FitResult.h"
#include "../../include/event/DVCSEvent.h"

/*
 * Parent class for those representing single bins.
 */
class Bin : public BaseObject{

public:

	//constructor
	Bin(const std::string& className, 
		const std::pair<double, double>& rangeXB, 
		const std::pair<double, double>& rangeQ2, 
		const std::pair<double, double>& rangeT,
		const std::pair<double, double>& rangePhi);

	//destructor
	virtual ~Bin();

	//reset
	virtual void reset();

	//store
	virtual void fill(DVCSEvent& event, double weight = 1.);

	//print
	virtual void print() const;

	//analyse
	virtual void analyse();

	//get fit result
	FitResult* getFitResult() const;

	//get mean xB
	double getMeanXB(KinematicsType::Type kinematicType = KinematicsType::Observed) const;

	//get mean Q2
	double getMeanQ2(KinematicsType::Type kinematicType = KinematicsType::Observed) const;

	//get mean t
	double getMeanT(KinematicsType::Type kinematicType = KinematicsType::Observed) const;

	//get mean phi
	double getMeanPhi(KinematicsType::Type kinematicType = KinematicsType::Observed) const;

	//get number of events
	size_t getNEvents(KinematicsType::Type kinematicType = KinematicsType::Observed) const;

	//get sum of weights
	double getSumWeights(KinematicsType::Type kinematicType = KinematicsType::Observed) const;

	//get range xB
	const std::pair<double, double>& getRangeXB() const;

	//get range Q2
	const std::pair<double, double>& getRangeQ2() const;

	//get range T
	const std::pair<double, double>& getRangeT() const;

	//get range phi
	const std::pair<double, double>& getRangePhi() const;

protected:

	//check if belongs to this bin for given kinematics
	bool belongsToThisBin(DVCSEvent& event, KinematicsType::Type kinematicType) const;

	//get mean
	double getMean(double sum, double sumOfWeights) const;

	//check range
	std::pair<double, double> checkRange(const std::pair<double, double>& range) const;

	//print histogram
	void printHistogram(TH1* h, const std::string& token) const;

	FitResult* m_fitResult; //fit result

	std::map<KinematicsType::Type, double> m_sumXB;		//sum of xB values
	std::map<KinematicsType::Type, double> m_sumQ2;		//sum of Q2 values
	std::map<KinematicsType::Type, double> m_sumT;		//sum of t values
	std::map<KinematicsType::Type, double> m_sumPhi;	//sum of phi values

	std::map<KinematicsType::Type, size_t> m_nEvents;	//number of events
	std::map<KinematicsType::Type, double> m_sumWeights;//sum of weights

	std::pair<double, double> m_rangeXB;	//xB range
	std::pair<double, double> m_rangeQ2;	//Q2 range
	std::pair<double, double> m_rangeT;		//t range
	std::pair<double, double> m_rangePhi;	//phi range
};

#endif