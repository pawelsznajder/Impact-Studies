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
	Bin(const std::string& className);

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

	//get number of stored events
	size_t getNEvents() const;

	//get sum of weights
	double getSumWeights() const;

	//get fit result
	FitResult* getFitResult() const;

protected:

	//get mean
	double getMean(double sum, double sumOfWeights) const;

	//check range
	std::pair<double, double> checkRange(const std::pair<double, double>& range) const;

	//print histogram
	void printHistogram(TH1* h, const std::string& token) const;

	size_t m_nEvents;		//number of stored events
	double m_sumWeights;	//sum of weights

	FitResult* m_fitResult; //fit result
};

#endif