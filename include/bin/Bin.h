#ifndef BIN_H
#define BIN_H

#include <cstring>

#include "../../include/other/BaseObject.h"
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
	virtual void fill(const DVCSEvent& event, double weight);

	//analyse (make fit) 
	//return histogram and result and uncertainty pair
	virtual std::pair<double, double> analyse();

	//get number of stored events
	size_t getNEvents() const;

	//get sum of weights
	double getSumWeights() const;

protected:

	//get mean
	double getMean(double sum, size_t nEvents) const;

	//check range
	std::pair<double, double> checkRange(const std::pair<double, double>& range) const;

	size_t m_nEvents;		//number of stored events
	double m_sumWeights;	//sum of weights
};

#endif