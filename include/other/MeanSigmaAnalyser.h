#ifndef MEAN_SIGMA_ANALYSER_H
#define MEAN_SIGMA_ANALYSER_H

#include <vector>
#include <TH1.h>

#include "../../include/other/BaseObject.h"

/*
 * Container storing result of a single fit.
 */
class MeanSigmaAnalyser : public BaseObject{

public:

	//constructor
	MeanSigmaAnalyser(const std::string& title, size_t nBins, const std::pair<double, double>& range, size_t maxNEvents = 0);

	//constructor
	MeanSigmaAnalyser(TH1* h, size_t maxNEvents = 0);

	//destructor
	virtual ~MeanSigmaAnalyser();

	//fill
	void fill(double x, double value, double weight);

	//reset
	void reset();

	//analyse
	TH1* getH();

private:

	//get mean
	double getMean(size_t nEvents, const std::vector<std::pair<double, double> >& v) const;

	//get rmse
	double getRMSE(double mean, size_t nEvents, const std::vector<std::pair<double, double> >& v) const;

	std::string m_title;						//title
	size_t m_nBins;								//number of bins
	std::pair<double, double> m_range;			//range
	size_t m_maxNEvents;						//maximum number of events in a given bin

	std::vector<size_t> m_nEvents;				//number of events in bins
	std::vector<std::vector<std::pair<double, double> > > m_events;	//events
};

#endif