#ifndef ANALYSIS_TSlope_H
#define ANALYSIS_TSlope_H

#include "../include/analysis/Analysis.h"
#include "../include/event/DVCSEvent.h"
#include "../include/bin/BinTSlope.h"

/*
 * Class to extract t-slope.
 */
class AnalysisTSlope : public Analysis{

public:

	//constructor
	AnalysisTSlope(double targetLuminosity);

	//destructor
	virtual ~AnalysisTSlope();

	//fill with events
	virtual void fill(DVCSEvent& event, double weight);

	//perform analysis
	virtual void analyse();

	//plot histograms to file
	virtual void plot(const std::string& path);

private:

	//set bin boundaries
	void setBinBoundaries();

	//initialise bins
	void initialiseBins();

	//bin boundaries
	std::vector<double> m_binBoundsXB;
	std::vector<double> m_binBoundsQ2;

	//bin ranges
	std::vector<std::pair<double, double> > m_binRangesXB;
	std::vector<std::pair<double, double> > m_binRangesQ2;

	//number of t bins
	size_t m_nBinsT;

	//analysis bins
	std::vector<BinTSlope> m_bins;

	//integrated luminosity
	double m_lumiM, m_lumiP;
};
#endif