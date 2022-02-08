#ifndef ANALYSIS_ALU_H
#define ANALYSIS_ALU_H

#include "../include/analysis/Analysis.h"
#include "../include/event/DVCSEvent.h"
#include "../include/bin/BinALU.h"

/*
 * Class to extract ALU asymmetries.
 */
class AnalysisALU : public Analysis{

public:

	//constructor
	AnalysisALU();

	//destructor
	virtual ~AnalysisALU();

	//fill with events
	virtual void fill(const DVCSEvent& event, double weight);

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
	std::vector<double> m_binBoundsT;

	//bin ranges
	std::vector<std::pair<double, double> > m_binRangesXB;
	std::vector<std::pair<double, double> > m_binRangesQ2;
	std::vector<std::pair<double, double> > m_binRangesT;

	//number of phi bins
	size_t m_nBinsPhi;

	//analysis bins
	std::vector<BinALU> m_bins;
};
#endif