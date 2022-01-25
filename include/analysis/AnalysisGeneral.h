#ifndef ANALYSIS_GENERAL_H
#define ANALYSIS_GENERAL_H

#include <TH2.h>

#include "../include/analysis/Analysis.h"
#include "../include/event/DVCSEvent.h"

/*
 * Class for general analysis.
 */
class AnalysisGeneral : public Analysis{

public:

	//constructor
	AnalysisGeneral();

	//destructor
	virtual ~AnalysisGeneral();

	//fill with events
	virtual void fill(const DVCSEvent& event, double weight);

	//perform analysis
	virtual void analyse();

	//plot histograms to file
	virtual void plot(const std::string& path);

private:

	//examplary histogram
	TH2* m_hXBvsQ2;
};
#endif