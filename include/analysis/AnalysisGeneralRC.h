#ifndef ANALYSIS_GENERALRC_H
#define ANALYSIS_GENERALRC_H

#include <TH2.h>
#include <TH1.h>

#include "../include/analysis/Analysis.h"
#include "../include/event/DVCSEvent.h"

/*
 * Class for general analysis.
 */
class AnalysisGeneralRC : public Analysis{

public:

	//constructor
	AnalysisGeneralRC();

	//destructor
	virtual ~AnalysisGeneralRC();

	//fill with events
	virtual void fill(DVCSEvent& event, double weight);

	//perform analysis
	virtual void analyse();

	//plot histograms to file
	virtual void plot(const std::string& path);

private:

	//2D histograms
	TH2* m_hXBvsXB;
	TH2* m_hQ2vsQ2;
	TH2* m_hTvsT;
	TH2* m_hPhivsPhi;
	TH2* m_hPhiSvsPhiS;
	TH2* m_hYvsY;

};
#endif
