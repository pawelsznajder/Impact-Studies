#ifndef ANALYSIS_GENERAL_H
#define ANALYSIS_GENERAL_H

#include <TH2.h>
#include <TH1.h>

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
	virtual void fill(DVCSEvent& event, double weight);

	//perform analysis
	virtual void analyse();

	//plot histograms to file
	virtual void plot(const std::string& path);

private:

	//examplary histogram
	TH2* m_hXBvsQ2;

	// 1D histogram
	TH1* m_hXB;
	TH1* m_hQ2;
	TH1* m_hT;
	TH1* m_hPhi;
	TH1* m_hPhiS;
	TH1* m_hY;
	TH1* m_hEtaEOut;
	TH1* m_hEtaPOut;
	TH1* m_hEtaGOut;
};
#endif
