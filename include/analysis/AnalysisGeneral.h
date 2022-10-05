#ifndef ANALYSIS_GENERAL_H
#define ANALYSIS_GENERAL_H

#include <TH2.h>
#include <TH1.h>
#include <TLegend.h>

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

	//reconstruction probabilities
	double m_resProbPOut[2];
	double m_resProbEOut[2];
	double m_resProbGOut[2];

	double m_resProbExcl[2];

	//2D histograms
	TH2* m_hXBvsQ2;
	TH2* m_hxBvsQ2[2];
	TH2* m_hXBvsT[2];
	TH2* m_hXBvsY[2];
	TH2* m_hQ2vsT[2];
	TH2* m_hYvsT[2];
	TH2* m_hQ2vsY[2];

	//1D histograms
	TH1* m_hXB[2];
	TH1* m_hQ2[2];
	TH1* m_hT[2];
	TH1* m_hPhi[2];
	TH1* m_hPhiS[2];
	TH1* m_hY[2];
	TH1* m_hEtaEOut[2];
	TH1* m_hEtaPOut[2];
	TH1* m_hEtaGOut[2];
	TH1* m_hPPOut[2];
	TH1* m_hPPtOut[2];
	TH1* m_hPThOut[2];
	TH1* m_hPPhOut[2];
	TH1* m_hEPOut[2];
	TH1* m_hEPtOut[2];
	TH1* m_hEThOut[2];
	TH1* m_hEPhOut[2];
	TH1* m_hGPOut[2];
	TH1* m_hGPtOut[2];
	TH1* m_hGThOut[2];
	TH1* m_hGPhOut[2];
	
	TH1* m_hRatio[21];
		
	TLegend* leg[6];
};
#endif
