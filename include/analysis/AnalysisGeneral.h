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
	AnalysisGeneral(double targetLuminosity);

	//destructor
	virtual ~AnalysisGeneral();

	//fill with events
	virtual void fill(DVCSEvent& event, double weight);

	//perform analysis
	virtual void analyse();

	//plot histograms to file
	virtual void plot(const std::string& path);

private:

	//evaluate acceptance
	TH1* evaluateAcceptance(TH1** h, bool isAcceptance) const;
	TH2* evaluateAcceptance(TH2** h, bool isAcceptance) const;

	//luminosity corresponding to accumulated luminosity
	double m_lumi;

	//reconstruction probabilities
	double m_resProbPOut[3];
	double m_resProbEOut[3];
	double m_resProbGOut[3];
	double m_resProbExcl[3];
	double m_resCut[2];

	//2D histograms
	TH2* m_hxBvsQ2[5];
	TH2* m_hxBvsTOverQ2[5];
	TH2* m_hXBvsT[5];
	TH2* m_hXBvsY[5];
	TH2* m_hQ2vsT[5];
	TH2* m_hYvsT[5];
	TH2* m_hQ2vsY[5];

	//1D histograms
	TH1* m_hXB[5];
	TH1* m_hQ2[5];
	TH1* m_hT[5];
	TH1* m_hPhi[5];
	TH1* m_hPhiS[5];
	TH1* m_hY[5];
	TH1* m_hEtaEOut[5];
	TH1* m_hEtaPOut[5];
	TH1* m_hEtaGOut[5];
	TH1* m_hPPOut[5];
	TH1* m_hPPtOut[5];
	TH1* m_hPThOut[5];
	TH1* m_hPPhOut[5];
	TH1* m_hEPOut[5];
	TH1* m_hEPtOut[5];
	TH1* m_hEThOut[5];
	TH1* m_hEPhOut[5];
	TH1* m_hGPOut[5];
	TH1* m_hGPtOut[5];
	TH1* m_hGThOut[5];
	TH1* m_hGPhOut[5];
			
	TLegend* leg[6];
};
#endif
