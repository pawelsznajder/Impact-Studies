#ifndef BINTSLOPE_H
#define BINTSLOPE_H

#include <TH1.h>
#include <TF1.h>
#include <utility>

#include "../../include/bin/Bin.h"

/*
 * Single xB, Q2 bin.
 */
class BinTSlope : public Bin{

public:

	//constructor
	BinTSlope(
		const std::pair<double, double>& rangeXB, 
		const std::pair<double, double>& rangeQ2, 
		size_t nTBins,
		const std::pair<double, double>& rangeT
	);

	//destructor
	virtual ~BinTSlope();

	//reset
	virtual void reset();

	//store
	virtual void fill(DVCSEvent& event, double weight);

	//print
	virtual void print() const;

	//analyse
	virtual void analyse();

	//get range xB
	const std::pair<double, double>& getRangeXB() const;

	//get range Q2
	const std::pair<double, double>& getRangeQ2() const;

	//get range T
	const std::pair<double, double>& getRangeT() const;

	//get mean xB
	double getMeanXB() const;

	//get mean Q2
	double getMeanQ2() const;

	//get mean t
	double getMeanT() const;

	//get histograms to store events
	TH1* getHDistribution() const;

	//get histograms for resulting t-slope
	TH1* getHTSlope() const;

private:

	std::pair<double, double> m_rangeXB;	//xB range
	std::pair<double, double> m_rangeQ2;	//Q2 range
	std::pair<double, double> m_rangeT;		//t range

	double m_sumXB;		//sum of xB values.
	double m_sumQ2;		//sum of Q2 values.
	double m_sumT;		//sum of t values.

	TH1* m_hDistribution;	//histograms to store events
	TH1* m_hTSlope;			//histogram for resulting t-slope
	TF1* m_fFit;			//function for fitting 
};

#endif