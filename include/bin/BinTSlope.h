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

	//get histograms to store events
	const std::vector<TH1*>& getHDistributions() const;

	//get histograms for resulting t-slope
	TH1* getHTSlope() const;

	//get histograms for acceptance
	TH1* getHTAcceptance() const;

	//get histograms for acceptance
	TH1* getHTRC() const;

private:

	std::vector<TH1*> m_hDistributions;	//histograms to store events
	TH1* m_hTSlope;			//histogram for resulting t-slope
	TH1* m_hTAcceptance;	//histogram for acceptance
	TF1* m_fFit;			//function for fitting 
};

#endif