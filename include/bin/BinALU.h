#ifndef BINALU_H
#define BINALU_H

#include <TH1.h>
#include <TF1.h>
#include <utility>

#include "../../include/bin/Bin.h"

/*
 * Single xB, Q2, t bin.
 */
class BinALU : public Bin{

public:

	//constructor
	BinALU(
		const std::pair<double, double>& rangeXB, 
		const std::pair<double, double>& rangeQ2, 
		const std::pair<double, double>& rangeT,
		size_t nPhiBins, 
		const std::pair<double, double>& rangePhi
	);

	//destructor
	virtual ~BinALU();

	//reset
	virtual void reset();

	//store
	virtual void fill(DVCSEvent& event, double weight);

	//print
	virtual void print() const;

	//analyse
	virtual void analyse();

	//get histograms to store events (Observed with lumi. limit)
	const std::pair<TH1*, TH1*>& getHDistributions() const;

	//get histograms to store events (Observed)
	const std::pair<TH1*, TH1*>& getHDistributionsObserved() const;

	//get histograms to store events (True)
	const std::pair<TH1*, TH1*>& getHDistributionsTrue() const;

	//get histograms to store events (Born)
	const std::pair<TH1*, TH1*>& getHDistributionsBorn() const;

	//get histograms to store events (Observed/True)
	const std::pair<TH1*, TH1*>& getHDistributionsAcceptance() const;

	//get histograms to store events (Born/True)
	const std::pair<TH1*, TH1*>& getHDistributionsRC() const;

	//get histograms to store events (Observed with lumi. limit corrected for acceptance)
	const std::pair<TH1*, TH1*>& getHDistributionsCorrected() const;

	//get histograms for resulting sum
	TH1* getHSum() const;

	//get histograms for resulting difference
	TH1* getHDifference() const;

	//get histograms for resulting asymmetry
	TH1* getHAsymmetry() const;

private:

	//get appropriate histogram
	TH1* getH(const std::pair<TH1*, TH1*>& histogramPair, int beamPolarisation) const;

	std::pair<TH1*, TH1*> m_hDistributions;		//histograms to store events
	std::pair<TH1*, TH1*> m_hDistributionsObserved;	//histograms to store events (True kinematics)
	std::pair<TH1*, TH1*> m_hDistributionsTrue;	//histograms to store events (True kinematics)
	std::pair<TH1*, TH1*> m_hDistributionsBorn;	//histograms to store events (Born kinematics)
	std::pair<TH1*, TH1*> m_hDistributionsAcceptance;	//histograms to Born/True ratio
	std::pair<TH1*, TH1*> m_hDistributionsRC;	//histograms to Born/True ratio
	std::pair<TH1*, TH1*> m_hDistributionsCorrected;	//histograms to store events

	TH1* m_hAsymmetry;							//histogram for resulting asymmetry
	TH1* m_hSum;								//histogram for sum
	TH1* m_hDif;								//histogram for difference
	TF1* m_fFit;								//function for fitting 
};

#endif