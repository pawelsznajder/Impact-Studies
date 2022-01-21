#ifndef BINALU_H
#define BINALU_H

#include <TH1.h>
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
	virtual void fill(const DVCSEvent& event, double weight);

	//analyse (make fit)
	//return histogram and result and uncertainty pair
	virtual std::pair<double, double> analyse();

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
	const std::pair<TH1*, TH1*>& getHDistributions() const;

	//get histograms for resulting asymmetry
	TH1* getHAsymmetry() const;

private:

	std::pair<double, double> m_rangeXB;	//xB range
	std::pair<double, double> m_rangeQ2;	//Q2 range
	std::pair<double, double> m_rangeT;		//t range
	std::pair<double, double> m_rangePhi;	//phi range

	double m_sumXB;	//sum of xB values.
	double m_sumQ2;	//sum of Q2 values.
	double m_sumT;	//sum of t values.

	std::pair<TH1*, TH1*> m_hDistributions;	//histograms to store events
	TH1* m_hAsymmetry;						//histogram for resulting asymmetry
};

#endif