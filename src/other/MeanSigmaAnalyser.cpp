#include "../../include/other/MeanSigmaAnalyser.h"

#include <iostream>
#include <cmath>

#include "../../include/other/HashManager.h"

MeanSigmaAnalyser::MeanSigmaAnalyser(const std::string& title, size_t nBins, const std::pair<double, double>& range, size_t maxNEvents) 
	: BaseObject("MeanSigmaAnalyser"){

	m_title = title;
	m_nBins = nBins;
	m_range = range;
	m_maxNEvents = maxNEvents;

	reset();
}

MeanSigmaAnalyser::MeanSigmaAnalyser(TH1* h, size_t maxNEvents) : BaseObject("MeanSigmaAnalyser") {

	m_title = h->GetTitle(); 
	m_nBins = h->GetNbinsX();
	m_range = std::make_pair(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
	m_maxNEvents = maxNEvents;

	reset();
}

MeanSigmaAnalyser::~MeanSigmaAnalyser(){
}

void MeanSigmaAnalyser::reset(){

	m_events.clear();
	m_nEvents.clear();

	if(m_maxNEvents > 0){

		m_events = std::vector<std::vector<std::pair<double, double> > >
			(m_nBins, std::vector<std::pair<double, double> >(m_maxNEvents));

		m_nEvents = std::vector<size_t>(m_nBins, 0);
	}else{
		m_events = std::vector<std::vector<std::pair<double, double> > >
			(m_nBins, std::vector<std::pair<double, double> >());

		m_nEvents = std::vector<size_t>(m_nBins, 0);
	}
}

void MeanSigmaAnalyser::fill(double x, double value, double weight){	

	int bin = std::floor(m_nBins * (x - m_range.first) / (m_range.second - m_range.first));

	if(bin < 0 || bin >= m_nBins) return;

	if(m_maxNEvents > 0){

		if(m_nEvents.at(bin) >= m_maxNEvents) return;

		m_events.at(bin).at(m_nEvents.at(bin)) = std::make_pair(value, weight);
		m_nEvents.at(bin)++;
	}else{

		m_events.at(bin).push_back(std::make_pair(value, weight));
		m_nEvents.at(bin)++;
	}
}


TH1* MeanSigmaAnalyser::getH(){

	TH1* h = new TH1D((HashManager::getInstance()->getHash()).c_str(), m_title.c_str(), m_nBins, m_range.first, m_range.second);

	for(size_t i = 0; i < m_nBins; i++){

		double mean = getMean(m_nEvents.at(i), m_events.at(i));
		double rmse = getRMSE(mean, m_nEvents.at(i), m_events.at(i));

		std::cout << mean << " " << rmse << std::endl;

		h->SetBinContent(i + 1, mean);
		h->SetBinError(i + 1, rmse);
	}

	return h;
}

double MeanSigmaAnalyser::getMean(size_t nEvents, const std::vector<std::pair<double, double> >& v) const{

	double result = 0.;
	double sum = 0.;

	for(size_t i = 0; i < nEvents; i++){

		result += v.at(i).first * v.at(i).second;
		sum += v.at(i).second;
	}

	return (result == 0.)?(0.):(result / sum);
}

double MeanSigmaAnalyser::getRMSE(double mean, size_t nEvents, const std::vector<std::pair<double, double> >& v) const{

	double result = 0.;
	double sum = 0.;

	for(size_t i = 0; i < nEvents; i++){

		result += pow(v.at(i).first - mean, 2) * v.at(i).second;
		sum += v.at(i).second;
	}

	return (result == 0.)?(0.):(result / sum);
}
