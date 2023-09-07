#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <string.h>

#include "../include/other/BaseObject.h"
#include "../include/event/DVCSEvent.h"

/*
 * Parent class to those used to perform analyses.
 */
class Analysis : public BaseObject{

public:

	//constructor
	Analysis(const std::string& className, double targetLuminosity) : BaseObject(className), m_targetLuminosity(targetLuminosity){
	}

	//destructor
	virtual ~Analysis(){
	}

	//fill with events
	virtual void fill(DVCSEvent& event, double weight) = 0;

	//perform analysis
	virtual void analyse() = 0;

	//plot histograms to file
	virtual void plot(const std::string& path) = 0;

protected:

	double m_targetLuminosity;	///< target luminosity in fb-1
};

#endif