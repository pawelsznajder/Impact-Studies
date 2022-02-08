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
	Analysis(const std::string& className) : BaseObject(className){
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
};

#endif