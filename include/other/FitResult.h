#ifndef FIT_RESULT_H
#define FIT_RESULT_H

#include <vector>

#include "../../include/other/BaseObject.h"

/*
 * Container storing result of a single fit.
 */
class FitResult : public BaseObject{

public:

	//constructor
	FitResult();

	//destructor
	virtual ~FitResult();

	//get status code 
	int getStatusCode() const; 		

	//get number of points
	size_t getNPoints() const;

	//get chi2	
	double getChi2() const;		

	//get normalised chi2	
	double getChi2Norm() const;		

	//set status code 
	void setStatusCode(int statusCode); 		

	//set number of points
	void setNPoints(size_t nPoints);

	//set chi2	
	void setChi2(double chi2);		

	//get ndf
	size_t getNDF() const;

	//add parameter
	void addParameter(const std::pair<double, double>& parameter);	

	//get number of parameters
	size_t getNParameters() const;

	//get parameter
	const std::pair<double, double>& getParameter(size_t i) const;

	//print
	void print() const;

private:

	int m_statusCode; 			//status code
	size_t m_nPoints;			//number of points
	double m_chi2;				//chi2
	std::vector<std::pair<double, double> >	m_parameters;	//values and unc. of fitted parameters
};

#endif