#include <cmath>

#include "../include/analysis/AnalysisALU.h"

int main(){

	AnalysisALU analysisALU;
	analysisALU.analyse();
	analysisALU.plot("test.pdf");

	return 0;
}