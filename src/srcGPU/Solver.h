#ifndef Solver_H_
#define Solver_H_


#include "commonData.h"
#include <fstream>
#include <string>
#include <cusparse_v2.h>
class Solver {
	public :
    vector<double> solve3Diag (const vector <double> & lDiag, 
	                           const vector <double> & Diag, 
							   const vector <double> & uDiag,
							   const vector <double> & rHS) ;

    vector<double> SOR3DiagPeriodic (const vector <double> & lDiag, 
	                                 const vector <double> & Diag, 
								     const vector <double> & uDiag,
								     const vector <double> & rHS,
									 vector <double> & firstGuess) ; 

}; 

#endif
