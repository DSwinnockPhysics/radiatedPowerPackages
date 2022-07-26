#include "QTNMApproxFields.h"
#include <iostream>
#include <iomanip>
#include <cmath>
//#include "McDonaldDerivatives.h"
#include "BFieldModels.h"

inline constexpr double pi{ 3.141592653589793238462643383279502884 };
double vacPermeability{ 4e-7*pi };


// Class methods
rad::McDonaldHarmonic::McDonaldHarmonic(const double radius, const double current, const double background) : coil( McD_Loop(7, radius, current, 0, 0,0 ) ) {
    btBkg = TVector3(0, 0, background);
}

TVector3 rad::McDonaldHarmonic::evaluate_field_at_point(const TVector3 vec)
{
	double carP[3] = {vec[0],vec[1],vec[2]}; 				// point to calculate field at in cartesian coordinates
	double BCarVec[3] = {0.,0.,0.}; 			            // placeholder for the field in cartesian coordinates
    coil.getB(carP,BCarVec);
    TVector3 totalField{btBkg[0]-BCarVec[0], btBkg[1]-BCarVec[1], btBkg[2]-BCarVec[2]};
    return totalField;
}

