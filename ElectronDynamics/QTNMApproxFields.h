#ifndef QTNM_APPROX_FIELDS_H
#define QTNM_APPROX_FIELDS_H

#include "TVector3.h"
#include "ElectronDynamics/BaseField.h"
#include "BasicFunctions/Constants.h"
#include "BFieldModels.h"

namespace rad
{
    class McDonaldHarmonic : public BaseField {
        public:
            // MyMcDonaldCoil coil;
            McD_Loop coil;
            TVector3 btBkg;                                         
            McDonaldHarmonic(const double radius, const double current, const double background);                                                  
            TVector3 evaluate_field_at_point(const TVector3 vec);
    };
}
#endif