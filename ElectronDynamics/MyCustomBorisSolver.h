/*
  MyCustomBorisSolver.h

  Class containing an implementation of a Boris push method for solving electron dynamics
*/

#ifndef MY_CUSTOM_BORIS_SOLVER_H
#define MY_CUSTOM_BORIS_SOLVER_H

#include "ElectronDynamics/BaseField.h"
#include "BasicFunctions/Constants.h"
#include "BasicFunctions/BasicFunctions.h"

#include "ElectronDynamics/BorisSolver.h"

#include "TMath.h"
#include "TVector3.h"

#include <tuple>


namespace rad 
{
  class CustomBorisSolver : public BorisSolver {
    public:
      double mass;
      double charge;
      double tau;
      BaseField* field;

     
      CustomBorisSolver(BaseField* field_v, const double charge_v, const double mass_v, const double tau_v) ;

      /// Alternative get_omega if the field is specified for speed up
      TVector3 get_omega(const TVector3 pos, TVector3 BField);


      /// Calculate radiation acceleration from position, velocity and B field. 
      /// \param pos Position of the charge
      /// \param vec Velocity of the charge
      /// \param BField B field 
      /// \Returns a 3-vector of the acceleration from the RR force
      TVector3 custom_radiation_acceleration(const TVector3 pos, const TVector3 vel, TVector3 BField);


      /// Advances position and velocity vector by a set time, uses the version of get_omega that uses the B field as an input to avoid extra field evaluations.
      /// \param time_step The time step over which to advance the charge's motion
      /// \param x0 Vector of the charge's starting position
      /// \param v0 Vector of the charge's starting velocity
      /// \returns Tuple containing (1) the output position vector and (2) the output velocity vector
      std::tuple<TVector3, TVector3> custom_advance_step(const double time_step, const TVector3 x0, const TVector3 v0);



  };


}

#endif
