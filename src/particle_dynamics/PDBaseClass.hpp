
# ifndef PDBASECLASS_H
# define PDBASECLASS_H

# include "../utils/CommonParams.h"
# include "../utils/GetPot"
using namespace std;

// ---------------------------------------------------------------------
// This is the base class for particle-dynamics classes in the PD App.
// This class serves as an interface, and contains a factory method.
// It also contains basic elements of particle simulation (e.g. velocity-
// verlet, etc.)
// ---------------------------------------------------------------------

class PDBaseClass {

public:

   // -------------------------------------------------------------------
   // Factory method that creates objects of sub-classes:
   // -------------------------------------------------------------------

   static PDBaseClass* PDFactory(const CommonParams&, const GetPot&);

   // -------------------------------------------------------------------
   // pure virtual functions:
   // -------------------------------------------------------------------

   virtual void initParticles() = 0;
   virtual void fijFunc(int,int) = 0;

   // -------------------------------------------------------------------
   // common functions:
   // -------------------------------------------------------------------

   PDBaseClass(const CommonParams&, const GetPot&);
 	 ~PDBaseClass();
   void updateParticles();
   void outputParticles();
   void setTimeStep(int step) {current_step = step;}

protected:

   int current_step;

private:

  double dt;
  double dtover2;

};

# endif  // PDBASECLASS_H
