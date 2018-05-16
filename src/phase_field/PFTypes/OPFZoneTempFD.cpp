
# include "OPFZoneTempFD.hpp"
# include <iostream>
# include <fstream>


// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

OPFZoneTempFD::OPFZoneTempFD(const CommonParams& pin,
							 const GetPot& input_params) : p(pin), c(p)
{

    // ---------------------------------------
    // set needed parameters:
    // ---------------------------------------

    nxyz = p.nx*p.ny*p.nz;
    nx = p.nx;
    ny = p.ny;
    nz = p.nz;
    deli = (nz+2)*(ny+2);
	delj = (nz+2);
	delk = 1;
    co = input_params("PFApp/co",0.5);
    M = input_params("PFApp/M",1.0);
    noiseStr = input_params("PFApp/noiseStr",0.1);
    wzone = input_params("PFApp/wzone",10.0);
	vzone = input_params("PFApp/vzone",1.0);
	S = input_params("PFApp/S",0.05);
	Tmin = input_params("PFApp/Tmin",313);
	Tmax = input_params("PFApp/Tmax",483);
	N = input_params("PFApp/N",288);
	alpha = input_params("PFApp/alpha",3.9);
	beta = input_params("PFApp/beta",0.028);
	templating = input_params("PFApp/templating",0);
	templateSpacing = input_params("PFApp/templateSpacing",10);
	templateSpacingY = input_params("PFApp/templateSpacingY",10);
}

// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

OPFZoneTempFD::~OPFZoneTempFD()
{

}

// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void OPFZoneTempFD::initPhaseField()
{

    //	---------------------------------------
    // initialize the concentration field:
    //	---------------------------------------

    srand(time(NULL)*(p.rank+1));   // set the random seed
    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                double r = (double)rand()/RAND_MAX;
                double val = co + 0.1*(r-0.5);
                c.setValue(ndx,val);
            }
        }
    }

}

// -------------------------------------------------------------------------
// Step forward in time the phase-field method:
// -------------------------------------------------------------------------

void OPFZoneTempFD::updatePhaseField()
{
	
	//allocate storage space
    SfieldFD mu(p);
    SfieldFD mob(p);
	SfieldFD C_5(p);
	SfieldFD C_7(p);
	
    // ---------------------------------------
    // calculate chemical potential & mobility
    // ---------------------------------------

    c.updatePBCFluxY();

    for (int i=1; i<nx+1; i++) {

		//calculate local mobility
		double Mc;
		if (vzone != 0) { //if zone annealing is active
			//calculate time 
			double time = current_step*p.dt;
		
			//calculate offset
			int xf = int(time*vzone - 0.5*wzone - p.xOff);
		
			//calculate local mobility
			Mc = 0.5*(1 - tanh(6.0*double(i-xf)/wzone));
		}//if zone annealing
		else { 
			//assign constant mobility
			Mc = M;
		}//else

		//calculate local Temperature
        double Temp = Tmin + Mc*(Tmax - Tmin);
		
		//calculate local chi-N
		double chiN  = (beta + alpha/Temp)*N;
		double chiN2 = chiN  * chiN;
		double chiN3 = chiN2 * chiN;
		double chiN4 = chiN2 * chiN2;
		
		//calculate local chemical potential coefficients
		double c_2 = S*(-0.00368*chiN2 - 1.964*chiN + 15.99); 
		double c_4 = S*( 0.01882*chiN2 + 3.176*chiN - 26.9); 
		double c_5 = S*(-0.000003452*chiN4 + 0.0004019*chiN3 - 0.01791*chiN2 + 0.3735*chiN - 1.684); 
		double c_7 = S*( 0.000004339*chiN4 -  0.000446*chiN3 + 0.01729*chiN2 - 0.3200*chiN + 11.59); 
		
		for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                
				//get current position
				int ndx = i*deli + j*delj + k*delk;
				
				//get current concentration 
				double cc = c.getValue(ndx);
				
				//calculate chemical potential
                double dFdc = 2*c_2*(cc-0.5) + 4*c_4*(cc-0.5)*(cc-0.5)*(cc-0.5); 
				
				//write mobility values
				mob.setValue(ndx, Mc);
				
				//write potential values
				mu.setValue(ndx,dFdc);
				
				//write coefficients needed outside current scope
				C_5.setValue(ndx,c_5);
				C_7.setValue(ndx,c_7);
            }//k
        }//j
    }//i
	
	//ensure all processors reach this point before proceeding
	MPI::COMM_WORLD.Barrier();
	
	//update boundary conditions for mu and mobility
    mu.updatePBCFluxY();
    mob.updatePBCFluxY();
	
	// ---------------------------------------
    // Apply surface templating:
    // ---------------------------------------

	if (templating != 0)
	{
		//determine z position for (1)2d or (0)3d simulation
		int k = (nz==1);
		
		//apply templating type
		switch(templating)
		{
			case 1: // horizontal striped templating (X-direction, y-spacing)
			{
				for (int j=1; j<ny+1; j+=templateSpacing) {
					for (int i=1; i<nx+1; i++) {
										
						//set current position
						int ndx = i*deli + j*delj + k*delk;
						c.setValue(ndx,1);
					}//j
				}//i
			}//case 1
			break;
			case 2: // vertical striped templating (Y-direction, x-spacing)
			{
				//determine starting point from x-offset
				int startAt;
				if (p.xOff%templateSpacing > 0) {
					startAt = templateSpacing - (p.xOff%templateSpacing) + 1;
				}//if
				else {
					startAt = 1;
				}//else
				for (int i=startAt; i<nx+1; i+=templateSpacing) {
					for (int j=1; j<ny+1; j++) {
					
						//set current position
						int ndx = i*deli + j*delj + k*delk;
						c.setValue(ndx,1);
					}//j
				}//i
			}//case 2
			break;
			case 3: // square dot templating
			{
				//determine starting point from x-offset
				int startAt;
				if (p.xOff%templateSpacing > 0) {
					startAt = templateSpacing - (p.xOff%templateSpacing) + 1;
				}//if
				else {
					startAt = 1;
				}//else
					
				for (int i=startAt; i<nx+1; i+=templateSpacing) {
					for (int j=1; j<ny+1; j+=templateSpacing) {
						//set current position
						int ndx = i*deli + j*delj + k*delk;
						c.setValue(ndx,1);
					}//j
				}//i
			}//case 3
			break;
			case 4: // uneven dot spacing in x and y
			{
				//determine starting point from x-offset
				int startAt;
				if (p.xOff%templateSpacing > 0) {
					startAt = templateSpacing - (p.xOff%templateSpacing) + 1;
				}//if
				else {
					startAt = 1;
				}//else
					
				for (int i=startAt; i<nx+1; i+=templateSpacing) {
					for (int j=1; j<ny+1; j+=templateSpacingY) {
						//set current position
						int ndx = i*deli + j*delj + k*delk;
						c.setValue(ndx,1);
					}//j
				}//i
			}//case 4
		}//switch templating
	}//if templating
	
		
    // ---------------------------------------
    // update CH equation:
    // ---------------------------------------

	//calculate first laplacian of CH right hand side
	SfieldFD inside = mu - C_5*c.Laplacian();
	inside.updatePBCFluxY();
	
	//calculate total right hand side of CH equation
	SfieldFD RHS = inside.Laplacian(mob); 	
	// Apply energy barrier for Di-block system:
	RHS -= C_7*(c-co);
	
	//add thermal noise
	for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
				
				// Add random fluctuations:
				double r = (double)rand()/RAND_MAX;
                double val = noiseStr*(r-0.5);
                RHS.addValue(ndx,val);
            }
        }
    }
	
	RHS.updatePBCFluxY();
	c += p.dt*RHS;
}//updatePhaseField

// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void OPFZoneTempFD::outputPhaseField()
{
    int iskip = p.iskip;
    int jskip = p.jskip;
    int kskip = p.kskip;
    c.writeVTKFile("c",current_step,iskip,jskip,kskip);
	
}
