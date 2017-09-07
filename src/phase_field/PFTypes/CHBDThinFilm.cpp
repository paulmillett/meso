/*
 * CHBDThinFilm.cpp
 * Copyright (C) 2017 joseph <joseph@JMC-WORKSTATION>
 */

#include "CHBDThinFilm.hpp"



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

CHBDThinFilm::CHBDThinFilm(const CommonParams& pin,
           const GetPot& input_params) : CHBD(pin,input_params)
{
    thickness = input_params("PFApp/thickness",1);
}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

CHBDThinFilm::~CHBDThinFilm()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void CHBDThinFilm::initPhaseField()
{

    //	---------------------------------------
    // initialize the particle suspension:
    //	---------------------------------------

    particles.initParticles();
    cp = particles.mapToGrid();

    //	---------------------------------------
    // initialize the concentration fields:
    //	---------------------------------------

    srand(time(NULL)*(p.rank+1));   // set the random seed
    for (int i=0; i<nxyz; i++) {
        double cpi = creal(cp.getValue(i));
        double r = (double)rand()/RAND_MAX;
        double val1 = co + 0.1*(r-0.5);
        double val2 = 1.0-val1;
        val1 *= 1.0 - cpi;
        val2 *= 1.0 - cpi;
        c1.setValue(i,val1);
        c2.setValue(i,val2);
    }
    initWalls();

    //	---------------------------------------
    // initialize the fourier wave-vectors:
    //	---------------------------------------

    calculateKfields();

    //	---------------------------------------
    // Sync the processors:
    //	---------------------------------------

    MPI::COMM_WORLD.Barrier();

}



// -------------------------------------------------------------------------
// Step forward in time the phase-field method:
// -------------------------------------------------------------------------

void CHBDThinFilm::updatePhaseField()
{

    //	---------------------------------------
    // Update particles:
    //	---------------------------------------

    particles.setTimeStep(current_step);
    if (current_step%part_step_skip == 0) {
        particles.calcCapillaryForce(cp,c1,c2);
        if (p.rank == 0) particles.updateParticles();
        cp = particles.mapToGrid();
    }

    //	---------------------------------------
    // Update Cahn-Hilliard:
    //	---------------------------------------

    mapWalls();
    Sfield dfdc1 = 12.0*w*(c1*c1*c1 - c1*c1 + c1*c2*c2 + c1*cp);
    Sfield dfdc2 = 12.0*w*(c2*c2*c2 - c2*c2 + c2*c1*c1 + c2*cp);
    c1.fft(p_forward);
    c2.fft(p_forward);
    dfdc1.fft(p_forward);
    dfdc2.fft(p_forward);
    if (eCH > 0.0)
    {
        // update with electric field in z-direction
        c1 = (c1 - p.dt*k2*dfdc1 - p.dt*eCH*(*kz)*(*kz)*c1)\
             /(1.0 + kap*p.dt*k4);
        c2 = (c2 - p.dt*k2*dfdc2 - p.dt*eCH*(*kz)*(*kz)*c2)\
             /(1.0 + kap*p.dt*k4);
    }
    else
    {
        // update with no electric field
        c1 = (c1 - p.dt*k2*dfdc1)/(1.0 + kap*p.dt*k4);
        c2 = (c2 - p.dt*k2*dfdc2)/(1.0 + kap*p.dt*k4);
    }
    c1.ifft(p_backward);
    c2.ifft(p_backward);

    //	---------------------------------------
    // Sync the processors:
    //	---------------------------------------

    MPI::COMM_WORLD.Barrier();

}



// -------------------------------------------------------------------------
// make top and bottom 2 layers (in z-dir) have values of 1 in the cp
// Sfield in order to simulate walls. Also zero out the concentration
// fields in these regions
// -------------------------------------------------------------------------

void CHBDThinFilm::initWalls()
{
    // create wall on top
    for (int i=0; i<p.nx; i++)
        for (int j=0; j<p.ny; j++)
            for (int k=p.nz-thickness; k<p.nz; k++)
            {
                int index = k+p.nz*j+p.nz*p.ny*i;
                cp.setValue(index,1.0);
                c1.setValue(index,0.0);
                c2.setValue(index,0.0);
            }
    // create wall on bottom
    for (int i=0; i<p.nx; i++)
        for (int j=0; j<p.ny; j++)
            for (int k=0; k<thickness; k++)
            {
                int index = k+p.nz*j+p.nz*p.ny*i;
                cp.setValue(index,1.0);
                c1.setValue(index,0.0);
                c2.setValue(index,0.0);
            }
}



// -------------------------------------------------------------------------
// make top and bottom 2 layers (in z-dir) have values of 1 in the cp
// Sfield in order to simulate walls.
// -------------------------------------------------------------------------

void CHBDThinFilm::mapWalls()
{
    // create wall on top
    for (int i=0; i<p.nx; i++)
        for (int j=0; j<p.ny; j++)
            for (int k=p.nz-thickness; k<p.nz; k++)
            {
                int index = k+p.nz*j+p.nz*p.ny*i;
                cp.setValue(index,1.0);
            }
    // create wall on bottom
    for (int i=0; i<p.nx; i++)
        for (int j=0; j<p.ny; j++)
            for (int k=0; k<thickness; k++)
            {
                int index = k+p.nz*j+p.nz*p.ny*i;
                cp.setValue(index,1.0);
            }
}
