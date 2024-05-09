//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/*
** This file is annotated with DoxyFile comments in order to show the code on
** the documentation - This is not needed for your real drivers.
** Please ignore these comments.
**
** For full documentation of this code, go to http://docs.mercurydpm.org/Alpha/d0/db0/BeginnerTutorials.html#T1
*/

//! [T1:headers]
#include <Mercury3D.h>
#include <Species/Species.h>
#include <Species/LinearViscoelasticSpecies.h>
//! [T1:headers]

//! Tutorial1: Simulates a single particle moving at constant velocity
//! [T1:class]
class Tutorial1 : public Mercury3D 
{

public:

    //! Use setupInitialConditions to define your particle and wall positions
    void setupInitialConditions() override {
        //! [T1:createParticle]
        // create a new particle
        SphericalParticle p0;
        // set species (material type)
        p0.setSpecies(speciesHandler.getObject(0));
        // set particle radius, position, velocity
        p0.setRadius(0.05);
        p0.setPosition(Vec3D(0.1 * getXMax(), 0.1 * getYMax(), 0.1 * getZMax()));
        p0.setVelocity(Vec3D(0.5, 0.1, 0.1));
        // pass particle to the particle handler (which contains all particles in the simulation)
        particleHandler.copyAndAddObject(p0);
        //! [T1:createParticle]
    }
};
//! [T1:class]

//! [T1:main]
int main(int argc, char* argv[])
{
    // Instantiate an object of class Tutorial1
    Tutorial1 problem;

    // Problem setup: set file name, gravity, domain size, simulation time
    //! [T1:problemSetup]
    problem.setName("Tutorial1");
    problem.setGravity(Vec3D(0.0, 0.0, 0.0));
    problem.setXMax(1.0);
    problem.setYMax(1.0);
    problem.setZMax(1.0);
    problem.setTimeMax(2.3);
    //! [T1:problemSetup]

    // Set the species (material properties such as density and stiffness) of particles and walls
    // (in this case, both are assumed to consist of the same material)
    //! [T1:speciesProp]
    LinearViscoelasticSpecies species;
    species.setDensity(2500.0); //sets the species type_0 density in kg/m3
    // The normal spring stiffness and normal dissipation is computed such that
    // collision time tc=0.005 and restitution coefficient rc=1.0
    species.setStiffness(258.5);//sets the spring stiffness in kN/m.
    species.setDissipation(0.0); //sets the dissipation.
    problem.speciesHandler.copyAndAddObject(species);
    //! [T1:speciesProp]

    // Define what output gets written
    //! [T1:output]
    // number of time steps skipped between saves (i.e. every 10-th time step is written to file)
    problem.setSaveCount(10);
    // creates file with particle positions, velocities, radii for multiple time steps (for plotting)
    problem.dataFile.setFileType(FileType::ONE_FILE);
    // file with contact forces, overlaps for multiple time steps (for plotting)
    problem.fStatFile.setFileType(FileType::NO_FILE);
    // file with all parameters of the last saved time step (for restarting)
    problem.restartFile.setFileType(FileType::ONE_FILE);
    // writes global properties like kinetic/potential energy and center of mass (for quick analysis)
    problem.eneFile.setFileType(FileType::NO_FILE);
    //! [T1:output]

    //! [T1:xballsOutput]
    // additional arguments passed into the .xballs file
    problem.setXBallsAdditionalArguments("-solidf -v0");
    //! [T1:xballsOutput]

    //! [T1:paraviewOutput]
    // whether the wall geometry is written into a vtu file (either once initially or for several time steps)
    problem.wallHandler.setWriteVTK(FileType::ONE_FILE);
    // whether the particle positions are written into a vtu file
    problem.setParticlesWriteVTK(true);
    //! [T1:paraviewOutput]

    //! [T1:solve]
    // set a time step 1/50th of collision time
    problem.setTimeStep(0.005 / 50.0);
    // start the solver (calls setupInitialConditions and initiates time loop)
    problem.solve(argc, argv);
    //! [T1:solve]

    return 0;
}
//! [T1:main]
