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

// Tutorial 3

/*
** This file is annotated with DoxyFile comments in order to show the code on
** the documentation - This is not needed for your real drivers.
** Please ignore these comments.
**
** For full documentation of this code, go to http://docs.mercurydpm.org/Alpha/d0/db0/BeginnerTutorials.html#T3
*/

//! [T3:headers]
#include <Species/LinearViscoelasticSpecies.h>
#include <Mercury3D.h>
#include <Walls/InfiniteWall.h>
//! [T3:headers]

/** Tutorial3 derives from Mercury3D (3D particle simulation with hGrid contact detection).
 *  It simulates the impact of a particle onto a wall under the influence of gravity.
 */
//! [T3:class]
class Tutorial3 : public Mercury3D
{
public:

    //! Use setupInitialConditions to define your particle and wall positions
    void setupInitialConditions() override {
        // add a particle of 1cm diameter at height zMax
        //! [T3:createParticle]
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(0.005);
        p0.setPosition(Vec3D(0.5 * getXMax(), 0.5 * getYMax(), getZMax()));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);
        //! [T3:createParticle]

        // add a bottom wall at zMin
        //! [T3:infiniteWall]
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        //! [T3:infiniteWall]
    }
    
};
//! [T3:class]

//! [T3:main]
int main(int argc, char* argv[])
{
    
    // Problem setup
    Tutorial3 problem;

    // name the output files
    problem.setName("Tutorial3");
    // remove old output files with the same name
    problem.removeOldFiles();
    // set gravivational acceleration
    problem.setGravity(Vec3D(0.0, 0.0, -9.81));
    // set domain size (xMin = yMin = zMin = 0 by default)
    problem.setXMax(0.5);
    problem.setYMax(0.5);
    problem.setZMax(0.5);
    problem.setTimeMax(2.0);

    //! [T3:speciesProp]
    // Sets a linear spring-damper contact law.
    // The normal spring stiffness and normal dissipation is computed and set as
    // For collision time tc=0.005 and restitution coefficient rc=1.0,
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    species->setDensity(2500.0);
    double mass = species->getMassFromRadius(0.005);
    double collisionTime = 0.005;
    double restitution = 1.0;
    species->setCollisionTimeAndRestitutionCoefficient(collisionTime,restitution,mass);
    logger(INFO, "Stiffness %, dissipation %", species->getStiffness(), species->getDissipation());
    //! [T3:speciesProp]

    //! [T3:output]
    problem.setSaveCount(25);
    //! [T3:output]

    // Output vtk files for ParaView visualisation
    //! [T3:visualOutput]
    problem.wallHandler.setWriteVTK(FileType::ONE_FILE);
    problem.setParticlesWriteVTK(true);
    //! [T3:visualOutput]

    //! [T3:solve]
    // Sets time step to 1/50th of the collision time
    problem.setTimeStep(0.005 / 50.0);
    // Start the solver (calls setupInitialConditions and initiates time loop)
    problem.solve(argc, argv);
    //! [T3:solve]
    return 0;
}
//! [T3:main]
