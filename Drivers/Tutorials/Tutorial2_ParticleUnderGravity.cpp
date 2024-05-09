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

// Tutorial 2

/*
** This file is annotated with DoxyFile comments in order to show the code on
** the documentation - This is not needed for your real drivers.
** Please ignore these comments.
**
** For full documentation of this code, go to http://docs.mercurydpm.org/Alpha/d0/db0/BeginnerTutorials.html#T2
*/

//! [T2:headers]
#include <Species/LinearViscoelasticSpecies.h>
#include <Mercury3D.h>
//! [T2:headers]

//! [T2:class]
class Tutorial2 : public Mercury3D
{
public:
    
    void setupInitialConditions() override {
        //! [T2:createParticle]
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(0.05);
        p0.setPosition(Vec3D(0.5 * getXMax(), 0.5 * getYMax(), getZMax()));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);
        //! [T2:createParticle]
    }
    
};
//! [T2:class]

//! [T2:main]
int main(int argc, char* argv[])
{
    // Problem setup
    Tutorial2 problem;

    //! [T2:problemSetup]
    problem.setName("Tutorial2");
    problem.setGravity(Vec3D(0.0, 0.0, -9.81));
    problem.setXMax(1.0);
    problem.setYMax(1.0);
    problem.setZMax(5.0);
    problem.setTimeMax(1.5);
    //! [T2:problemSetup]

    //! [T2:speciesProp]
    // The normal spring stiffness and normal dissipation is computed and set as
    // For collision time tc=0.005 and restitution coefficient rc=1.0,
    LinearViscoelasticSpecies species;
    species.setDensity(2500.0); //sets the species type_0 density
    species.setStiffness(258.5);//sets the spring stiffness.
    species.setDissipation(0.0); //sets the dissipation.
    problem.speciesHandler.copyAndAddObject(species);
    //! [T2:speciesProp]

    //! [T2:output]
    problem.setSaveCount(10);
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);
    //! [T2:output]

    //! [T1:paraviewOutput]
    // whether the wall geometry is written into a vtu file (either once initially or for several time steps)
    problem.wallHandler.setWriteVTK(FileType::ONE_FILE);
    // whether the particle positions are written into a vtu file
    problem.setParticlesWriteVTK(true);
    //! [T1:paraviewOutput]


    //! [T2:solve]
    problem.setTimeStep(.005 / 50.0);// (collision time)/50.0
    problem.solve(argc, argv);
    //! [T2:solve]
    
    return 0;
}
//! [T2:main]
