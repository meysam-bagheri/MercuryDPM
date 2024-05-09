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

// Tutorial 9

/*
** This file is annotated with DoxyFile comments in order to show the code on
** the documentation - This is not needed for your real drivers.
** Please ignore these comments.
**
** For full documentation of this code, go to http://docs.mercurydpm.org/Alpha/d0/db0/BeginnerTutorials.html#T9
*/

//! [T9:headers]
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Walls/InfiniteWall.h>
//! [T9:headers]

//! [T9:class]
class Tutorial9 : public Mercury3D
{
public:
    
    void setupInitialConditions() override {
        SphericalParticle p0;

        //sets the particle to species type-1
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(0.005);
        p0.setPosition(Vec3D(0.05 * getXMax(), 0.25 * getYMax(), getZMin() + p0.getRadius()));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);
        
        //sets the particle to species type-2
        p0.setSpecies(speciesHandler.getObject(1));
        p0.setRadius(0.005);
        p0.setPosition(Vec3D(0.05 * getXMax(), 0.5 * getYMax(), getZMin() + p0.getRadius()));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);
        
        //sets the particle to species type-3
        p0.setSpecies(speciesHandler.getObject(2));
        p0.setRadius(0.005);
        p0.setPosition(Vec3D(0.05 * getXMax(), 0.75 * getYMax(), getZMin() + p0.getRadius()));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);

        //! [T9:infiniteWalls]
        InfiniteWall w0;
        
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        //! [T9:infiniteWalls]
    }
};
//! [T9:class]

//! [T9:main]
int main(int argc, char* argv[])
{

    // Instantiate an object of class Tutorial 9
    Tutorial9 problem;

    // Problem setup
    problem.setName("Tutorial9");
    problem.removeOldFiles();
    double angle = constants::pi / 180.0 * 20.0;
    problem.setGravity(Vec3D(sin(angle), 0.0, -cos(angle)) * 9.81);
    problem.setXMax(0.3);
    problem.setYMax(0.3);
    problem.setZMax(0.05);
    problem.setTimeMax(0.4);

    // Set the species (material properties such as density and stiffness) of particles and walls
    // (in this case, particle-0 and the wall are made from species-0, particle-1 belongs to species-1, particle-2 belongs to species-2)

    //! [T9:speciesProp]
    // Species-0 properties
    LinearViscoelasticFrictionSpecies species0;
    // The normal spring stiffness and normal dissipation is computed and set as
    // For collision time tc=0.005 and restitution coefficient rc=0.88,
    species0.setDensity(2500.0);//sets the species type_0 density
    species0.setStiffness(259.018);//sets the spring stiffness.
    species0.setDissipation(0.0334);//sets the dissipation.
    species0.setSlidingStiffness(2.0 / 7.0 * species0.getStiffness());
    species0.setRollingStiffness(2.0 / 5.0 * species0.getStiffness());
    species0.setSlidingFrictionCoefficient(0.0);
    species0.setRollingFrictionCoefficient(0.0);
    auto ptrToSp0=problem.speciesHandler.copyAndAddObject(species0);

    // Species-1 properties
    LinearViscoelasticFrictionSpecies species1;
    species1.setDensity(2500.0);//sets the species type-1 density
    species1.setStiffness(259.018);//sets the spring stiffness
    species1.setDissipation(0.0334);//sets the dissipation
    species1.setSlidingStiffness(2.0 / 7.0 * species1.getStiffness());
    species1.setRollingStiffness(2.0 / 5.0 * species1.getStiffness());
    species1.setSlidingFrictionCoefficient(0.5);
    species1.setRollingFrictionCoefficient(0.0);
    auto ptrToSp1=problem.speciesHandler.copyAndAddObject(species1);

    // Properties of contacts between species-0 and species-1
    // (no new species is defined since the mixed species is automatically created)
    auto species01 = problem.speciesHandler.getMixedObject(ptrToSp0,ptrToSp1);
    species01->setStiffness(259.018);//sets the spring stiffness
    species01->setDissipation(0.0334);//sets the dissipation
    species01->setSlidingStiffness(2.0 / 7.0 * species01->getStiffness());
    species01->setRollingStiffness(2.0 / 5.0 * species01->getStiffness());
    species01->setSlidingFrictionCoefficient(0.5);
    species01->setRollingFrictionCoefficient(0.0);

    // Species 2 properties
    LinearViscoelasticFrictionSpecies species2;
    species2.setDensity(2500.0);//sets the species type-2 density
    species2.setStiffness(258.5);//sets the spring stiffness
    species2.setDissipation(0.0);//sets the dissipation
    species2.setSlidingStiffness(2.0 / 7.0 * species2.getStiffness());
    species2.setRollingStiffness(2.0 / 5.0 * species2.getStiffness());
    species2.setSlidingFrictionCoefficient(0.5);
    species2.setRollingFrictionCoefficient(0.5);
    auto ptrToSp2 = problem.speciesHandler.copyAndAddObject(species2);

    // Properties of contacts between species-0 and species-2
    auto species02 = problem.speciesHandler.getMixedObject(ptrToSp0, ptrToSp2);
    species02->setStiffness(259.018);//sets the stiffness
    species02->setDissipation(0.0334);//sets the dissipation
    species02->setSlidingStiffness(2.0 / 7.0 * species02->getStiffness());
    species02->setRollingStiffness(2.0 / 5.0 * species02->getStiffness());
    species02->setSlidingFrictionCoefficient(0.5);
    species02->setRollingFrictionCoefficient(0.5);

    // Note: Properties of contacts between species-0 and species-1 are never defined, since no contacts between those
    // two materials occurs in this simulation; note: if no properties are defined, then default properties are assumed,
    // mostly harmonic means of the properties of both individual species
    //! [T9:speciesProp]

    //Output
    problem.setSaveCount(25);
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::ONE_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);

    //! [T9:paraviewOutput]
    // whether the wall geometry is written into a vtu file (either once initially or for several time steps)
    problem.wallHandler.setWriteVTK(FileType::ONE_FILE);
    // whether the particle positions are written into a vtu file
    problem.setParticlesWriteVTK(true);
    //! [T9:paraviewOutput]

    //![T9: solve]
    // set a time step 1/50th of collision time
    problem.setTimeStep(0.005 / 50.0);
    // start the solver (calls setupInitialConditions and initiates time loop)
    problem.solve(argc, argv);
    //![T9: solve]
    return 0;
}
//! [T9:main]
