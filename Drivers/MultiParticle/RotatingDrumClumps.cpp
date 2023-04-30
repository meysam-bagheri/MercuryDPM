//Copyright (c) 2013-2022, The MercuryDPM Developers Team. All rights reserved.
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

// Example 4 - Gomboc (rolly - polly out of simply-commected shape of constant density)

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/MultiParticle.h"
#include "clump/ClumpIO.h"
#include "clump/Mercury3DClump.h"
# include <stdlib.h>
#include <Boundaries/CubeInsertionBoundary.h>

//Mdouble f_min = -2; Mdouble f_max = 2;
int SAVECOUNT = 400;

class multiParticleT1 : public Mercury3Dclump
{

    // Group id of the rotating geometry
    unsigned rotatingWallID = 0;

    // Pointers to the insertion boundaries
    CubeInsertionBoundary* insertionBoundary;

public:
    explicit  multiParticleT1()
    {
        setGravity(Vec3D(0,-9.8,0));
        // Set name of output files
        setName("RotatingDrumClumps");
        setXBallsAdditionalArguments("-solidf -v0");
        // Set domain size
        setMin({-3.75,-3.75,0});
        setMax({3.75,3.75,4.25});
        // Output files: wall-vtu
        setWallsWriteVTK(true);

        /*
        setXMax(2*f_max);
        setYMax(2*f_max);
        setZMax(f_max);
        setXMin(2*f_min);
        setYMin(2*f_min);
        setZMin(f_min);
         */
        load_clumps(data);
        setClumpIndex(0);
        clump_mass = data.mass[clump_index];
    }

    void setClumpDamping(Mdouble damp){ clump_damping = damp;}

    void setClumpIndex(Mdouble index){ clump_index = index;}

    Mdouble getClumpMass(){return clump_mass;}

    void setupInitialConditions() override
    {
        //Contact law and density
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        species->setDensity(1.0); // sets the species type-0 density
        species->setDissipation(50.0);
        species->setStiffness(1e6);
        species->setSlidingFrictionCoefficient(0.6);
        species->setSlidingStiffness(5e5);
        species->setRollingFrictionCoefficient(0.0);
        species->setRollingStiffness(5e5);
        const Mdouble collisionTime = species->getCollisionTime(getClumpMass());
        setClumpDamping(0);

        // Set time step
        setTimeStep(collisionTime / 50.0);

        // Generate gomboc
        setClumpIndex(2);
        MultiParticle p0;
        p0.setSpecies(speciesHandler.getObject(0)); // Assign the material type to MultiParticle 1
        p0.setMaster();
        //data = rotate_clump(data, clump_index, uniform_random_pds(seed)); // here you can try different seeds
        p0.setRadius(data.pebbles_r[clump_index][0]);
        Vec3D pos = Vec3D(0, 0, 1);
        p0.setPosition(pos);
        for (int j = 0; j < data.pebbles_r[clump_index].size(); j++) {
            p0.addSlave(Vec3D(data.pebbles_x[clump_index][j],
                              data.pebbles_y[clump_index][j],
                              data.pebbles_z[clump_index][j]),
                        data.pebbles_r[clump_index][j]);
        }
        p0.setPrincipalDirections(
                Matrix3D(data.pd[clump_index][0], data.pd[clump_index][1], data.pd[clump_index][2],
                         data.pd[clump_index][3], data.pd[clump_index][4], data.pd[clump_index][5],
                         data.pd[clump_index][6], data.pd[clump_index][7], data.pd[clump_index][8]));
        p0.setInitInertia(
                MatrixSymmetric3D(data.toi[clump_index][0], data.toi[clump_index][1], data.toi[clump_index][2],
                                  data.toi[clump_index][4], data.toi[clump_index][5],
                                  data.toi[clump_index][8]));
        p0.setMassMultiparticle(data.mass[clump_index]);
        double mag = 0;
        p0.setAngularVelocity(Vec3D(0,0,0));
        p0.setVelocity(Vec3D(0,0,0));

        p0.setDamping(clump_damping);
        particleHandler.copyAndAddObject(p0);

        // Add particles
        double insertionVolume = 0.5/species->getDensity();  //mass of particles: 0.5 Kg
        logger(INFO,"Insertion volume: %", insertionVolume);
        CubeInsertionBoundary boundary;

        boundary.set(p0, 1000, Vec3D(-2.5,-2.5,0), Vec3D(2.5,2.5,4.0),Vec3D(0,0,0),Vec3D(0,0,0));
        boundary.setInitialVolume(insertionVolume);
        insertionBoundary = boundaryHandler.copyAndAddObject(boundary);

        // Introduce a rotating wall
        Mdouble wallScaleFactor = 1e-3; // Scale used in the stl file (mm)
        Vec3D shift = {0,0,0};
        Vec3D velocity = {0,0,0};
        rotatingWallID = wallHandler.readTriangleWall("RotatingDrum.stl",speciesHandler.getObject(0), wallScaleFactor,shift,velocity,Vec3D(0,0,0));

        /*
        // Rectangular box
        wallHandler.clear();
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0, 0, getZMax()));
        wallHandler.copyAndAddObject(w0);
         */
    }

    void actionsAfterTimeStep() override {
        if (getNumberOfTimeSteps()/SAVECOUNT == getNumberOfTimeSteps()/(double) SAVECOUNT ){
            for (std::vector<BaseParticle*>::iterator it= particleHandler.begin(); it!=particleHandler.end(); ++it){
                if ((*it)->IsMaster()) {

                    std::cout<<"Saving timestep "<<getNumberOfTimeSteps() <<std::endl;

                    // Add velocity to log
                    std::ofstream funct("clump_seq.txt", std::ios_base::app | std::ios_base::out);
                    funct << getNumberOfTimeSteps() <<" "
                    << static_cast<MultiParticle*>(*it)->getPosition()<<" "
                    << static_cast<MultiParticle*>(*it)->getPrincipalDirections_e1()<<" "
                    << static_cast<MultiParticle*>(*it)->getPrincipalDirections_e2()<<" "
                    << static_cast<MultiParticle*>(*it)->getPrincipalDirections_e3()<<" "
                    <<"\n";
                    funct.close();

                }
            }

        }
        // Start the rotation after the particles have been inserted
        static bool rotationStarted = false;
        Vec3D angularVelocity = Vec3D(0,0,1.0/12.0*constants::pi);
        if (!rotationStarted and insertionBoundary->getInitialVolume()<=insertionBoundary->getVolumeOfParticlesInserted()) {
            rotationStarted = true;
            logger(INFO,"Starting rotation");
            for (const auto wall : wallHandler) {
                if (wall->getGroupId()==rotatingWallID) {
                    wall->setAngularVelocity(angularVelocity);
                }
            }
        }
    }


private:
    int clump_index;
    clump_data data;
    Mdouble clump_mass;
    Mdouble clump_damping = 1;
};

int main(int argc, char* argv[])
{
    multiParticleT1 problem;
     // Set final time and how often to output
    problem.setSaveCount(SAVECOUNT);
    problem.setTimeMax(50.0);
    // Remove existing vtk files
    problem.removeOldFiles();
    // Start solving in time
    problem.solve();
    return 0;
}
