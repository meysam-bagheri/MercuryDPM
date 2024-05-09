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

// Example 4 - Gomboc (rolly - polly out of simply-commected shape of constant density)

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/ClumpParticle.h"
#include "../ClumpHeaders/ClumpInput.h"
#include "../ClumpHeaders/Mercury3DClump.h"
# include <stdlib.h>
#include <Boundaries/CubeInsertionBoundary.h>

Mdouble f_min = -4; Mdouble f_max = 4; Mdouble z_shift = 2;
Mdouble margin = 2;
Mdouble z_margin = 2.8;

Mdouble av_min = 0; // range of angular velocities
Mdouble av_max = 0;
Mdouble tv_min = 0; // range of translational velocities
Mdouble tv_max = 0;
int N_att = 1000;   // Number of attempts to add particle


int SAVECOUNT = 400;

class ChangingTOIParticle : public Mercury3Dclump
{

    // Group id of the rotating geometry
    unsigned rotatingWallID = 0;

    // Pointers to the insertion boundaries
    CubeInsertionBoundary* insertionBoundary;

public:
    explicit  ChangingTOIParticle()
    {
        setGravity(Vec3D(0,-9.8,0));
        // Set name of output files
        setName("RotatingDrumClumps");
        setXBallsAdditionalArguments("-solidf -v0");
        // Set domain size
        setMin({f_min,f_min,f_min});
        setMax({f_max,f_max,f_max});
        // Output files: wall-vtu

        LoadClumps(data);
        setClumpIndex(2);
        clump_mass = data.mass[clump_index];
    }

    void setClumpDamping(Mdouble damp){ clump_damping = damp;}

    void setClumpIndex(Mdouble index){ clump_index = index;}

    Mdouble getClumpMass(){return clump_mass;}

    void setupInitialConditions() override
    {

        setClumpIndex(2);
        int N_created = 0;
        for (int part = 0; part<N_att; part++) {

            ClumpParticle p0;
            p0.setSpecies(speciesHandler.getObject(0)); // Assign the material type to Clump 1
            p0.setClump();
            ClumpData rdata = RotateClump(data, clump_index, UniformRandomPDs()); // Rotate clump arbitrarily



            p0.setRadius(rdata.pebblesR[clump_index][0]);

            for (int j = 0; j < rdata.pebblesR[clump_index].size(); j++) {
                p0.addPebble(Vec3D(rdata.pebblesX[clump_index][j],
                                   rdata.pebblesY[clump_index][j],
                                   rdata.pebblesZ[clump_index][j]),
                             rdata.pebblesR[clump_index][j]);
            }
            p0.setPrincipalDirections(
                    Matrix3D(rdata.pd[clump_index][0], rdata.pd[clump_index][1], rdata.pd[clump_index][2],
                             rdata.pd[clump_index][3], rdata.pd[clump_index][4], rdata.pd[clump_index][5],
                             rdata.pd[clump_index][6], rdata.pd[clump_index][7], rdata.pd[clump_index][8]));
            p0.setInitInertia(
                    MatrixSymmetric3D(rdata.toi[clump_index][0], rdata.toi[clump_index][1], rdata.toi[clump_index][2],
                                      rdata.toi[clump_index][4], rdata.toi[clump_index][5],
                                      rdata.toi[clump_index][8]));
            p0.setClumpMass(rdata.mass[clump_index]);

            p0.setDamping(clump_damping);


            Vec3D pos = Vec3D(f_min + margin + RandomDouble(f_max - f_min - 2 * margin),
                              f_min + margin + RandomDouble(f_max - f_min - 2 * margin),
                              z_shift + f_min + z_margin + RandomDouble(f_max - f_min - 2 * z_margin));

            p0.setPosition(pos);

            Vec3D angVel = Vec3D(av_min + RandomDouble(av_max - av_min),
                                 av_min + RandomDouble(av_max - av_min),
                                 av_min + RandomDouble(av_max - av_min));

            Vec3D vel = Vec3D(tv_min + RandomDouble(tv_max - tv_min),
                              tv_min + RandomDouble(tv_max - tv_min),
                              tv_min + RandomDouble(tv_max - tv_min));

            p0.setAngularVelocity(angVel);
            p0.setVelocity(vel);


            if (checkClumpForInteractionPeriodic(p0)) {
                particleHandler.copyAndAddObject(p0);
                N_created++;
            }
        }
        std::cout<<"N_created = "<<N_created<<std::endl;

        // Introduce a rotating wall
        Mdouble wallScaleFactor = 1e-3; // Scale used in the stl file (mm)
        Vec3D shift = {0,0,0};
        Vec3D velocity = {0,0,0};
        rotatingWallID = wallHandler.readTriangleWall(getMercuryDPMSourceDir() + "/Drivers/Clump/RotatingDrum/RotatingDrum.stl",speciesHandler.getObject(0), wallScaleFactor,shift,velocity,Vec3D(0,0,0));
        wallHandler.setWriteVTK(true);

    }

    void actionsAfterTimeStep() override {
        Vec3D angularVelocity = Vec3D(0,0,1.0/12.0*constants::pi);
        for (const auto wall : wallHandler) {
                if (wall->getGroupId()==rotatingWallID) {
                    wall->setAngularVelocity(angularVelocity);
                }
            }
    }


private:
    int clump_index;
    ClumpData data;
    Mdouble clump_mass;
    Mdouble clump_damping = 1;
};

int main(int argc, char* argv[])
{
    ChangingTOIParticle problem;
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    species->setDensity(1.0); // sets the species type-0 density
    species->setDissipation(50.0);
    species->setStiffness(1e6);

    species->setSlidingFrictionCoefficient(0.6);
    species->setSlidingStiffness(5e5);
    species->setRollingFrictionCoefficient(0.0);
    species->setRollingStiffness(5e5);

    const Mdouble collisionTime = species->getCollisionTime(problem.getClumpMass());
    problem.setClumpDamping(0);
    problem.setTimeStep(collisionTime / 50.0);
    problem.setSaveCount(SAVECOUNT);
    problem.setHGridMaxLevels(1);
    //problem.setTimeMax(48.0);
    problem.setTimeMax(0.08);

    problem.removeOldFiles();
    problem.solve();
    return 0;
}
