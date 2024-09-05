//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include "Species/Species.h"
#include "Species/HertzianViscoelasticMindlinLiquidBridgeClassicalWilletSpecies.h"

/// In this file two particles are symmetrically placed in a bi-axial box are allowed to jump around under gravity. It tests walls gravity and symmetry.

class TwoParticleClassicalWilletCollision : public Mercury2D
{


    void setupInitialConditions() override {
        setMax({0.01,0.01,0.0});
        setGravity({0.0,-9.8,0.0});
        setParticlesWriteVTK(writeVTK);
        wallHandler.setWriteVTK(writeVTK);
        
        SphericalParticle p0, p1;
        p0.setSpecies(speciesHandler.getObject(0));
        p1.setSpecies(speciesHandler.getObject(0));

        p0.setPosition(Vec3D(0.006, 0.005, 0.0));
        p1.setPosition(Vec3D(0.004, 0.005, 0.0));

        p0.setVelocity(Vec3D(-0.1, 0.0, 0.0));
        p1.setVelocity(Vec3D(0.1, 0.0, 0.0));

        p0.setRadius(0.0005);
        p1.setRadius(0.0012);
        particleHandler.copyAndAddObject(p0);
        particleHandler.copyAndAddObject(p1);

        wallHandler.clear();
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1, 0, 0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1, 0, 0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, -1, 0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, 1, 0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);
    }
    bool writeVTK = true;

};

int main(int argc UNUSED, char* argv[] UNUSED)
{

    TwoParticleClassicalWilletCollision twoParticleClassicalWilletCollisionProblem;
    twoParticleClassicalWilletCollisionProblem.setName("TwoParticleClassicalWilletCollisionSelfTest");

    HertzianViscoelasticMindlinLiquidBridgeClassicalWilletSpecies species;
    species.setDensity(2000);
    species.setDissipation(0.16);
    species.setSlidingDissipation(0.16);
    Mdouble particleElasticModulus = 5e6; // taking 100 times that of soft material
    Mdouble particlePoissonRatio = 0.35;
    Mdouble particleSlidingFrictionCoefficient = 0.50;
    //Mdouble particleRadius = 0.0011;
    Mdouble particleShearModulus = particleElasticModulus/ 2 /(1+particlePoissonRatio);
    Mdouble particleEffectiveShearModulus = particleShearModulus/(2-particlePoissonRatio);
    Mdouble particleEffectiveElasticModulus = particleElasticModulus/(1-particlePoissonRatio*particlePoissonRatio);
    species.setEffectiveElasticModulus(particleEffectiveElasticModulus/2);
    species.setEffectiveShearModulus(particleEffectiveShearModulus/2);
    species.setSlidingFrictionCoefficient(particleSlidingFrictionCoefficient);    
    species.setLiquidBridgeVolume(8e-11);
    species.setContactAngle(0);
    species.setSurfaceTension(0.070);

    Mdouble radius = 0.0005;
    Mdouble rayleighTime = helpers::getRayleighTime(radius,species.getEffectiveShearModulus(), particlePoissonRatio,species.getDensity());

    twoParticleClassicalWilletCollisionProblem.setTimeStep(rayleighTime*0.1);
    twoParticleClassicalWilletCollisionProblem.speciesHandler.copyAndAddObject(species);

    twoParticleClassicalWilletCollisionProblem.setTimeMax(0.25);
    twoParticleClassicalWilletCollisionProblem.setSaveCount(twoParticleClassicalWilletCollisionProblem.getTimeMax()/twoParticleClassicalWilletCollisionProblem.getTimeStep());
    twoParticleClassicalWilletCollisionProblem.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    twoParticleClassicalWilletCollisionProblem.solve();
}
