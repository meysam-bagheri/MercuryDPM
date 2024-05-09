//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and MarbleRun forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in MarbleRun form must reproduce the above copyright
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

#include<iostream>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Walls/TriangleWall.h>
#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"

//#define DEBUG_OUTPUT

class MarbleRun : public Mercury3D 
{
    // pointer to the particle
    SphericalParticle* particle;
    // pointer to the species storing the particle/contact properties
    LinearViscoelasticFrictionSpecies* species;
    // restitution is stored separately, because it is not a direct contact property
    double restitutionCoefficient;

public:

    // create default particle and species (properties to be filled by the user)
    MarbleRun() {
        species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        particle = particleHandler.copyAndAddObject(SphericalParticle(species));
    }

    void setParticlePosition(Vec3D pos) {
        particle->setPosition(pos);
    }

    void setParticleRadius(double radius) {
        particle->setRadius(radius);
    }

    void setParticleDensity(double density) {
        species->setDensity(density);
    }

    void setSlidingFrictionCoefficient(double coeff) {
        species->setSlidingFrictionCoefficient(coeff);
    }

    void setRollingFrictionCoefficient(double coeff) {
        species->setRollingFrictionCoefficient(coeff);
    }

    void setTorsionFrictionCoefficient(double coeff) {
        species->setTorsionFrictionCoefficient(coeff);
    }

    double getParticleMass() {
        return species->getMassFromRadius(particle->getRadius());
    }
    
    void setRestitutionCoefficient(double coeff) {
        restitutionCoefficient = coeff;
    }
    
    void loadSTLFile(std::string stlFile) {
        // read assuming the file is written in mm units
        wallHandler.readTriangleWall(stlFile,speciesHandler.getLastObject(),1e-3);
    }

    void includeInDomain(Vec3D pos) {
        if (pos.X<getXMin()) {
            setXMin(pos.X);
        } else if (pos.X>getXMax()) {
            setXMax(pos.X);
        }
        if (pos.Y<getYMin()) {
            setYMin(pos.Y);
        } else if (pos.Y>getYMax()) {
            setYMax(pos.Y);
        }
        if (pos.Z<getZMin()) {
            setZMin(pos.Z);
        } else if (pos.Z>getZMax()) {
            setZMax(pos.Z);
        }
    }

    void setupInitialConditions() override
    {
        // set domain such that both particle and walls are included
        setDomain(particle->getPosition(),particle->getPosition());
        for (const auto wall : wallHandler) {
            auto triangle = static_cast<TriangleWall*>(wall);
            includeInDomain(triangle->getVertices()[0]);
            includeInDomain(triangle->getVertices()[1]);
            includeInDomain(triangle->getVertices()[2]);
        }
        logger(INFO,"Simulation domain set to [%,%]x[%,%],[%,%]",
               getXMin(),getXMax(),getYMin(),getYMax(),getZMin(),getZMax());
        // set restitution and contact time
        double gravityTimeScale = sqrt(2.0*particle->getRadius()/getGravity().getLength());
        double collisionTime = gravityTimeScale/20; //ensures small overlaps
        species->setCollisionTimeAndRestitutionCoefficient(collisionTime,restitutionCoefficient, getParticleMass());
        species->setSlidingDissipation(2./7.*species->getDissipation());
        species->setRollingDissipation(2./5.*species->getDissipation());
        species->setTorsionDissipation(2./5.*species->getDissipation());
        species->setSlidingStiffness(2./7.*species->getStiffness());
        species->setRollingStiffness(2./5.*species->getStiffness());
        species->setTorsionStiffness(2./5.*species->getStiffness());
        // set time step
        setTimeStep(collisionTime/25.0);
        logger(INFO,"Simulating for %s with a time step of %s",getTimeMax(),getTimeStep());
        // output
        removeOldFiles();
        setSaveCount(0.05/getTimeStep());
        fStatFile.setFileType(FileType::NO_FILE);
        restartFile.writeFirstAndLastTimeStep();
        setParticlesWriteVTK(true);
        wallHandler.setWriteVTK(FileType::ONE_FILE);
    }
};

int main() {

    // Set up a problem of type MarbleRun
    MarbleRun dpm;
    // Set name of output files
    dpm.setName("MarbleRunConcept");
    // Set name of output files
    dpm.loadSTLFile("MarbleRun_Concept_STL_DwarshuisScholten.STL");
    // Set physical particle properties
    dpm.setParticlePosition(Vec3D(0.03,1,0.4)); // (0.03,1,0.9)
    dpm.setParticleRadius(0.01);
    // Set material particle properties
    dpm.setParticleDensity(1000);
    // Set contact properties
    dpm.setSlidingFrictionCoefficient(0.5);
    dpm.setRollingFrictionCoefficient(1e-2); // (1e-4)
    dpm.setTorsionFrictionCoefficient(0.0);
    dpm.setRestitutionCoefficient(0.5);
    // set gravity direction
    dpm.setGravity(Vec3D(0,-9.8,0));
    // Set simulation time
    dpm.setTimeMax(10.0);
    // start the solver
    dpm.solve();
}
