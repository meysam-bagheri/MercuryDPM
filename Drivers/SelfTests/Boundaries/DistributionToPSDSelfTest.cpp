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

#include <iostream>
#include <CMakeDefinitions.h>
#include "Mercury3D.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"


/**
 * An example code on how to add particles of a given particle size distribution and flow rate.
 */

class DistributionToPSDSelfTest : public Mercury3D
{
public:
    
    /**
     * Use setupInitialConditions to define a.o. an insertion boundary, which will insert particles in a given region, at a given flow rate and of a given particle size distribution
     */
    void setupInitialConditions() override
    {
        setName("DistributionToPSDSelfTest");
        setSystemDimensions(3);
        setGravity(Vec3D(0, 0, 0));
        setTimeStep(1e-4);
        dataFile.setSaveCount(10);
        setTimeMax(2e-2);
        setHGridMaxLevels(2);
        
        setMin(Vec3D(0, 0, 0));
        setMax(Vec3D(.1, .1, .1));
        
        LinearViscoelasticSpecies species;
        species.setDensity(2000);
        species.setStiffness(10000);
        speciesHandler.copyAndAddObject(species);
        
        // define a cube insertion boundary, which inserts particles in a cuboid region
        CubeInsertionBoundary insertionBoundary;
        
        // define the type of particle you want to insert (radius, position and velocity will be set by the insertion boundary)
        SphericalParticle templateParticle(speciesHandler.getObject(0));
        
        // insert particles in the whole domain (between getMin and getMax) with initial velocity 0 and and a radius between 1 and 2 (uniform number distribution)
        Vec3D posMin = getMin();
        Vec3D posMax = getMax();
        Vec3D velMin = {0, 0, 0};
        Vec3D velMax = {0, 0, 0};
        unsigned maxFail = 1; //insert as quick as possible: try every time step, until you maxFail=1 particle fails to be insertable (overlaps with another particle or wall)
        double radMin = 1;
        double radMax = 2;
        insertionBoundary.set(&templateParticle, maxFail, posMin, posMax, velMin, velMax);
    
        //create a predefined discrete uniform distribution between radMin and radMax with a certain resolution
        // (numberOfBins)
        PSD psd;
        // laser diffraction measurements usually give volumetric responses. This PSD was already converted to a
        // CUMULATIVE_NUMBER_DISTRIBUTION.
//        psd.setPSDFromCSV(getMercuryDPMSourceDir() + "/Drivers/SelfTests/Boundaries/InputData/APAPM.csv",
//                          PSD::TYPE::CUMULATIVE_NUMBER_DISTRIBUTION);
        psd.setDistributionNormal(0.001, 0.0001, 50);
        insertionBoundary.setPSD(psd);
    
        //instead of inserting 1 particle per timestep, insert at a given flow rate, such as 0.001 m^3/s
        // insertionBoundary.setVolumeFlowRate(1e-3);
    
        //add the insertion boundary to the handler
        boundaryHandler.copyAndAddObject(insertionBoundary);
    }
    
    void printTime() const override
    {
        logger(INFO, "t=%, tMax=%, N=%", getTime(), getTimeMax(), particleHandler.getSize());
    }
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    DistributionToPSDSelfTest problem;
    problem.solve();
}
