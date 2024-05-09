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
** For full documentation of this code, go to http://docs.mercurydpm.org/Alpha/d0/db0/BeginnerTutorials.html#T10
*/

// Tutorial 10: This tutorial shows how to load restart and data files

//! [T10:headers]
#include <Mercury3D.h>
#include <Math/Helpers.h>
//! [T10:headers]

class DPM : public Mercury3D {
    void setupInitialConditions () override {
        //for MPI cases, particles have to be inserted in setupInitialConditions (more specifically, after decompose is called). You cannot insert particles in main()
        readDataFile("Tutorial10MPI.ini.data");
    }
};


//! [T10:main]
int main(int argc, char* argv[])
{
    //writeToFile is used here to create a restart and a data file, which will be loaded below.
    helpers::writeToFile("Tutorial10MPI.ini.restart",
        "restart_version 1.0 name Tutorial10MPI\n"
        "dataFile fileType ONE_FILE saveCount 10 counter 0 nextSavedTimeStep 0\n"
        "fStatFile fileType NO_FILE saveCount 10 counter 0 nextSavedTimeStep 0\n"
        "eneFile fileType ONE_FILE saveCount 10 counter 0 nextSavedTimeStep 0\n"
        "restartFile fileType ONE_FILE saveCount 10 counter 0 nextSavedTimeStep 0\n"
        "statFile fileType ONE_FILE saveCount 10 counter 0 nextSavedTimeStep 0\n"
        "xMin 0 xMax 2 yMin 0 yMax 2 zMin 0 zMax 2\n"
        "timeStep 1e-03 time 0 ntimeSteps 0 timeMax 10\n"
        "systemDimensions 3 particleDimensions 3 gravity 0 0 -1\n"
        "Species 1\n"
        "LinearViscoelasticSpecies id 0 density 1.9098593 stiffness 2000 dissipation 1\n"
        "Walls 1\n"
        "InfiniteWall id 0 indSpecies 0 position 0 0 0 orientation 0 0 0 1 velocity 0 0 0 angularVelocity 0 0 0 0 force 0 0 0 torque 0 0 0 normal 0 0 -1 factor 1\n"
        "Boundaries 0\n"
        "Particles 0\n"
        "Interactions 0\n"
    );

    helpers::writeToFile("Tutorial10MPI.ini.data",
        "1 0 0 0 0 2 2 2\n"
        "1 1 1.5  0 0 0  0.5  0 0 0  0 0 0  0\n"
        //"1 1 0.5  0 0 0  0.5  0 0 0  0 0 0  0\n"
    );

    DPM Tutorial10;
    Tutorial10.setNumberOfDomains({1,1,2});
    
    //use readRestartFile to load information from a restart file
    Tutorial10.readRestartFile("Tutorial10MPI.ini.restart");
    Tutorial10.setRestarted(false);
    
    //now start the calculations
    Tutorial10.solve();

    return 0;
}
//! [T10:main]
