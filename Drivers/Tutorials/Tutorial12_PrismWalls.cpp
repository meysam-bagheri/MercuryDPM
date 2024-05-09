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

// Tutorial 12

/*
** This file is annotated with DoxyFile comments in order to show the code on
** the documentation - This is not needed for your real drivers.
** Please ignore these comments.
**
** For full documentation of this code, go to http://docs.mercurydpm.org/Alpha/d0/db0/BeginnerTutorials.html#T8
*/

//! [T12:headers]
#include <Mercury3D.h>
#include <Walls/InfiniteWall.h>
#include <Walls/IntersectionOfWalls.h>
#include <Species/LinearViscoelasticSpecies.h>
//! [T12:headers]

//! [T12:class]
class Tutorial12 : public Mercury3D
{
public:
    //! [T12:initialConditions]
    void setupInitialConditions() override {
		//! [T12: Define particle]
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(0.005);
        p0.setPosition(Vec3D(0.15 * getXMax(), 0.3335 * getYMax(), 0.0));
        p0.setVelocity(Vec3D(1.2, 0.2, 0.0));
        particleHandler.copyAndAddObject(p0);
        //! [T12: Define particle]
		
		//! [T12: Defining outer walls]
		// Creating outer walls
        InfiniteWall w0;
        
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        
        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        
        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0.0, getYMax(), 0.0));
        wallHandler.copyAndAddObject(w0);
        
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0, getYMin(), 0.0));
        wallHandler.copyAndAddObject(w0);
		//! [T12: Defining outer walls]

		//! [T12: Create object]
		// Defining the object in het center of the box
		// Create an intersection of walls object w1
        IntersectionOfWalls w1;
        // Set the species of the wall
        w1.setSpecies(speciesHandler.getObject(0));
        std::vector<Vec3D> Points(4);
        // Define the vertices of the insert and ensure that points are in clockwise order
        Points[0] = Vec3D(0.25 * getXMax(), -0.25 * getYMax(), 0.0);
        Points[1] = Vec3D(-0.25 * getXMax(), -0.25 * getYMax(), 0.0);
        Points[2] = Vec3D(-0.25 * getXMax(), 0.25 * getYMax(), 0.0);
        Points[3] = Vec3D(0.25 * getXMax(), 0.25 * getYMax(), 0.0);
        // Creating the object from the vector containing points
        w1.createPrism(Points);
        //! [T12: Create object]
        /* Define the angular velocity of the object (optional).
         * Setting the angular velocity of the object results in rotation of the object around
         * the origin (0.0,0.0,0.0). If the object is build such that the center of rotation of the object
         * coincides with the origin the object will rotate around its own axis. */
        //! [T12: Set Angular Velocity]
        w1.setAngularVelocity(Vec3D(0.0,0.0,1.0));
        //! [T12: Set Angular Velocity]
        /* Define the translational velocity of the object (optional) */
        //! [T12: Set Velocity]
        w1.setVelocity(Vec3D(0.0,0.0,0.0));
        //! [T12: Set Velocity]
        /* Set the desired position of center of the object (optional).
         * With this command you can position your object (with the set angular velocity)
         * at a specific location in the domain.*/
        //! [T12: Set Position]
        w1.setPosition(Vec3D(0.5 * getXMax(),0.5 * getYMax(),0.0));
        //! [T12: Set Position]
        // Add the object to wall w1
        wallHandler.copyAndAddObject(w1);
        
    }
    //! [T12:initialConditions]
};
//! [T12:class]

//! [T12:main]
int main(int argc, char* argv[])
{
    //! [T12: constructor]
    // Problem setup
    Tutorial12 problem; // instantiate an object of class Tutorial 12
    
    problem.setName("Tutorial12");
    problem.setSystemDimensions(2);
    problem.setGravity(Vec3D(0.0, 0.0, 0.0));
    problem.setXMax(0.5);
    problem.setYMax(0.5);
    problem.setTimeMax(5.0);

    //! [T12:speciesProp]
    // The normal spring stiffness and normal dissipation is computed and set as
    // For collision time tc=0.005 and restitution coefficient rc=1.0,
    LinearViscoelasticSpecies species;
    species.setDensity(2500.0); //sets the species type_0 density
    species.setStiffness(258.5);//sets the spring stiffness.
    species.setDissipation(0.0); //sets the dissipation.
    problem.speciesHandler.copyAndAddObject(species);

    //! [T12:speciesProp]
    
    problem.setSaveCount(10);
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);
    
    problem.setXBallsAdditionalArguments("-solidf -v0 -s .85");
    
    problem.setTimeStep(.005 / 50.0); // (collision time)/50.0
    //! [T12: constructor]
    
    //! [T12: Solve Problem]
    problem.solve(argc, argv);
    //! [T12: Solve Problem]
    
    return 0;
}
//! [T12:main]
