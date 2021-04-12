//Copyright (c) 2013-2021, The MercuryDPM Developers Team. All rights reserved.
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
#include "MercuryOS.h"
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Walls/TriangleWall.h>
#include <Boundaries/DeletionBoundary.h>

/**
 * Granular flow from a silo. Four cases are considered, cases 1 and 2 have a small orifice, cases 3 and 4 a large orifice, with either particles of type M1 or M2. 
 * An initial placement of particles is read from the files Silo/PartCoordinates*.txt. This placement is generated by filling the silo with 100 000 particles of diameter 4 mm of materials M1, and in the second case study of
material M2.
 * The initial placement of the walls is read from the file Silo/walls*.txt. One can change to a smooth wall implementation (non-triangulated) by setting useMercuryWalls_ to true.
 * A deletion zone is placed placed below the outlet – so the particles are removed directly. 
 * The calculations are performed for a time interval of 5 seconds. 
 * The simulation time step for case study with M1 is 1.5*10 -6 and for case study with M2 2*10 -6 s.
 */
class Silo : public MercuryOS
{
    /// Whether the small or large orifice is simulated
    bool useSmallOrifice_ = false;
    /// Whether particles of type M1 or M2 are simulated
    bool useM1_ = false;

public:
    
    // sets the variable useMercuryWalls_
    void useSmallOrifice(bool useSmallOrifice)
    {
        useSmallOrifice_ = useSmallOrifice;
    }
    
    // returns the variable useMercuryWalls_
    bool useSmallOrifice() const
    {
        return useSmallOrifice_;
    }
    
    // sets the variable useMercuryWalls_
    void useM1(bool useM1)
    {
        useM1_ = useM1;
    }
    
    // returns the variable useMercuryWalls_
    bool useM1() const
    {
        return useM1_;
    }
    
    // returns the height of the orifice (needed to define smooth walls)
    double getOrificeHeight() const
    {
        return useSmallOrifice() ? -0.0593008 : -0.0554222;
    }
    
    // returns the radius of the orifice (needed to define smooth walls)
    double getOrificeRadius() const
    {
        return useSmallOrifice() ? 0.02 : 0.03;
    }
    
    // used to set the initial conditions of the particles, walls, species, etc
    void setupInitialConditions() override
    {
        {
            // name of the output files
            std::string size = useSmallOrifice() ? "Small" : "Large";
            std::string material = useM1() ? "M1" : "M2";
            setName("Silo" + size + material);
        }
        
        // turn on gravity
        setGravity({0, 0, -9.8});
        
        // set time step and maximum simulation time
        setTimeStep(useM1() ? 1.5e-6 : 2e-6);
        setTimeMax(5.0);
        
        // set output frequency
        setSaveCount(static_cast<unsigned>(0.1 / getTimeStep()));
        
        //remove files from previous run
        removeOldFiles();
        
        // determine which output files to write
        if (writeOutput()) {
            // write more frequently
            setSaveCount(static_cast<unsigned>(0.01 / getTimeStep()));
            // write ene, data, fstat, restart and vtu files
            setParticlesWriteVTK(true);
            setWallsWriteVTK(true);
            fStatFile.writeFirstAndLastTimeStep();
            restartFile.writeFirstAndLastTimeStep();
        } else {
            // only .ene files are written
            setFileType(FileType::NO_FILE);
            eneFile.setFileType(FileType::ONE_FILE);
        }
    
        //use a shortened simulation in test mode
        if (test()) {
            setSaveCount(20);
            setTimeMax(200*getTimeStep());
            logger(INFO,"Test mode, reduced timeMax to %",getTimeMax());
        }
        
        // set domain for visualisation
        setMax(Vec3D(0.1, 0.1, 0.18));
        setMin(Vec3D(-0.1, -0.1, -0.1));
        
        // define the material properties of M1, M2, steel (see MercuryOS.h)
        setMaterialProperties();
        
        // read particle positions from file
        {
            //open file
            std::string size = useSmallOrifice() ? "small" : "large";
            std::string material = useM1() ? "M1" : "M2";
            std::string fileName = "Silo/PartCoordinates " + material + " (" + size + " orifice).txt";
            std::ifstream file(fileName);
            logger.assert_always(file.is_open(), "File % could not be opened", fileName);
            // read particle positions
            SphericalParticle particle;
            particle.setSpecies(useM1() ? m1 : m2);
            particle.setRadius(2e-3);
            Vec3D pos;
            while (file >> pos) {
                particle.setPosition(pos);
                particleHandler.copyAndAddObject(particle);
            }
            logger(INFO, "Read % particles from %", particleHandler.getSize(), fileName);
        }
        
        // add walls
        if (useMercuryWalls()) {
            // use smooth walls
            AxisymmetricIntersectionOfWalls cylinder;
            cylinder.setSpecies(steel);
            cylinder.setAxis(Vec3D(0, 0, 1));
            cylinder.addObject(Vec3D(1, 0, 0), Vec3D(0.1, 0, 0));
            wallHandler.copyAndAddObject(cylinder);
            
            AxisymmetricIntersectionOfWalls cone;
            cone.setSpecies(steel);
            cone.setAxis(Vec3D(0, 0, 1));
            cone.addObject(Vec3D(1, 0, -1), Vec3D(getOrificeRadius(), 0, getOrificeHeight()));
            cone.addObject(Vec3D(0, 0, 1), Vec3D(getOrificeRadius(), 0, getOrificeHeight()));
            wallHandler.copyAndAddObject(cone);
            logger(INFO, "Created % walls", wallHandler.getSize());
        } else {
            // read triangle walls from file
            // open file
            std::string size = useSmallOrifice() ? "small" : "large";
            std::string fileName = "Silo/walls (" + size + " orifice).txt";
            std::ifstream file(fileName);
            logger.assert_always(file.is_open(), "File % could not be opened", fileName);
            // read header line
            std::string header;
            getline(file, header);
            // read triangle vertices
            TriangleWall wall;
            wall.setSpecies(steel);
            Vec3D a, b, c;
            while (file >> a >> b >> c) {
                wall.setVertices(a, b, c);
                wallHandler.copyAndAddObject(wall);
            }
            logger(INFO, "Read % walls from %", wallHandler.getSize(), fileName);
        }
        
        // delete particles ~2cm below silo
        DeletionBoundary boundary;
        boundary.set(Vec3D(0, 0, -1), 0.02 - getOrificeHeight());
        boundaryHandler.copyAndAddObject(boundary);
    }
    
    // Write requested output to the ene file
    void writeEneTimeStep(std::ostream &os) const override
    {
        //compute the number of particles in the silo
        double particlesInSilo = 0;
        for (auto p : particleHandler) {
            if (p->getPosition().Z > getOrificeHeight()) {
                ++particlesInSilo;
            }
        }
        if (eneFile.getCounter() == 1) os << "Time ParticleNumber\n";
        os << getTime() << ' ' << particlesInSilo << std::endl;
    }
    
    // Also write the ene information to the screen
    void printTime() const override {
        writeEneTimeStep(std::cout);
    }
};

int main(int argc, char **argv)
{
    // create an instance of the Silo class
    Silo dpm;
    // command line arguments:
    dpm.setNumberOfOMPThreads(helpers::readFromCommandLine(argc, argv, "-omp", 1));
    // turn on additional output files for viewing/analysing the data
    dpm.writeOutput(helpers::readFromCommandLine(argc, argv, "-writeOutput"));
    // turn on additional output files for viewing/analysing the data
    dpm.test(helpers::readFromCommandLine(argc, argv, "-test"));
    // read from command line whether large or small orifice is used
    dpm.useSmallOrifice(helpers::readFromCommandLine(argc, argv, "-useSmallOrifice"));
    // read from command line whether large or small orifice is used
    dpm.useM1(helpers::readFromCommandLine(argc, argv, "-useM1"));
    // call the solve routine
    dpm.solve();
    return 0;
}
