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

#include "Mercury3D.h"
#include "Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h"
#include "Walls/InfiniteWall.h"
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Boundaries/CubeDeletionBoundary.h>

/* This code reads the information generated by LB_S1_insertion.
 * A plane is located at the top surface to compress the particles.
*/

// Main class:
class S2_Compression : public Mercury3D{

private:
    InfiniteWall* lid;

    Mdouble meanCoordinationNumber = 0.0;
    Mdouble scalarNormalForce = 0.0;

    Mdouble newHeight = 0.0;
    Mdouble pressure_ = 0.0; //pressure on the lid [N/m]

    bool moveUpward = false;

public:

    //Constructor
    S2_Compression(Mdouble setPressure)
    {
        pressure_ = setPressure;
        lid = new InfiniteWall;

//        particleSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
//        particleSpecies = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
//        lidSpecies = new ThermalSinterLinFrictionReversibleAdhesiveSpecies(*particleSpecies);
    }

    void setupInitialConditions() override
    {
        lid->setSpecies(speciesHandler.getObject(1));
        lid->set(Vec3D(0.0,0.0,1.0), Vec3D(0.0,0.0,1.2*getZMax()));
        lid = wallHandler.copyAndAddObject(lid);

        logger(INFO, "NumWalls=%",wallHandler.getNumberOfObjects());
    }

    //override continueSolve function such that the code stops
    //when the packing is relaxed (Ekin<1e-2*Eela) and
    //the pressure on the lid is correct (|p/p_lid-1|<1e-1)
    //override continueSolve function such that the code stops
    //when the packing is relaxed (Ekin<1e-5*Eela) and
    //the pressure on the lid is correct (|p/p_lid-1|<1e-3)
    bool continueSolve() const override
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;

            static Mdouble lidArea = constants::pi * mathsFunc::square(getXMax()- getXMin());
            static Mdouble particleArea = constants::pi * mathsFunc::square(particleHandler.getObject(0)->getRadius());
            static Mdouble stiffness = dynamic_cast<const ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(1))->getLoadingStiffness();
            // amount by which the pressure has to be increased
            Mdouble dPressure = lid->getForce().Z/lidArea - pressure_;
            // amount by which position should be changed to achieve the right pressure
            Mdouble dV = dPressure * particleArea /stiffness/getTimeStep();
            //std::cout << "dP/P" << dPressure/pressure << " Z" << dZ << std::endl;

            if (std::abs(dPressure)<1e-5*pressure_ or  moveUpward)
            {
                logger(INFO,"moving lid upward (Releasing) ...");
                lid->setVelocity(Vec3D(0.0,0.0, -dV/100.0));
//                lid->setVelocity(Vec3D(0.0,0.0, 0.0));

                if(getKineticEnergy()<1e-12*getElasticEnergy()) {
                    printTime();
                    logger(INFO,"Stop simulation ...");
                    return false;
                }
                return true;

            } else
            {
                logger(INFO,"moving lid downward (compression)...");
                lid->setVelocity(Vec3D(0.0,0.0,dV/10.0));
                return true;
            }
        }
        return true;
    }

    //--------------------------------------------------
    // To define the height of the compression plane.
    Mdouble getSurfaceHeight() const
    {
        Mdouble height = 0.0;

        for (const BaseParticle* p : particleHandler)
        {
            double newHeight = p->getPosition().Z;
            if (height<newHeight) height = newHeight;
        }
        return height+particleHandler.getMeanRadius();
    }

    // To compute the mean coordination number.
    void actionsAfterTimeStep() override
    {
        newHeight = getSurfaceHeight();

        //To compute the normal force
        for (auto i : interactionHandler)
        {
            scalarNormalForce += Vec3D::dot(i->getForce(), i->getNormal());
        }

        //To compute the coordination number
        for (int i = particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
        }
        meanCoordinationNumber /= particleHandler.getNumberOfObjects();

        //Minimum height to decrease the lid.
        if(lid->getPosition().Z <= getZMax()+ 0.1*particleHandler.getObject(0)->getRadius())
        {
            moveUpward = true;
        }
    }

    //--------------------------------------------------
    //Function to override the output file with specific parameters
    void writeFstatHeader(std::ostream &os) const override
    {
        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();

        Mdouble volParticle = (4.0/3.0)*constants::pi*(std::pow(particleHandler.getMeanRadius(),3.0));
        Mdouble densityParticle = speciesHandler.getObject(0)->getDensity();

        Mdouble volTotalParticles = particleHandler.getVolume();
        Mdouble volumeFraction = volTotalParticles/volSystem;

        Mdouble massParticle = densityParticle*volParticle;
        Mdouble massTotalParticles = particleHandler.getMass();

        //Mdouble volParticlesPlusVoidInSystem = constants::pi*(std::pow(getXMax()/2,2.0))*getSurfaceHeight();
        //Mdouble bulkDensity = massTotalParticles/volParticlesPlusVoidInSystem;
        //Mdouble TheoDensity = 550.0;

        static Mdouble lidArea = constants::pi * mathsFunc::square(getXMax()- getXMin());
        Mdouble dPressure = lid->getForce().Z/lidArea - pressure_;

        os << getTime() //[1] Current time
           << " " << getTimeMax() //[2] Max time
           << " " << volSystem //[3] Volume of the system
           << " " << getZMax() //[4] Height of the system
           << " " << particleHandler.getNumberOfObjects() //[5] Number of particles inserted
           << " " << particleHandler.getMeanRadius() // [6] Mean radius of a particle
           << " " << volParticle //[7] Mean volume of a particle
           << " " << densityParticle //[8] Mean density of a particle
           << " " << volTotalParticles //[9] Volume of particles inserted
           << " " << volumeFraction //[10] Volume fraction
           << " " << massTotalParticles //[11] total mass of particle inserted
           << " " << getSurfaceHeight() // [12] Max height reached by particles.
           << " " << getKineticEnergy() // [13] kinetic energy
           << " " << getElasticEnergy() //[14] Elastic energy
           << " " << scalarNormalForce //[15] Normal force
           << " " << meanCoordinationNumber // [16] meanCoordinationNumber
           << " " << lid->getPosition().Z //[17] Lid position
           << " " << lidArea // [18] Lid area
           << " " << lid->getForce().Z //[19] Lid force
           << " " <<  dPressure //[20] Differencial of pressure
           << std::endl;
    }

    //--------------------------------------------------
    void printTime() const override
    {
        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();
        static Mdouble lidArea = constants::pi * mathsFunc::square(getXMax()- getXMin());
        Mdouble dPressure = lid->getForce().Z/lidArea - pressure_;

        std::cout
                << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
                << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                << ", No. Particles=" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getNumberOfObjects()
                << ", KineticEnergy=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()
                << ", Force=" << std::setprecision(3) << std::left<< std::setw(3)<< lid->getForce().Z
                //                << ", dPressure=" << std::setprecision(3) << std::left<< std::setw(3)<< dPressure
                << ", lidZ=" << std::setprecision(3) << std::left<< std::setw(3)<< lid->getPosition().Z
                << ", ParticlePositionZmax=" << std::setprecision(3) << std::left<< std::setw(3)<< getSurfaceHeight()
                << ", Volume Fraction=" << particleHandler.getVolume()/volSystem
                << std::endl;
    }
};

// Main Function:
int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    //Read file name:
    std::string readFilename = "S1_PA12_500R_200L";
    //Output file name:
    std::string setFilename = "S2_PA12_500R_200L";
    //Set pressure:
    Mdouble setPressure = 1.0e13;

    //Object:
    S2_Compression oTest(setPressure);

    oTest.setName(readFilename); //Name of the file to read.
    oTest.readRestartFile();
    oTest.setRestarted(false);
    oTest.setName(setFilename); //Name of the output file.

    //Simulation parameters:
    oTest.setTimeMax(0.5);

    oTest.removeOldFiles();
    oTest.setParticlesWriteVTK(false);
//    oTest.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);

    oTest.setFileType(FileType::ONE_FILE);
    oTest.solve();

    return 0;
}