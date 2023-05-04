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
//#include "CG/CG.h"

/* This code reads the information generated by S2_Compression file.
 * The temperature increases and heats up all particles. Thus, the Sintering process
 * occurs.
*/
// Main class:
class S3_HomogeneousSintering : public Mercury3D{

private:
    Mdouble thermalConductivity_ = 0.0;
    Mdouble heatCapacity_ = 0.0;

    Mdouble T0_ = 0.0;
    Mdouble maxTemp_ = 0.0;
    Mdouble gradientTemp_ = 0.0; //[C/s] Heating rate
    Mdouble holdingTime_ = 0.0;
    Mdouble meltingTemperature_ = 0.0;

    Mdouble meanCoordinationNumber = 0.0;
    Mdouble deltaR_ = 0.0; //expansion coefficient

    Mdouble scalarNormalForce = 0.0;
    Mdouble newHeight = 0.0;

    Mdouble deltaC_ = 0.0;
    Mdouble fluidity_ = 0.0;
    Mdouble compliance0_ = 0.0;
    Mdouble surfaceTension_ = 0.0;
    Mdouble sinterAdhesion_ = 0.0;

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies;
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint;

public:
    //Constructor
    S3_HomogeneousSintering(Mdouble deltaR, Mdouble startingTemperature, Mdouble  maxTemp, Mdouble gradientTemp,
                            Mdouble meltingTemperature, Mdouble holdingTime, Mdouble setHCapacity, Mdouble setTConductivity)
    {
        deltaR_ = deltaR;

        maxTemp_ = maxTemp;
        T0_ = startingTemperature;
        gradientTemp_ = gradientTemp;
        meltingTemperature_ = meltingTemperature;
        holdingTime_ = holdingTime;

        thermalConductivity_ = setTConductivity;
        heatCapacity_ = setHCapacity;

        //Read file from S2_Compression
        setName("S2_PA12_1000R_350L");
        readRestartFile();
        //Set output file
        setRestarted(false);
        setName("S3_H_PA12_1000R_350L");

        //Remove Lid from S2_Compression
//        wallHandler.removeLastObject();

        //-------------------------------------------------
        //Get the particle species to set extra parameters to the particles
        particleSpecies = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
        //-------------------------------------------------
        wallHandler.getLastObject()->setVelocity(Vec3D(0.0,0.0, 0.00));
//        wallHandler.getLastObject()->setSpecies(particleSpecies);

        //Set a "virtual" species to allow change particle properties close to the melting point
        particleSpeciesAtMeltingPoint = new ThermalSinterLinFrictionReversibleAdhesiveSpecies(*particleSpecies);
    }
    //--------------------------------------------------
    //Destructure
    ~S3_HomogeneousSintering() override
    {
        delete particleSpeciesAtMeltingPoint;
    }

    //--------------------------------------------------
    //Set sintering parameters
    void sinteringParameters(Mdouble deltaC,Mdouble fluidity,Mdouble compliance0,Mdouble surfaceTension, Mdouble sinterAdhesion)
    {
        deltaC_ = deltaC;
        fluidity_ = fluidity;
        compliance0_ = compliance0;
        surfaceTension_ = surfaceTension;
        sinterAdhesion_ = sinterAdhesion;
    }
    //--------------------------------------------------
    //Initial conditions
    void setupInitialConditions() override
    {
//        Mdouble sinterAdhesion = 0.01*particleSpecies->getLoadingStiffness();
        //Set thermal parameters
        particleSpecies->setThermalConductivity(thermalConductivity_);
        particleSpecies->setHeatCapacity(heatCapacity_);
        particleSpecies->setSinterAdhesion(sinterAdhesion_);

        particleSpecies->setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT

        particleSpecies->setPenetrationDepthMax(1.5); //Allow particles penetrate each other
        particleSpecies->setComplianceZero(compliance0_); //Book: Adhesive Particles
        particleSpecies->setSurfTension(surfaceTension_); //

        //To control:
        particleSpecies->setFluidity(fluidity_);
        particleSpecies->setSeparationDis(deltaC_);

        //+++++++++++++++++++++++++++++++++++++
        particleSpeciesAtMeltingPoint->setLoadingStiffness(particleSpecies->getLoadingStiffness());
        particleSpeciesAtMeltingPoint->setSinterAdhesion(particleSpecies->getSinterAdhesion());
        particleSpeciesAtMeltingPoint->setSurfTension(particleSpecies->getSurfTension());
        particleSpeciesAtMeltingPoint->setFluidity(particleSpecies->getFluidity());


        //Set initial temperature to all particles
        for (BaseParticle* p : particleHandler)
        {
            auto* tp = dynamic_cast<ThermalParticle*>(p);
            tp->setTemperature(T0_);
            tp->setVelocity(Vec3D(0.0,0.0,0.0));
        }
        logger(INFO, "Set uniform temperature to all particles [C]= %", T0_);

        if(holdingTime_>getTimeMax()){

            logger(INFO, "*******************************");
            logger(INFO, "Set correctly the holding time");
            logger(INFO, "*******************************");
            exit(-1);
        }
    }
    //--------------------------------------------------
    Mdouble getSurfaceHeight() const
    {
        double height = 0;

        for (const BaseParticle* p : particleHandler)
        {
            double newHeight = p->getPosition().Z;
            if (height<newHeight) height = newHeight;
        }
        return height+particleHandler.getMeanRadius();
    }

    //--------------------------------------------------
    Mdouble getTemperature()const
    {
        Mdouble heatingTemperature = T0_ + gradientTemp_*getTime();
        Mdouble coolingTemperature = maxTemp_;
        Mdouble isTimeForColling =  getTimeMax()- holdingTime_;

        if(heatingTemperature >= maxTemp_)
        {
             heatingTemperature = maxTemp_;

             if(getTime() >= isTimeForColling){
                 coolingTemperature = T0_ + (maxTemp_-T0_)/(getTimeMax()-isTimeForColling)*(getTimeMax() - getTime());

             }

        }
        return std::min(std::min(heatingTemperature,coolingTemperature),maxTemp_);
    }

    //--------------------------------------------------
    void setTemperature(Mdouble temperature)
    {
        //Variables to update
        static Mdouble oldTemperature = T0_;
        Mdouble deltaTemperature = temperature-oldTemperature;
        Mdouble deltaCooling = temperature - maxTemp_;

        //add temperature above a certain height
//        Mdouble lb_z = getSurfaceHeight() - 2.0 * particleHandler.getMeanRadius();

        for (BaseParticle* p : particleHandler)
        {
            auto *tp = dynamic_cast<ThermalParticle *>(p);

            if (p->getSpecies()==particleSpecies)
            {
                tp->setTemperature(temperature); //Change the temperature of the particle
            }
        }
        //Set changes of particles heated close to the melting point
        setHeatEffect(temperature,deltaTemperature,deltaCooling);

        //Update temperature
        oldTemperature = temperature;
    }

    //Function to change particle properties based on the temperature
    void setHeatEffect(Mdouble temperature,Mdouble deltaTemperature, Mdouble deltaCooling)
    {
        //Thermal expansion
        Mdouble factorRadius = 1.0 + (deltaR_ * deltaTemperature);

        for (BaseParticle* p : particleHandler)
        {
            if (p->getSpecies()==particleSpecies)
            {
                auto* tp = dynamic_cast<ThermalParticle*>(p);
                Mdouble pT = tp->getTemperature();
                //Check which particles are 95% closer to the melting point, then the particle radius can decrease
                if(pT >= 0.95*meltingTemperature_)
                {
                    p->setRadius((factorRadius * p->getRadius()));
                }
            }
        }

        //change species properties
        Mdouble stiffnessFactor = 0.5*(1.0+tanh((meltingTemperature_-temperature)/gradientTemp_));
        Mdouble adhesionFactor = 0.5*(1.0-tanh((meltingTemperature_-temperature)/gradientTemp_));
        Mdouble complianceFactor = (2.0-tanh((maxTemp_-temperature)/gradientTemp_));
        Mdouble oldLoadingStiffness = particleSpecies->getLoadingStiffness();

        particleSpecies->setLoadingStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint->getLoadingStiffness());
//        particleSpecies->setSinterAdhesion(adhesionFactor*particleSpeciesAtMeltingPoint->getSinterAdhesion());
//        particleSpecies->setFluidity(complianceFactor*particleSpeciesAtMeltingPoint->getFluidity());
//        particleSpecies->setSurfTension(stiffnessFactor*particleSpeciesAtMeltingPoint->getSurfTension());

        //for decreasing temperature, change the maxOverlap
        if (deltaCooling<0.0)
        {
            for (BaseInteraction* cBase : interactionHandler)
            {
                auto c =  dynamic_cast<SinterLinInteraction*>(cBase);
                Mdouble unloadingStiffness = c->getUnloadingStiffness();
                c->setMaxOverlap(c->getMaxOverlap()
                                 *(unloadingStiffness-oldLoadingStiffness)
                                 /(unloadingStiffness-particleSpecies->getLoadingStiffness())
                );
            }
        }
    }
    //--------------------------------------------------
    void actionsAfterTimeStep() override
    {
        setTemperature(getTemperature());
        newHeight = getSurfaceHeight();


        //To compute the normal force
        for (auto i : interactionHandler)
        {
            scalarNormalForce += Vec3D::dot(i->getForce(),i->getNormal());
        }

        //To measure the mean coordination number.
        for (int i = particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
        }
        meanCoordinationNumber /= particleHandler.getNumberOfObjects();
    }

    //System
    Mdouble getThermalEnergy() const
    {
        Mdouble eneEnergy;

        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;

            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            eneEnergy += (tp->getTemperature()-T0_) * tp->getMass() * particleSpecies->getHeatCapacity();
//            eneEnergy += deltaTemperature * tp->getMass() * particleSpecies->getHeatCapacity();
        }
        return eneEnergy;
    }

    //--------------------------------------------------
    //Prints the temperature into the data file, such that you can plot it in paraview
    Mdouble getInfo(const BaseParticle& p) const override
    {
        return (dynamic_cast<const ThermalParticle &>(p).getTemperature()-T0_)/(maxTemp_-T0_);
    }

    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();
        return sqrt(meanOverlap/meanRadius); //Todo: Check if it is mean Radius or initial radius
    }
    //--------------------------------------------------
    //Function to override the output file with specific parameters
    void writeFstatHeader(std::ostream &os) const override
    {
        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();
        Mdouble volSintering = constants::pi*(std::pow(getXMax()/2,2))*getSurfaceHeight();

        Mdouble volParticle = (4.0/3.0)*constants::pi*(std::pow(particleHandler.getMeanRadius(),3.0));
        Mdouble densityParticle = speciesHandler.getObject(0)->getDensity();

        Mdouble volTotalParticles = particleHandler.getVolume();
        Mdouble volumeFraction = volTotalParticles/volSystem;

        Mdouble massParticle = densityParticle*volParticle;
        Mdouble massTotalParticles = particleHandler.getMass();

        Mdouble meanRadius = particleHandler.getMeanRadius();
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;

        Mdouble volParticlesPlusVoidInSystem = constants::pi*(std::pow(getXMax()/2,2.0))*getSurfaceHeight();
        Mdouble bulkDensity = massTotalParticles/(volParticlesPlusVoidInSystem);

        Mdouble X = getMeanRelativeContactRadius()/2;
        Mdouble rho_p = 1.0/(1 + std::pow(1 - X,3.0));

        Mdouble K_s = 2.0*meanRadius/(std::pow(meanContactRadius,2.0));
        Mdouble Shrinkage = (1.0/5.0)*std::pow(getMeanRelativeContactRadius(),2.0);

        Mdouble tangentialStress = K_s * particleSpecies->getSurfTension();

        Mdouble rateContactGrowth = getMeanRelativeContactRadius()/getTime();

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
           << " " << getTemperature()//[17]Current temperature
           << " " << particleSpecies->getLoadingStiffness()// [18] Loading stiffness
           << " " << particleSpeciesAtMeltingPoint->getLoadingStiffness()// [19] Loading stiffness
           << " " << getMeanRelativeContactRadius()//[20] Neckgrowth
           << " " << meanOverlap//[21]
           << " " << meanContactRadius //[22]
           << " " << getThermalEnergy() //[23]
           << " " << volSintering //[24]
           << " " << bulkDensity //[25]
           << " " << particleSpecies->getSinterAdhesion()// [26] sintering Adhesion
           << " " << particleSpecies->getFluidity()// [27] fluidity
           << " " << particleSpecies->getSurfTension()// [28] fluidity
           << " " << rho_p //29
           << " " << Shrinkage //30
           << " " << tangentialStress //31
           << " " << rateContactGrowth //32
           << std::endl;
    }
    //--------------------------------------------------
    //Function to display data at the console
    void printTime() const override
    {
        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();
        Mdouble volTotalParticles = particleHandler.getVolume();

        std::cout << "t=" << std::setprecision(3)<< std::left << std::setw(3) << std::left<< getTime()
                  << ", tmax= " <<  std::setprecision(3) << std::left << std::setw(3) << getTimeMax()
                  << ", lidPosition=" << wallHandler.getLastObject()->getPosition().Z
                  << ", particlePosition= "<< getSurfaceHeight()
                  << ", Temp=" << std::setw(3)<< std::left << std::setw(3) << getTemperature()
                  << ", KineticEnergy=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()
                  << ", ThermalEnergy="<< std::setw(3)<< std::left << std::setw(3) << getThermalEnergy()
                  << ", loading stiffness= " << std::setw(3)<< std::left << std::setw(3) << particleSpecies->getLoadingStiffness()
                  << ", MeanOverlap= "<< std::setw(3)<< std::left << std::setw(3) << getMeanRelativeContactRadius()
                  << std::endl;
    }
};

// Main function:
int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Material information
    Mdouble deltaR = 0.0001; //Expansion coefficient close to the melting point
    Mdouble heatCapacity = 1185.0; //[J/(kgK)]
    Mdouble ThermalConductivity = 0.230;  //0.25[W/(mK)] //Todo:Check

    Mdouble E_Modulus = 1.650e9; //[Pa] Young's Modulus of polyamide12

    // define temperature profile
    Mdouble meltingTemperature = 185.0; // [C] melting point of the material
    Mdouble startingTemperature = 170.0; // [C] initial temperature
    Mdouble gradientTemp = 840.0; //[C/s]
    Mdouble maxTemp = 185.0; // [C] max temperature

    Mdouble holdingTime = 0.5; //[s] After this time temperature starts decreasing
    Mdouble maxTime = 1.0;

    // sintering parameters
    Mdouble deltaC = 1.8e-7; //separation distance [m]
    Mdouble fluidity = 10.5; //fluidity [Pa s]^{-1}
    Mdouble compliance0 = 1.0/(2.0*E_Modulus);

    Mdouble surfaceTension = 0.04;

    Mdouble meanRadius = 3.0e-5; //[mm]//
    Mdouble k1 = 4.0/3.0 * E_Modulus*sqrt(meanRadius)*1e-8;
    Mdouble sinterAdhesion = 0.0001*k1;

    S3_HomogeneousSintering oTest(deltaR,startingTemperature,maxTemp,gradientTemp,meltingTemperature,holdingTime,heatCapacity,ThermalConductivity);

    oTest.sinteringParameters(deltaC, fluidity,compliance0,surfaceTension,sinterAdhesion);

    oTest.setParticlesWriteVTK(true);
    oTest.wallHandler.setWriteVTK(true);

    oTest.setTimeMax(maxTime); //[s]
    oTest.setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");

    //--------------------------------------------------
    //This is for Coarse-Graining live.
    //Live CG to get results over time. For instance, plot stress vs time.
    //    //define coarse-graining resolved in z
    //    CG<CGCoordinates::Z> cgZ; //declare a cg object
    //    cgZ.setN(20);  //set number of mesh points in each spatial direction
    //    cgZ.setWidth(0.5); //set cg width
    //    s.cgHandler.copyAndAddObject(cgZ); // add the CG object to the cgHandler

    oTest.removeOldFiles();
//    oTest.setSaveCount(1000);
//    oTest.setTimeStep(0.001);
    oTest.solve();

    //Coarse-Graining:
    //--------------------------------------------------
    // It creates the coarse-graining output file at the last iteraction.
    logger(INFO,"Execute 'source S3_Sintering.sh' to get coarse-grained statistics at specific time step");
    helpers::writeToFile("S3_Sintering.sh","../../../../../MercuryCG/fstatistics S3_Sintering -stattype XZ -w 1.0e-5 -h 0.1e-5 -tmin 1.0 -tmax 1.1 -o S3_SinteringScaledMass.XZ.stat");

    //--------------------------------------------------
    // It creates the matlab visualization.
    logger(INFO,"Run 'S3_Sintering.m' in MATLAB/octave to visualise the statistical output");
    helpers::writeToFile("S3_Sintering.m","clc;clear all;close all\n"
                                          "addpath('../../../../../MercuryCG/')\n"
                                          "data = loadStatistics('S3_Sintering.XZ.stat');\n"
                                          "colormap(1-gray)\n"
                                          "contourf(data.x,data.z,data.Density,20,'EdgeColor','none')\n"
                                          "c = colorbar\n"
                                          "c.Label.String = '\\rho';\n"
                                          "title('Density')\n"
                                          "xlabel('x')\n"
                                          "ylabel('z');\n"
                                          "axis equal\n"
                                          "%%\n"
                                          "addpath('/home/juan/Softwares/MercuryDPM_Branch/Trunk/Matlab');\n"
                                          "particles=read_data('S3_Sintering.data');\n"
                                          "Cell = particles{637}(1,1); %Specific particle position at 1.119, which is in the cell 24\n"
                                          "NumPart= Cell.N;\n"
                                          "Pradius = Cell.Radius;\n"
                                          "Position = Cell.Position;\n"
                                          "a=linspace(0,2*pi,40);\n"
                                          "xCircle = sin(a);\n"
                                          "zCircle = cos(a);\n"
                                          "hold on;\n"
                                          "for i=1:length(Pradius)\n"
                                          "  plot(Position(i,1) + Pradius(i)*xCircle,Position(i,3)+Pradius(i)*zCircle,'Color',.8*[1 1 1])\n"
                                          "end\n"
                                          "hold off");
    return 0;
}