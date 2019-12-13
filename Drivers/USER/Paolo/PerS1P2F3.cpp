//
// Created by paolo on 27/08/19.
//
//Ciao!

#include <iostream>
#include <Mercury3D.h>
#include <Particles/SphericalParticle.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Boundaries/StressStrainControlBoundary.h>
#include "ClusterGenerator.h"

class Compression : public Mercury3D
{
public:



    //! Initializing variables, boxSize, species, radii, walls, particles, timeStep, output, name
    void setupInitialConditions() override {

        /*
         * ----------------------------------------------------------------------
         *                        DEFINING VARIABLES
         * ----------------------------------------------------------------------
         */

        random.randomise();

        // Strain
        isoStrainDot = -1e-8 / getTimeStep(); // Fast for the first part


        std::cout << std::endl << std::endl << "CREATING WALLS" << std::endl << std::endl;
        createWalls();

        std::cout << std::endl << std::endl << "SETTING CLUSTERS SPECIES" << std::endl << std::endl;
        setClusterSpecies();

        std::cout << std::endl << std::endl << "SETTING CLUSTERS POSITIONS" << std::endl << std::endl;
        setClusterPositions();

        std::cout << std::endl << std::endl << "INSERTING CLUSTERS" << std::endl << std::endl;

        for (int i = 0; i < nClusters; ++i) {

            clusterGenerator = new ClusterGenerator;

            clusterGenerator->clusterProperties.random.randomise();

            clusterGenerator->clusterProperties.setClusterId(111);

            clusterGenerator->clusterProperties.setPosition(clusterPositions[i]);

            clusterGenerator->clusterProperties.setNumberOfParticles(nParticles);

            clusterGenerator->clusterProperties.setSizeDispersityParticle(dispersity);

            clusterGenerator->clusterProperties.setRadiusParticle(radiusParticle);

            clusterGenerator->clusterProperties.setParticleSpecies(species);

            clusterGenerator->clusterProperties.doVtkOutput(false);

            clusterGenerator->clusterProperties.doIntStrucOutput(false);

            clusterGenerator->clusterProperties.doAmatOutput(false);

            clusterGenerator->clusterProperties.doOverlOutput(false);

            clusterGenerator->clusterProperties.doCdatOutput(false);

            clusterGenerator->clusterProperties.doEneOutput(false);

            clusterGenerator->clusterProperties.doFStatOutput(false);

            clusterGenerator->clusterProperties.doRestartOutput(false);

            std::cout << "Inserting cluster n: " << i+1 << "/" << nClusters << std::endl;

            clusterGenerator->create();

            particleHandler.copyContentsFromOtherHandler(clusterGenerator->clusterProperties.particleHandler);

        }

        em = 0.418;


        std::cout << std::endl << std::endl << "COMPUTING TOTAL PARTICLE VOLUME" << std::endl << std::endl;

        totalParticleVolume = particleHandler.getVolume();

        t0 = getTime();

        makeIsoCFile();

        stage = 0;

    }

    //! Makes data analysis and writes to isotropic compression file
    void actionsBeforeTimeStep() override
    {
        if(stage == 0){
            LinearPlasticViscoelasticInteraction* interaction;
            Mdouble r;
            Mdouble deltaStar = 0.125;
            Mdouble delta;
            BaseInteractable* P;
            BaseInteractable* I;
            for (int i = 0; i < interactionHandler.getSize(); ++i) {
                interaction =  dynamic_cast<LinearPlasticViscoelasticInteraction*>(interactionHandler.getObject(i));
                P=interactionHandler.getObject(i)->getP();
                I=interactionHandler.getObject(i)->getI();
                r = 2 * particleHandler.getObject(P->getIndex())->getRadius() * particleHandler.getObject(I->getIndex())->getRadius() / ( particleHandler.getObject(P->getIndex())->getRadius() + particleHandler.getObject(I->getIndex())->getRadius() );
                delta = interactionHandler.getObject(i)->getOverlap() / r;
                interaction->setMaxOverlap( 0.5*interaction->getMaxOverlap()*(1 + sqrt( 1 + deltaStar / delta ) ) );
            }
            stage = 1;
        }

        if (stage == 1) {
            makeDataAnalysis();
            moveWalls();

            //dampVelocities();
            if (fmod(getTime() - t0, cdatOutputTimeInterval) < getTimeStep()) {
                writeToIsoCFile();
            }

            if ( uniStress > 5e5 && eM < 0.5 ) {
                writeToIsoCFile();
                writeToIsoCFile();
                stage = 2;
            }
        }

        if (stage == 2) {
            makeDataAnalysis();
            moveWallsVertically();
            //moveParticles();

            //dampVelocities();
            if (fmod(getTime() - t0, cdatOutputTimeInterval) < getTimeStep()) {
                writeToIsoCFile();
            }

            if ( uniStress > 2e6 )
                setTimeMax(getTime());
        }
    }

    //! Prints time and values of interest
    void printTime() const override
    {

        switch (stage)
        {
            case 1: std::cout << "Isotropic compression: ";
                break;

            case 2: std::cout << "Uniaxial compression: ";
                break;

            default: std::cout << "Final values: ";
                break;
        }
        std::cout <<

                  "t= " << std::scientific << std::setprecision(2) << std::setw(9) << getTime() <<
                  ", tMax: " << std::scientific << std::setprecision(2) << std::setw(5) << getTimeMax() <<
                  ", Stress: " << std:: fixed << std::setprecision(2) << std::setw(9) << totalStress <<
                  ", uniStress: " << std:: fixed << std::setprecision(2) << std::setw(9) << uniStress <<
                  ", E_ratio = " << std::scientific << std::setprecision(2) << std::setw(8) << getKineticEnergy()/getElasticEnergy() <<
                  ", cN = " << std::fixed << std::setw(5) << meanCoordinationNumber <<
                  ", vR = " << std::fixed << std::setprecision(3) << voidRatio <<
                  ", Porosity = " << std::fixed << std::setprecision(3) << e <<
                  ", dMin = " << std::scientific << std::setw(7) << std::setprecision(5) << minRelativeOverlap <<
                  ", dMean = " << std:: fixed << std::setw(7) << meanRelativeOverlap <<
                  ", dAvCu = " << std:: fixed << std::setw(7) << averageCubicRelativeOverlap <<
                  ", dMax = " << maxRelativeOverlap

                  << std::endl;
    }

    //! inserts the particles in the simulation domainss
    void setClusterPositions()
    {
        clusterPositions = {
                {-4.12471e-06, 1.50160e-10, -4.04206e-06},
                {1.17580e-05, -1.09613e-06, -2.49174e-06},
                {-1.17524e-05, 1.17537e-05, -1.17564e-05},
                {-1.11593e-05, -1.17501e-05, 4.23077e-06},
                {-2.22116e-07, -1.32469e-06, 1.14160e-05},
                {4.21677e-06, 1.17523e-05, 3.33631e-06},
                {-1.17517e-05, 1.13704e-05, 4.23706e-06},
                {4.24337e-06, -1.17532e-05, -1.17548e-05},
                {-1.17544e-05, -1.17532e-05, -1.17562e-05},
                {1.17504e-05, -1.17523e-05, 9.43951e-06},
                {9.52649e-06, 1.17528e-05, -1.17554e-05}
        };
    }

    //! Sets cluster species
    void setClusterSpecies()
    {
        /*
        speciesVector.reserve(nClusters);
        // Single Cluster Species
        for (int i = 0; i < nClusters; ++i) {
            speciesVector[i] = new LinearPlasticViscoelasticFrictionSpecies;
            speciesVector[i] -> setDensity(densityParticle);
            speciesVector[i] -> setConstantRestitution(true);
            speciesVector[i] -> setCollisionTimeAndRestitutionCoefficient(collisionTimeSmallestMass, restitutionCoefficient, 1);
            speciesVector[i] -> setUnloadingStiffnessMax(species->getLoadingStiffness() * unLoadingStiffnessMax / loadingStiffness);
            speciesVector[i] -> setCohesionStiffness(species-> getLoadingStiffness() * cohesionStiffness / loadingStiffness);
            speciesVector[i] -> setPenetrationDepthMax(penetrationDepthMax);

            speciesVector[i] -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
            speciesVector[i] -> setSlidingStiffness(species -> getLoadingStiffness()*2.0/7.0);
            speciesVector[i] -> setSlidingDissipation(species -> getDissipation()*2.0/7.0);
            speciesVector[i] -> setRollingFrictionCoefficient(slidingRollingCoefficient);
            speciesVector[i] -> setRollingStiffness(species -> getLoadingStiffness()*2.0/7.0);
            speciesVector[i] -> setRollingDissipation(species -> getDissipation()*2.0/7.0);
            speciesVector[i] -> setTorsionFrictionCoefficient(slidingTorsionCoefficient);
            speciesVector[i] -> setTorsionStiffness(species -> getLoadingStiffness()*2.0/7.0);
            speciesVector[i] -> setTorsionDissipation(species -> getDissipation()*2.0/7.0);
            speciesHandler.addObject(speciesVector[i]);
            speciesVector.push_back(speciesVector[i]);

        }
        std::cout << "Set single cluster species" << std::endl;

        // Cluster mixed species
        for (int i = 0; i < nClusters; ++i) {
            for (int j = 0; j < nClusters; ++j) {
                if (i!=j){
                    std::cout << i << std::endl;
                    std::cout << j << std::endl;
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setConstantRestitution(true);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setCollisionTimeAndRestitutionCoefficient(collisionTimeSmallestMass, restitutionCoefficient, 1);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setUnloadingStiffnessMax(species->getLoadingStiffness() * unLoadingStiffnessMax / loadingStiffness);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setCohesionStiffness(species-> getLoadingStiffness() * cohesionStiffness / loadingStiffness);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setPenetrationDepthMax(penetrationDepthMax);

                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setSlidingFrictionCoefficient(slidingFrictionCoefficient);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setSlidingStiffness(species -> getLoadingStiffness()*2.0/7.0);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setSlidingDissipation(species -> getDissipation()*2.0/7.0);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setRollingFrictionCoefficient(slidingRollingCoefficient);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setRollingStiffness(species -> getLoadingStiffness()*2.0/7.0);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setRollingDissipation(species -> getDissipation()*2.0/7.0);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setTorsionFrictionCoefficient(slidingTorsionCoefficient);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setTorsionStiffness(species -> getLoadingStiffness()*2.0/7.0);
                    speciesHandler.getMixedObject(speciesVector[i], speciesVector[j]) -> setTorsionDissipation(species -> getDissipation()*2.0/7.0);
                }
            }
        }
         */
        speciesHandler.addObject(species);
        std::cout << "Per ora solo add object. Qui si dovrebbero inserire le mixed species." << std::endl;

    }

    //! Creating stress strain control walls
    void createWalls()
    {
        wallPositionX = 3.95e-05/2;
        wallPositionY = 3.95e-05/2;
        wallPositionZ = 3.95e-05/2;

        wall.set(Vec3D(0, 0, 1), -wallPositionZ, wallPositionZ);
        WNS = boundaryHandler.copyAndAddObject(wall);

        wall.set(Vec3D(1, 0, 0), -wallPositionX, wallPositionX);
        WEW = boundaryHandler.copyAndAddObject(wall);

        wall.set(Vec3D(0, 1, 0), -wallPositionY, wallPositionY);
        WFB = boundaryHandler.copyAndAddObject(wall);

    }

    //! Computes values of interest, such as: coordination number, overlaps, volumeBox, stress
    void makeDataAnalysis()
    {
        //! Computes coordination number and overlaps.
        computeCoordinationNumberAndOverlaps();

        //! Computes void ratio (and volumeBox)
        computeVoidRatioAndPorosity();

        //! Computes the stress.
        computeStress();
    }

    void moveWalls()
    {
        wallPositionX *= 1+isoStrainDot*getTimeStep();
        wallPositionY *= 1+isoStrainDot*getTimeStep();
        wallPositionZ *= 1+isoStrainDot*getTimeStep();

        WNS->set(Vec3D(0, 0, 1), -wallPositionZ, wallPositionZ);

        WEW->set(Vec3D(1, 0, 0), -wallPositionX, wallPositionX);

        WFB->set(Vec3D(0, 1, 0), -wallPositionY, wallPositionY);
    }

    void moveWallsVertically()
    {
        wallPositionZ *= 1+isoStrainDot*getTimeStep();

        WNS->set(Vec3D(0, 0, 1), -wallPositionZ, wallPositionZ);
    }

    //! Computes coordination number and overlaps
    void computeCoordinationNumberAndOverlaps()
    {
        meanForceOnInteraction = 0.0;
        meanCoordinationNumber = 0.0;
        maxRelativeOverlap = 0.0;
        averageCubicRelativeOverlap = 0.0;
        meanRelativeOverlap = 0.0;
        minRelativeOverlap = 2 * ( 1 + dispersity );
        Mdouble relativeOverlap;

        //! loops over each particle to compute mean coordination number and center of mass.
        for (unsigned int i=0; i < particleHandler.getNumberOfRealObjects(); i++)
        {
            meanCoordinationNumber += (particleHandler.getObject(i)->getInteractions()).size();
        }
        meanCoordinationNumber /= particleHandler.getNumberOfRealObjects();

        //! loops over every interaction to compute mean force acting on interaction, maximum, mean and minimum relative particle overlap.
        for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i)
        {
            meanForceOnInteraction += ((*i) -> getForce()).getLength();

            /*!
             * \details the relative overlap is computed as an average of the relative overlap on the two particles.
             *          rO = ( O/R1 + O/R2 ) / 2.
             */
            relativeOverlap = ((*i) -> getOverlap()) / (particleHandler.getObject((*i) -> getP() -> getIndex()) -> getRadius()) +
                              ((*i) -> getOverlap()) / (particleHandler.getObject((*i) -> getI() -> getIndex()) -> getRadius());
            relativeOverlap /= 2;
            meanRelativeOverlap += relativeOverlap;
            averageCubicRelativeOverlap += pow( relativeOverlap, 3 );
            if (relativeOverlap > maxRelativeOverlap)
                maxRelativeOverlap = relativeOverlap;

            if (relativeOverlap < minRelativeOverlap)
                minRelativeOverlap = relativeOverlap;
        }
        meanForceOnInteraction/= interactionHandler.getSize();
        meanRelativeOverlap /= interactionHandler.getSize();
        averageCubicRelativeOverlap = cbrt( averageCubicRelativeOverlap);
        averageCubicRelativeOverlap /= interactionHandler.getSize();
    }

    //! computes voidRatio and porosity (e)
    void computeVoidRatioAndPorosity()
    {
        //volumeBox = (getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());
        volumeBox = 8 * wallPositionX * wallPositionY * wallPositionZ,
                voidRatio = 1 - totalParticleVolume / volumeBox;
        e = (volumeBox - totalParticleVolume) / totalParticleVolume;
        eM = (e-em)/(1+em);
    }

    //! computes stress
    void computeStress()
    {
        stressXX = 0.0;
        stressYY = 0.0;
        stressZZ = 0.0;
        stressXY = 0.0;
        stressXZ = 0.0;
        stressYZ = 0.0;

        for (std::vector<BaseInteraction*>::const_iterator i = interactionHandler.begin(); i != interactionHandler.end(); ++i) {
            stressXX += (*i)->getForce().X * (*i)->getNormal().X * (*i)->getDistance();
            stressYY += (*i)->getForce().Y * (*i)->getNormal().Y * (*i)->getDistance();
            stressZZ += (*i)->getForce().Z * (*i)->getNormal().Z * (*i)->getDistance();
            stressXY += (*i)->getForce().X * (*i)->getNormal().Y * (*i)->getDistance();
            stressXZ += (*i)->getForce().X * (*i)->getNormal().Z * (*i)->getDistance();
            stressYZ += (*i)->getForce().Y * (*i)->getNormal().Z * (*i)->getDistance();
        }

        stressXX /= volumeBox;
        stressYY /= volumeBox;
        stressZZ /= volumeBox;
        stressXY /= volumeBox;
        stressXZ /= volumeBox;
        stressYZ /= volumeBox;

        totalStress = (stressXX+stressYY+stressZZ)/3;

        uniStress = stressZZ;

    }

    //! makes a .isoc file in wich all ifo of isotropic compression are saved
    void makeIsoCFile()
    {

        std::ostringstream cdatName;
        cdatName << getName() << ".isoc";

        isoCFile.open(cdatName.str(), std::ios::out);

        isoCFile << "CLUSTER DATA AND INFORMATION" << std::endl << std::endl;
        isoCFile << "radiusParticle: " << std:: scientific << std::setprecision(2) << radiusParticle << std::endl;
        isoCFile << "sizeDispersityParticle: " << std::defaultfloat << dispersity << std::endl;
        isoCFile << "boxSize: " << std::defaultfloat << (getXMax()-getXMin()) << std::endl;
        isoCFile << "initial solid fraction: " << std::defaultfloat << initialSolidFraction << std::endl;
        // If constantRestitution(true) loading, unloading, and cohesion stiffness are multiplied by the mass of a particle whose radius is radiusParticle_*(1-sizeDispersity),
        // which is the mass that should be used to compute collision time.

        isoCFile << "COMPRESSION INFORMATION" << std::endl << std::endl;
        isoCFile <<

                 "   Time" << std::setw(13) <<
                 "Time Max" << std::setw(16) <<
                 "Stress" << std::setw(16) <<
                 "uniStress" << std::setw(15) <<
                 "E_ratio" << std::setw(12) <<
                 "cN" << std::setw(12) <<
                 "vR" << std::setw(12) <<
                 "e" << std::setw(14) <<
                 "em" << std::setw(12) <<
                 "eM" << std::setw(14) <<
                 "dMin" << std::setw(14) <<
                 "dMean" << std::setw(14) <<
                 "dMax"

                 << std::endl;
    }

    //! writes on the .isoc file
    void writeToIsoCFile()
    {
        isoCFile <<

                 std::scientific << std::setprecision(2) << std::setw(7) << getTime() <<
                 std::scientific << std::setprecision(2) << std::setw(13) << getTimeMax() <<
                 std::fixed << std::setprecision(2) << std::setw(16) << totalStress <<
                 std::fixed << std::setprecision(2) << std::setw(16) << uniStress <<
                 std::scientific << std::setprecision(2) << std::setw(15) << getKineticEnergy()/getElasticEnergy() <<
                 std::fixed << std::setw(12) << meanCoordinationNumber <<
                 std::fixed << std::setprecision(3) << std::setw(12) << voidRatio <<
                 std::fixed << std::setprecision(3) << std::setw(12) << e <<
                 std::fixed << std::setprecision(3) << std::setw(12) << em <<
                 std::fixed << std::setprecision(3) << std::setw(12) << eM <<
                 std::fixed << std::setprecision(5) << std::setw(14) << minRelativeOverlap <<
                 std::setw(14) << meanRelativeOverlap <<
                 std::setw(14) << maxRelativeOverlap

                 << std::endl;

    }

    //! Moves particles in X and Y direction during stage 3
    void moveParticles(){
        Vec3D centerBox;
        centerBox = Vec3D( (getXMax()+getXMin())/2 , (getYMax()+getYMin())/2 , (getYMax()+getYMin())/2 );
        for (auto& p : particleHandler)
        {
            relativeToCenter.X = p->getPosition().X - centerBox.X;
            relativeToCenter.Y = p->getPosition().Y - centerBox.Y;
            relativeToCenter.Z = p->getPosition().Z - centerBox.Z;
            p->move( Vec3D(-isoStrainDotFast * getTimeStep() * relativeToCenter.X,
                           -isoStrainDotFast * getTimeStep() * relativeToCenter.Y,
                           0) );
        }
    }

    /*
     * -------------------------- VARIABLES ------------------------------------
     */


    // Simulation
    Mdouble cdatOutputTimeInterval;

    // domain
    Mdouble volumeBox;
    Mdouble initialSolidFraction;


    // Particles
    Mdouble radiusParticle;
    Mdouble radiusCluster;
    Mdouble volumeParticle;
    Mdouble volumeCluster;
    int nParticles;
    int nClusters;
    Mdouble dispersity;

    // species
    LinearPlasticViscoelasticFrictionSpecies* species;
    Mdouble collisionTimeSmallestMass;

    Mdouble densityParticle;

    Mdouble loadingStiffness;
    Mdouble unLoadingStiffnessMax;
    Mdouble cohesionStiffness;
    Mdouble slidingFrictionCoefficient;
    Mdouble slidingRollingCoefficient;
    Mdouble slidingTorsionCoefficient;
    Mdouble restitutionCoefficient;
    Mdouble penetrationDepthMax;


private:

    // Simulation
    int stage;
    Mdouble t0;

    // Domain
    Mdouble isoStrainDotFast;
    Mdouble isoStrainDot;

    //Clusters & particles
    Mdouble totalParticleVolume;
    std::vector<Mdouble> radii;
    std::vector<Vec3D> clusterPositions;
    Vec3D relativeToCenter;

    // Species
    std::vector<LinearPlasticViscoelasticFrictionSpecies*> speciesVector;

    // Walls
    Mdouble wallPositionX;
    Mdouble wallPositionY;
    Mdouble wallPositionZ;
    PeriodicBoundary wall;
    PeriodicBoundary* WNS;
    PeriodicBoundary* WEW;
    PeriodicBoundary* WFB;

    // Stress
    Mdouble stressXX;
    Mdouble stressYY;
    Mdouble stressZZ;
    Mdouble stressXY;
    Mdouble stressXZ;
    Mdouble stressYZ;
    Mdouble totalStress;
    Mdouble uniStress;

    // Data analysis
    Mdouble  meanForceOnInteraction;
    Mdouble  meanCoordinationNumber;
    Mdouble  maxRelativeOverlap;
    Mdouble  averageCubicRelativeOverlap;
    Mdouble  meanRelativeOverlap;
    Mdouble  minRelativeOverlap;
    Mdouble  voidRatio;
    Mdouble  e;
    Mdouble  em;
    Mdouble  eM;

    // Output
    //!\brief cluster data file.
    std::ofstream isoCFile;



};



int main(){

    Compression script;

    script.dataFile.setFileType(FileType::ONE_FILE);
    script.restartFile.setFileType(FileType::ONE_FILE);
    script.fStatFile.setFileType(FileType::NO_FILE);
    script.eneFile.setFileType(FileType::NO_FILE);
    logger(INFO, "run number: %", script.dataFile.getCounter());

    Mdouble boxSize = 3.95e-05; // Le posizioni delle particelle sono sul file info.txt
    script.setXMin(-0.5*boxSize);
    script.setYMin(-0.5*boxSize);
    script.setZMin(-0.5*boxSize);
    script.setXMax( 0.5*boxSize);
    script.setYMax( 0.5*boxSize);
    script.setZMax( 0.5*boxSize);


    // Particles properties
    script.nParticles = 300;
    script.dispersity = 0.1;
    script.radiusParticle =1e-6;
    script.volumeParticle = 4 * constants::pi * pow(script.radiusParticle, 3) / 3;
    script.radiusCluster = 8 * script.radiusParticle; // Empiricamente: SETTARLO CON SICUREZZA!!!!!!!!!!
    script.volumeCluster = 4 * constants::pi * pow(script.radiusCluster, 3) / 3;

    // Computing the number of particles (if every particle had radius = radiusParticle)
    script.initialSolidFraction = 0.3;
    script.volumeBox = pow( script.getXMax() - script.getXMin() - 2*script.radiusCluster, 3 ); // settato così per evitare overlap iniziali
    //script.nClusters = floor( script.volumeBox * script.initialSolidFraction / script.volumeCluster );
    script.nClusters = 11;

    script.densityParticle = 2825; // ottenuta considerando la densità della bentonite = 2150 Kg/m^3 ed una vF di 0.37 (corrispondente al momento in cui E_ratio crolla)

    script.loadingStiffness = 1e3;
    script.unLoadingStiffnessMax = 5 * script.loadingStiffness;
    script.cohesionStiffness = 0.5 * script.loadingStiffness;

    script.slidingFrictionCoefficient = 0.5;
    script.slidingRollingCoefficient = 0.5 * script.slidingFrictionCoefficient;
    script.slidingTorsionCoefficient = 0.0;

    script.restitutionCoefficient = 0.5;

    script.penetrationDepthMax = 0.2;


    Mdouble smallestMass = script.densityParticle * 4 * constants::pi * pow(script.radiusParticle*(1-script.dispersity), 3) / 3;
    script.collisionTimeSmallestMass = sqrt( smallestMass * ( pow(constants::pi, 2) + pow(log(script.restitutionCoefficient), 2) ) / ( 2 * script.loadingStiffness ) );


    // SINGLE_PARTICLE-SINGLE_PARTICLE INTRACLUSTER (set here, added in class)

    script.species = new LinearPlasticViscoelasticFrictionSpecies;
    script.species -> setDensity(script.densityParticle);
    script.species -> setConstantRestitution(true);
    script.species -> setCollisionTimeAndRestitutionCoefficient(script.collisionTimeSmallestMass, script.restitutionCoefficient, 1);
    script.species -> setUnloadingStiffnessMax(script.species->getLoadingStiffness() * script.unLoadingStiffnessMax / script.loadingStiffness);
    script.species -> setCohesionStiffness(script.species-> getLoadingStiffness() * script.cohesionStiffness / script.loadingStiffness);
    script.species -> setPenetrationDepthMax(script.penetrationDepthMax);

    script.species -> setSlidingFrictionCoefficient(script.slidingFrictionCoefficient);
    script.species -> setSlidingStiffness(script.species -> getLoadingStiffness()*2.0/7.0);
    script.species -> setSlidingDissipation(script.species -> getDissipation()*2.0/7.0);
    script.species -> setRollingFrictionCoefficient(script.slidingRollingCoefficient);
    script.species -> setRollingStiffness(script.species -> getLoadingStiffness()*2.0/7.0);
    script.species -> setRollingDissipation(script.species -> getDissipation()*2.0/7.0);
    script.species -> setTorsionFrictionCoefficient(script.slidingTorsionCoefficient);
    script.species -> setTorsionStiffness(script.species -> getLoadingStiffness()*2.0/7.0);
    script.species -> setTorsionDissipation(script.species -> getDissipation()*2.0/7.0);



    Mdouble contactTimeOverTimeStep = 51;
    std::cout << "collision time: " << script.species -> getCollisionTime(smallestMass) << std::endl;
    script.setTimeStep(script.species -> getCollisionTime(smallestMass)/contactTimeOverTimeStep);
    std::cout << "timeStep: " << std::setprecision(4) << script.getTimeStep() << std::endl;
    std::cout << "tC_PP/dt, at least: " << std::setprecision(4) << script.species -> getCollisionTime(smallestMass)/script.getTimeStep() << std::endl << std::endl;


    script.cdatOutputTimeInterval = script.getTimeStep() * 5000;
    script.setSaveCount( floor(script.cdatOutputTimeInterval/script.getTimeStep()) );
    script.setTimeMax(1000);
    script.setXBallsAdditionalArguments("-v0 -p 10");

    // NAME SETTING
    std::ostringstream name;
    name << "PerS1P2F3";
    script.setName(name.str());


    script.solve();

    return 0;
}




