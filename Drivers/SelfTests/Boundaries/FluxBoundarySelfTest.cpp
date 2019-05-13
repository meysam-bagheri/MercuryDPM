#include "Mercury2D.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/DeletionBoundary.h"
#include "Boundaries/FluxBoundary.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"

FILE * flux_f; 

class FluxBoundarySelfTest : public Mercury2D
{
    public:
        FluxBoundarySelfTest()
        {
            LinearViscoelasticFrictionSpecies* speciesP;

            setTimeStep(1e-3);
            setTimeMax(1.5);
            setSaveCount(100);
            
            speciesP = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
            speciesP->setDensity(1);

            speciesP->setCollisionTimeAndRestitutionCoefficient(
                    5e-1,
                    0.5,
                    (4.0/3.0)*constants::pi*pow(0.05, 3)*1
            );

            speciesP->setSlidingFrictionCoefficient(tan(34 * constants::pi / 180));
            speciesP->setSlidingStiffness(2.0/7.0 * speciesP->getStiffness());
            speciesP->setSlidingDissipation(2.0/7.0 * speciesP->getDissipation());
            speciesP->setRollingFrictionCoefficient(0);
            speciesP->setRollingStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setRollingDissipation(2.0/5.0 * speciesP->getDissipation());
            speciesP->setTorsionFrictionCoefficient(0);
            speciesP->setTorsionStiffness(2.0/5.0 * speciesP->getStiffness());
            speciesP->setTorsionDissipation(2.0/5.0 * speciesP->getDissipation());

            insbA = new CubeInsertionBoundary;
            insbB = new CubeInsertionBoundary;
            BaseParticle* insertionBoundaryParticle = new SphericalParticle; // Possibly evil!
            insertionBoundaryParticle->setSpecies(speciesP);
            insbA->set(
                    insertionBoundaryParticle,
                    1,
                    Vec3D(-1, 0, 0),
                    Vec3D(-1, 0.2, 0),
                    Vec3D(3,0,0),
                    Vec3D(4,0,0),
                    0.05 * 0.9, 
                    0.05 * 1.1 
            );
            insbA = boundaryHandler.copyAndAddObject(insbA);
            insbB->set(
                    insertionBoundaryParticle,
                    1,
                    Vec3D(1, 0.0, 0),
                    Vec3D(1, 0.2, 0),
                    Vec3D(-5,0,0),
                    Vec3D(-4,0,0),
                    0.05 * 0.9, 
                    0.05 * 1.1 
            );
            insbB = boundaryHandler.copyAndAddObject(insbB);
            
            delbA = new DeletionBoundary;
            delbA->set(Vec3D(-1,0,0), 1.1);
            delbA = boundaryHandler.copyAndAddObject(delbA);
            delbB = new DeletionBoundary;
            delbB->set(Vec3D(1,0,0), 1.1);
            delbB = boundaryHandler.copyAndAddObject(delbB);

            fluxb = new FluxBoundary;
            fluxb->set(Vec3D(1,0,0), 0);
            fluxb = boundaryHandler.copyAndAddObject(fluxb);

            /* Output options */
            setXBallsAdditionalArguments("-cmode 7 -v0");
            setXMin(0);
            setXMax(1);
            setYMin(-1);
            setYMax(1);
            setZMin(-1);
            setZMax(1);

            dataFile.setFileType(FileType::ONE_FILE);
            restartFile.setFileType(FileType::ONE_FILE);
            fStatFile.setFileType(FileType::ONE_FILE);

            setGravity(Vec3D(0, 0, 0));
        }

        ~FluxBoundarySelfTest() override = default;

    void actionsAfterTimeStep() override
        {
            fprintf(flux_f, "%f %f %d %f %f %d %f %f %d %f %f %d %f %f\n", 
                    getTime(),
                    particleHandler.getMass(),
                    insbA->getNumberOfParticlesInserted(),
                    insbA->getMassOfParticlesInserted(),
                    insbA->getVolumeOfParticlesInserted(),
                    insbB->getNumberOfParticlesInserted(),
                    insbB->getMassOfParticlesInserted(),
                    insbB->getVolumeOfParticlesInserted(),
                    fluxb->getNumberOfParticlesCrossedForw(),
                    fluxb->getMassOfParticlesCrossedForw(),
                    fluxb->getVolumeOfParticlesCrossedForw(),
                    fluxb->getNumberOfParticlesCrossedBack(),
                    fluxb->getMassOfParticlesCrossedBack(),
                    fluxb->getVolumeOfParticlesCrossedBack() );
        }

    private:
        CubeInsertionBoundary* insbA;
        CubeInsertionBoundary* insbB;
        DeletionBoundary* delbA;
        DeletionBoundary* delbB;
        FluxBoundary* fluxb;
};

int main(int argc, char *argv[])
{
    flux_f = fopen("FluxBoundarySelfTest.flux", "w");
    fprintf(flux_f, "time mass insbA_num insbA_mass insbA_vol insbB_num insbB_mass insbB_vol ");
    fprintf(flux_f, "flux_num_forw flux_mass_forw flux_vol_forw flux_num_back flux_mass_back flux_vol_back\n");
    setbuf(flux_f, nullptr);
    FluxBoundarySelfTest problem;
    problem.setName("FluxBoundarySelfTest");
    fprintf(stdout, "Initialising the problem\n");
    problem.solve();
    fclose(flux_f);
}
