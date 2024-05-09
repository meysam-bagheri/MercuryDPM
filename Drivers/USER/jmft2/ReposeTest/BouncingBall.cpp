#include "Mercury3D.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"
#include <iostream>
#include <cmath>
#include <map>

#include "muICal3DParameters.h"

#define MAX_STRLEN 1024


class BouncingBall : public Mercury3D {
public:
    BouncingBall(std::string paramsFileName) {
      params = muICal3DParameters::read(paramsFileName);
      theta = 0;
      params->theta = theta;

      setName("BouncingBall");

      std::cout << theta << std::endl;

      double betaslide = params->betaslide;
      double betaroll = params->betaroll;
      double betators = params->betators;

      g = params->g;

      particleRadius = params->particleRadius;
      double rho = params->rho;

      setTimeMax(params->timeMax);
      setTimeStep(params->timeStep);
      setSaveCount(params->saveEvery);

      setXMin(-params->length / 2);
      setXMax(params->length / 2);
      setYMin(-params->width / 2);
      setYMax(params->width / 2);
      setZMin(0);
      setZMax(params->height);

      setNumberOfDomains({4, 1, 1});

      // dataFile.setFileType(FileType::NO_FILE);
      // fStatFile.setFileType(FileType::NO_FILE);
      dataFile.setFileType(FileType::MULTIPLE_FILES);
      fStatFile.setFileType(FileType::MULTIPLE_FILES);

      /* Define species */
      species = speciesHandler.copyAndAddObject(new LinearViscoelasticFrictionSpecies());
      species->setDensity(rho);
      species->setCollisionTimeAndRestitutionCoefficient(
        params->collisionTime,
        params->restitutionCoefficient,
        (4./3.) * constants::pi * pow(particleRadius, 3) * rho
      );
      species->setSlidingFrictionCoefficient(tan(betaslide));
      species->setSlidingStiffness(2.0 / 7.0 * species->getStiffness());
      species->setSlidingDissipation(2.0 / 7.0 * species->getDissipation());
      species->setRollingFrictionCoefficient(tan(betaroll));
      species->setRollingStiffness(2.0 / 5.0 * species->getStiffness());
      species->setRollingDissipation(2.0 / 5.0 * species->getDissipation());
      species->setTorsionFrictionCoefficient(tan(betators));
      species->setTorsionStiffness(2.0 / 5.0 * species->getStiffness());
      species->setTorsionDissipation(2.0 / 5.0 * species->getDissipation());

      logger(INFO, "Maximum collision velocity is %", species->getMaximumVelocity(
        params->particleRadius,
        species->getMassFromRadius(params->particleRadius)
      ));

      /* Walls */
      InfiniteWall bottomwall;
      bottomwall.setSpecies(species);
      bottomwall.set(Vec3D(0, 0, -1), Vec3D(0, 0, 0));
      wallHandler.copyAndAddObject(bottomwall);

      /* Periodic boundaries */
      PeriodicBoundary *xbounds = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
      xbounds->set(Vec3D(1, 0, 0), getXMin(), getXMax());
      PeriodicBoundary *ybounds = boundaryHandler.copyAndAddObject(new PeriodicBoundary);
      ybounds->set(Vec3D(0, 1, 0), getYMin(), getYMax());
    }

    void setupInitialConditions() override {
      setGravity(Vec3D(g * sin(theta), 0, -g * cos(theta)));


      SphericalParticle particle;
      particle.setSpecies(species);
      particle.setVelocity(Vec3D(0, 0, 0));

      particle.setRadius(particleRadius);
      Mdouble gap = particleRadius * (1 + params->dispersity);
      particle.setPosition(Vec3D(0, 0, getZMax()));
      particleHandler.copyAndAddObject(particle);

      logger(INFO, "Particle volume is %, domain volume is %",
             particleHandler.getVolume(), getTotalVolume());
    }

    void writeOutputFiles() override {
      Mercury3D::writeOutputFiles();

      if (eneFile.getLastSavedTimeStep() != getNumberOfTimeSteps()) {
        return;
      }
    }

    void actionsAfterSolve() override {
      dataFile.setFileType(FileType::MULTIPLE_FILES);
      writeDataFile();
      fStatFile.setFileType(FileType::MULTIPLE_FILES);
      writeFStatFile();
    }

private:
    muICal3DParameters *params;
    RNG generator;
    LinearViscoelasticFrictionSpecies *species;
    CubeInsertionBoundary *insb;

    double g;
    double theta;  // slope angle in radians
    double particleRadius;
};


int main(int argc, char **argv) {
  auto problem = new BouncingBall(argv[1]);
  argv[1] = argv[0];
  problem->solve(argc - 1, argv + 1);
  delete problem;
  return 0;
}

