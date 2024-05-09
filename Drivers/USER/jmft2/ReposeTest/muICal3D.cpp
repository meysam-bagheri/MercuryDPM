/* muICal3D - Drop some particles onto a base z=0, sloped at an angle
 * theta(t) to the horizontal. Let theta(t) vary slowly with time, so that the
 * flow may adjust to an equilibrium. Measure I (by measuring u and h) at
 * different times, and plot I against theta.
 * Then obtain the properties mu_1, mu_2 and I_0 by fitting an appropriate
 * curve. (This will be done externally, using Matlab.)
 * */

#include "Mercury3D.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Math/RNG.h"
#include "Math/ExtendedMath.h"
#include "File.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <map>

#include "muICal3DParameters.h"
#include "Mixins/PrintWallTimeMixin.h"

#define MAX_STRLEN 1024


class muICal3D : public Mercury3D, public PrintWallTimeMixin {
public:
    muICal3D(std::string paramsFileName, std::string thetaInDegreesStr) {
      params = muICal3DParameters::read(paramsFileName);
      float thetaInDegrees = std::stof(thetaInDegreesStr);
      theta = thetaInDegrees * DEGREES;
      params->theta = theta;

      setName(paramsFileName.erase(paramsFileName.find_last_of('.')) + "-" + thetaInDegreesStr);

      std::cout << theta << std::endl;

      double betaslide = params->betaslide;
      double betaroll = params->betaroll;
      double betators = params->betators;

      g = params->g;

      particleRadius = params->particleRadius;
      baseRadius = params->baseRadius;
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

      dataFile.setFileType(FileType::NO_FILE);
      fStatFile.setFileType(FileType::NO_FILE);

      /* Define species */
      species = speciesHandler.copyAndAddObject(new LinearViscoelasticFrictionSpecies());
      species->setDensity(rho);
      species->setCollisionTimeAndRestitutionCoefficient(
        params->collisionTime,
        params->restitutionCoefficient,
        (4. / 3.) * constants::pi * pow(particleRadius, 3) * rho
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

      /* Rough base */
      createSpacedGridRoughBase();
    }

    void createSpacedGridRoughBase() {
      if (params->baseConc <= 0) return;
      double spacing = 4 * baseRadius / sqrt(params->baseConc);
      for (double xpos = getXMin(); xpos <= getXMax(); xpos += spacing) {
        for (double ypos = getYMin(); ypos <= getYMax(); ypos += spacing) {
          double zpos = 0;

          SphericalParticle basalParticle;
          basalParticle.setSpecies(species);
          basalParticle.setRadius(baseRadius *
                                  (1 + params->baseDispersity * generator.getRandomNumber(-1, 1)));
          basalParticle.setPosition(Vec3D(xpos, ypos, zpos));
          basalParticle.fixParticle();
          particleHandler.copyAndAddObject(basalParticle);
        }
      }
    }

    void createRandomisedGridRoughBase() {
      if (params->baseConc <= 0) return;
      double spacing = 4 * baseRadius / sqrt(params->baseConc);
      for (double xpos = getXMin(); xpos <= getXMax(); xpos += spacing) {
        for (double ypos = getYMin(); ypos <= getYMax(); ypos += spacing) {
          double zpos = 0;

          double actualXpos = xpos + spacing * generator.getRandomNumber(-0.5, 0.5);
          double actualYpos = ypos + spacing * generator.getRandomNumber(-0.5, 0.5);

          SphericalParticle basalParticle;
          basalParticle.setSpecies(species);
          basalParticle.setRadius(baseRadius *
                                  (1 + params->baseDispersity * generator.getRandomNumber(-1, 1)));
          basalParticle.setPosition(Vec3D(actualXpos, actualYpos, zpos));
          basalParticle.fixParticle();
          particleHandler.copyAndAddObject(basalParticle);
        }
      }
    }

    void createRandomRoughBase() {
      if (params->baseConc <= 0) return;

      Mdouble domainArea = (getXMax() - getXMin()) * (getYMax() - getYMin());
      Mdouble areaPerParticle = M_PI * std::pow(params->baseRadius, 2);
      int numberOfBasalParticles = params->baseConc * domainArea / areaPerParticle;

      for (int i = 0; i < numberOfBasalParticles; i++) {
        double xpos = generator.getRandomNumber(getXMin(), getXMax());
        double ypos = generator.getRandomNumber(getYMin(), getYMax());
        double zpos = 0;

        SphericalParticle basalParticle;
        basalParticle.setSpecies(species);
        basalParticle.setRadius(baseRadius *
                                (1 + params->baseDispersity * generator.getRandomNumber(-1, 1)));
        basalParticle.setPosition(Vec3D(xpos, ypos, zpos));
        basalParticle.fixParticle();
        particleHandler.copyAndAddObject(basalParticle);
      }
    }

    void setupInitialConditions() override {
      /* Start writing to the .muI file. */
      char muICal3D_fn[MAX_STRLEN];
      snprintf(muICal3D_fn, MAX_STRLEN, "%s.muI", getName().c_str());
      muICal3D_f = fopen(muICal3D_fn, "w");
      setbuf(muICal3D_f, nullptr);
      fprintf(muICal3D_f, "time theta n depth mass xmom ke basfx basfy basfz\n");
      fprintf(stderr, "Started writing to .muI file\n");

      setGravity(Vec3D(g * sin(theta), 0, -g * cos(theta)));

      Mdouble maxRadius = particleRadius * (1 + params->dispersity);
      Mdouble maxBaseRadius = baseRadius * (1 + params->baseDispersity);

      SphericalParticle particle;
      particle.setSpecies(species);
      particle.setVelocity(Vec3D(0, 0, 0));

      double zpos = maxRadius + maxBaseRadius;
      while (particleHandler.getVolume() < getTotalVolume()) {
        zpos += 2 * maxRadius;
        for (double xpos = getXMin() + 1 * maxRadius; xpos <= getXMax() - maxRadius; xpos += 2 * maxRadius) {
          for (double ypos = getYMin() + 1 * maxRadius; ypos <= getYMax() - maxRadius; ypos += 2 * maxRadius) {
            particle.setRadius(particleRadius *
                               (1 + params->dispersity * generator.getRandomNumber(-1, 1)));
            particle.setPosition(Vec3D(xpos, ypos, zpos));
            particleHandler.copyAndAddObject(particle);
          }
        }
      }

      logger(INFO, "Particle volume is %, domain volume is %",
             particleHandler.getVolume(), getTotalVolume());
    }

    void actionsOnRestart() override {
      /* Continue writing to the .muI file. */
      char muICal3D_fn[MAX_STRLEN];
      snprintf(muICal3D_fn, MAX_STRLEN, "%s.muI", getName().c_str());
      muICal3D_f = fopen(muICal3D_fn, "a");
      setbuf(muICal3D_f, nullptr);
      logger(INFO, "Continuing writing to .muI file\n");
    }

    Vec3D calculateBasalForce() {
      /* Calculate the forces on the basal particles
       * and the basal wall. */
      Vec3D basalForce;
      for (auto p: particleHandler) {
        if (p->isFixed()) {
          basalForce += p->getForce();
        }
      }
      basalForce += wallHandler.getObject(0)->getForce();

      return basalForce;
    }

    void writeOutputFiles() override {
      Mercury3D::writeOutputFiles();

      if (eneFile.getLastSavedTimeStep() != getNumberOfTimeSteps()) {
        return;
      }

      // Calculate the forces on the basal particles and wall.
      Vec3D basalForce = calculateBasalForce();

      // Write all these details to the .muI file.
      fprintf(
        muICal3D_f,
        "%g %g %d %g %g %g %g %g %g %g\n",
        getTime(),
        theta / DEGREES,
        particleHandler.getNumberOfObjects(),
        2 * getCentreOfMass().Z,
        getTotalMass(),
        getTotalMomentum().X,
        getKineticEnergy(),
        basalForce.X,
        basalForce.Y,
        basalForce.Z
      );
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
    FILE *muICal3D_f;

    double g;
    double theta;  // slope angle in radians
    double particleRadius;
    double baseRadius;
};


int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Example: " << argv[0] << " params.txt 18" << std::endl;
    return 1;
  }

  auto problem = new muICal3D(argv[1], argv[2]);
  argv[2] = argv[0];
  problem->solve(argc - 2, argv + 2);
  delete problem;
  return 0;
}
