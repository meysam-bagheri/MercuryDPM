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

#include <cmath>

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"

/*!
 * Class to test if two particles that have the same initial overlap with a wall have consistent movement. Note that
 * not "correct" movement is checked, purely consistency.
 * Could be broken if changes are made in time-integrator or if not all fields are initialised correctly
 * (e.g. tangential spring), currently only test LinearViscoelasticFrictionSpecies.
 */
class ConsistencyInitialParticleWallOverlapUnitTest : public Mercury3D
{
public:
    void setupInitialConditions() override {
        LinearViscoelasticFrictionSpecies* s = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        s->setDensity(1000.);
        double mass = s->getMassFromRadius(getRadius());
        s->setCollisionTimeAndRestitutionCoefficient(0.005, 0.8, mass);

        s->setSlidingDissipation(s->getDissipation()*2./7.);
        s->setSlidingStiffness(s->getStiffness()*2./7.);
        s->setSlidingFrictionCoefficient(0.4);
        s->setRollingStiffness(s->getStiffness()*2.0/7.0);
        s->setRollingFrictionCoefficient(0.2);
        s->setRollingDissipation(s->getDissipation()*2./7.);
        s->setTorsionStiffness(s->getStiffness()*2.0/7.0);
        s->setTorsionFrictionCoefficient(0.1);
        s->setTorsionDissipation(s->getDissipation()*2./7.);
        speciesHandler.copyAndAddObject(s);

        SphericalParticle p0;
        p0.setSpecies(s);
        p0.setRadius(getRadius());
        p0.setPosition(getInitPos0());
        particleHandler.copyAndAddObject(p0);
        p0.setPosition(getInitPos1());
        particleHandler.copyAndAddObject(p0);

        InfiniteWall w;
        w.setSpecies(s);
        w.set({0., 0., -1.0},{0., 0., 0.0});
        wallHandler.copyAndAddObject(w);
    }

    void setInitOverlap(const double& d) {initOverlap_ = d;}
    double getInitOverlap() { return initOverlap_;}
    void setRadius(const double& r) {pRad_ = r;}
    double getRadius() { return pRad_;}
    void setInitPos0(const Vec3D& x) {initPos0_ = x;}
    Vec3D getInitPos0() { return initPos0_;}
    void setInitPos1(const Vec3D& x) {initPos1_ = x;}
    Vec3D getInitPos1() { return initPos1_;}

private:
    double initOverlap_;
    double pRad_;
    Vec3D initPos0_;
    Vec3D initPos1_;
};

int main(int argc, char *argv[])
{
    // Problem setup
    ConsistencyInitialParticleWallOverlapUnitTest problem;
    problem.setXMax(0.4);
    problem.setYMax(0.6);
    problem.setZMax(0.1);
    problem.setXMin(0.0);
    problem.setYMin(0.2);
    problem.setZMin(0.0);

    problem.setGravity(Vec3D{0.,0.,-9.81});

    problem.setRadius(1.e-2);
    problem.setInitOverlap(1.e-3);
    problem.setInitPos0({0., 0.4, problem.getRadius() - problem.getInitOverlap()});
    problem.setInitPos1({0.2, 0.4, problem.getRadius() - problem.getInitOverlap()});

    problem.setTimeMax(1.e-1);
    problem.setTimeStep(0.005/50.);

    problem.setName("ConsistencyInitialParticleWallOverlapUnitTest");
    // Uncomment to check vtu files
    //problem.setParticlesWriteVTK(true);
    //problem.setWallsWriteVTK(true);
    problem.dataFile.setFileType(FileType::NO_FILE);
    problem.restartFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);

    problem.solve();

    Vec3D initDiff = problem.getInitPos1() - problem.getInitPos0();
    Vec3D currDiff = problem.particleHandler.getObject(1)->getPosition() - problem.particleHandler.getObject(0)->getPosition();

    // Check if particle 0 and 1 only moved in z-dir
    helpers::check(problem.getInitPos0().X, problem.particleHandler.getObject(0)->getPosition().X, 1.e-15, "Particle 0 x-movement");
    helpers::check(problem.getInitPos0().Y, problem.particleHandler.getObject(0)->getPosition().Y, 1.e-15, "Particle 0 y-movement");

    helpers::check(problem.getInitPos1().X, problem.particleHandler.getObject(1)->getPosition().X, 1.e-15, "Particle 1 x-movement");
    helpers::check(problem.getInitPos1().Y, problem.particleHandler.getObject(1)->getPosition().Y, 1.e-15, "Particle 1 y-movement");

    helpers::check((currDiff.X-initDiff.X),0.,1.e-15,"Difference in x-movement");
    helpers::check((currDiff.Y-initDiff.Y),0.,1.e-15,"Difference in y-movement");
    helpers::check((currDiff.Z-initDiff.Z),0.,1.e-15,"Difference in z-movement");

    return 0;
}
