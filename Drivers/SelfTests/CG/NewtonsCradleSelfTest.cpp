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
#include "StatisticsVector.h"
#include "Walls/InfiniteWall.h"
#include <array>
#include <Species/LinearViscoelasticSpecies.h>
#include <CG/TimeAveragedCG.h>

using constants::pi;

class NewtonsCradleSelfTest : public DPMBase {
	
public:
	
	void setupInitialConditions() override {
		// make sure the number of particles is right
		logger(INFO,"Creating a cubic packing of %^3 particles",N);
		//set_N(mathsFunc::cubic(N));
		Mdouble Radius = 0.5;
		
		//set Particles' position, radius, velocity and bounding box
		setXMin(0);
		setYMin(0);
		setZMin(0);
		
		setXMax(2.*Radius);
		setYMax(2.*Radius);
		setZMax(N*2.*Radius);
		
		SphericalParticle P0;
        P0.setSpecies(speciesHandler.getObject(0));
        P0.setPosition(Radius * Vec3D(1.,1.,0.));
        for (int k = 0; k < N; k++)
        {
            P0.setRadius(Radius * (1. - .5 * ((Mdouble)k / N)));
            P0.setPosition(P0.getPosition()+Vec3D(0,0,P0.getRadius()));
            particleHandler.copyAndAddObject(P0);
            P0.setPosition(P0.getPosition()+Vec3D(0,0,P0.getRadius()));
        }
        setZMax(P0.getPosition().Z);
		//set walls
		InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(0));
		w.set(Vec3D(0.0,0.0,-1.0),Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w);
	}
	
	int N;	
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	logger(INFO,"Simulating a stack of particles under gravity on a wall until it is relaxed");
	NewtonsCradleSelfTest problem;
	problem.setName("NewtonsCradleSelfTest");
    auto species=problem.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    problem.N=5;	//set the number of particles
    problem.removeOldFiles();
    problem.setSystemDimensions(3);
	problem.setParticleDimensions(3);
	species->setDensity(6./pi);
	problem.setGravity(Vec3D(0.,0.,-1.));
	species->setCollisionTimeAndRestitutionCoefficient(.01,.1,1.);
    problem.setTimeStep(0.0001);
	problem.setTimeMax(8.0);
	problem.setSaveCount(1000);
	problem.solve();

	Mdouble w=0.46/3.;

	logger(INFO,"\nTesting CGHandler::restartAndEvaluateDataFiles");

    Mercury3D cg;
    auto c0 = cg.cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::Z,CGFunctions::Gauss>());
    c0->setNZ(400);
    c0->setWidth(w);
    c0->setTimeMin(8+1e-12);
    c0->statFile.setName(problem.getName()+".fromData.stat");
    cg.cgHandler.restartAndEvaluateDataFiles("NewtonsCradleSelfTest");

    logger(INFO,"\nTesting CGHandler::restartAndEvaluateRestartFiles");

    c0->statFile.setName(problem.getName()+".fromRestart.stat");
    cg.cgHandler.restartAndEvaluateRestartFiles("NewtonsCradleSelfTest");

    ///\todo TW the output of the old and new stat doesn't compare perfectly, check.
    helpers::writeToFile("NewtonsCradleSelfTest.gnu",
                         "p 'NewtonsCradleSelfTest.0.stat' u 2:22 w lp, 'NewtonsCradleSelfTest.stat' u 3:41 w l");

	StatisticsVector<Z> stats1("NewtonsCradleSelfTest");
    stats1.setNZ(400);
    stats1.setCGWidth(w);
	stats1.setSuperExact(false);
	stats1.setCGShape("Gaussian");
	stats1.setCGTimeMin(problem.getTimeMax()*1.00000999999);
	stats1.setTimeMaxStat(1e20);
	//stats1.statFile.setName("NewtonsCradleSelfTest.old.stat");
    logger(INFO,"Creating %",stats1.statFile.getName());
    stats1.statistics_from_fstat_and_data();

	logger(INFO,"Test StatisticsVector for different CG functions (Gauss, Heaviside, Lucy)");
	int i=0;
	std::array<std::string,3> str = {"Gaussian","HeavisideSphere","Lucy"};

	for (int i=0; i<3; i++)	{
		double w=0.125;
		if (i==1) w*=2;
		if (i==2) w*=2*sqrt(3);

		StatisticsVector<Z> stats0("NewtonsCradleSelfTest");
		stats0.setZMinStat(-0.5);
	    stats0.statFile.setName(stats0.getName() + "_" + str[i] + "Z.stat");
		stats0.set_h(0.02);
        stats0.verbose();
		stats0.setCGWidth(w);
		stats0.setSuperExact(true);
		stats0.setCGShape(str[i].c_str());
		stats0.setCGTimeMin(problem.getTimeMax()*.999999);
		stats0.setTimeMaxStat(1e20);
        logger(INFO,"Creating %",stats0.statFile.getName());
		stats0.statistics_from_fstat_and_data();

		StatisticsVector<XYZ> stats2("NewtonsCradleSelfTest");
		stats2.setZMinStat(-0.5);
	    stats2.statFile.setName(stats2.getName() + "_" + str[i] + "XYZ.stat");
		stats2.set_h(0.02);
		stats2.setNY(1);
		stats2.setCGWidth(w);
		stats2.setSuperExact(true);
		stats2.setCGShape(str[i].c_str());
		stats2.setCGTimeMin(problem.getTimeMax()*.999999);
		stats2.setTimeMaxStat(1e20);
        logger(INFO,"Creating %",stats2.statFile.getName());
		stats2.statistics_from_fstat_and_data();

	}

	return(0);
}

