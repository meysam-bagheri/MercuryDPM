//Copyright (c) 2013-2024, The MercuryDPM Developers Team. All rights reserved.
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

// A body with changing tensor of inertia featuring a controlled
// angular maneuver

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/ClumpParticle.h"
#include "../ClumpHeaders/ClumpInput.h"
#include "../ClumpHeaders/Mercury3DClump.h"
# include <stdlib.h>
#include<CMakeDefinitions.h>
#include <memory>

Mdouble f_min = -10; Mdouble f_max = 10;

// This function executes the OS command 'cmd' and returns the 'result'
std::string ExecCommand(const char* cmd) {
    std::array<char, 256> buffer;
    std::string result;

    // Prepare a pipe to write a to execute a system command and read the result
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

class ChangingTOIParticle : public Mercury3Dclump
{
public:

    // Helpers functions for matrix manipulations
    MatrixSymmetric3D MtoS( Matrix3D M){ return MatrixSymmetric3D(M.XX, M.XY, M.XZ, M.YY, M.YZ, M.ZZ);}
    Matrix3D StoM( MatrixSymmetric3D M){ return Matrix3D(M.XX, M.XY, M.XZ, M.XY, M.YY, M.YZ, M.XZ, M.YZ, M.ZZ);}
    Matrix3D transpose(Matrix3D M){ return Matrix3D(M.XX, M.YX, M.ZX, M.XY, M.YY, M.ZY, M.XZ, M.YZ, M.ZZ);}

    explicit  ChangingTOIParticle()
    {
        setGravity(Vec3D(0.0, 0.0, 0.0));
        setName("ChangingTOI");
        setXBallsAdditionalArguments("-solidf -v0");
        setXMax(f_max);
        setYMax(f_max);
        setZMax(f_max);
        setXMin(f_min);
        setYMin(f_min);
        setZMin(f_min);
        LoadClumps(data);
        clump_mass = data.mass[clump_index];
    }

    void setClumpDamping(Mdouble damp){ clump_damping = damp;}

    void setClumpIndex(Mdouble index){ clump_index = index;}

    Mdouble getClumpMass(){return clump_mass;}

    void setupInitialConditions() override
    {
        // Generate single clump
        ClumpParticle p0;
        p0.setSpecies(speciesHandler.getObject(0)); // Assign the material type to ClumpParticle
        p0.setClump();
        p0.setRadius(0.5);
        p0.setPosition(Vec3D(0, 0, 0));

        // Pebble particles are used for visualisation purposes only
        // (sphere/arrow glyphs are employed in Paraview to visualize the corresponding directions)
        // Pebble's inertial properties are not employed in calculations
        p0.addPebble(Vec3D(0,0,0),1);
        p0.addPebble(Vec3D(1,0,0),3.5);
        p0.addPebble(Vec3D(-1,0,0),4);
        p0.addPebble(Vec3D(0,1,0),2.5);
        p0.addPebble(Vec3D(0,-1,0),3);
        p0.addPebble(Vec3D(0,0,1),1.5);
        p0.addPebble(Vec3D(0,0,-1),2);

        // Extra direction markers (for Paraview arrow glyphs)
        p0.addPebble(Vec3D(10e-8,10e-8,10e-8),2); // Marker of point 4
        p0.addPebble(Vec3D(10e-8,0,10e-8),2); // marker of point 3
        p0.addPebble(Vec3D(0,10e-8,10e-8),2); // marker of point 2
        p0.addPebble(Vec3D(10e-8,10e-8,0),2); // marker of point 1
        p0.addPebble(Vec3D(10e-8,0,0),2); // marker of point 6
        p0.addPebble(Vec3D(0,10e-8,0),2); // marker of point 7
        p0.addPebble(Vec3D(0,0,10e-8),2); // marker of point 5
        p0.addPebble(Vec3D(0,0,-10e-8),2); // marker of point 8


        p0.setPrincipalDirections(Matrix3D(1,0,0, 0,1,0, 0,0,1));
        
    	MatrixSymmetric3D cInertia;   
        cInertia.XX = inertia_profiles[0][1];
        cInertia.YY = inertia_profiles[0][2];
        cInertia.ZZ = inertia_profiles[0][3];
        cInertia.XY = 0;
        cInertia.XZ = 0;
        cInertia.YZ = 0;
        p0.setInitInertia(cInertia);
                                                                          
        p0.setClumpMass(5);

        p0.setAngularVelocity(baseAngVel * init_orientation);
        particleHandler.copyAndAddObject(p0);
    }


    void actionsAfterTimeStep() override {
        MatrixSymmetric3D cInertia;
        MatrixSymmetric3D inertiaRate;
        MatrixSymmetric3D fInertia;
        fInertia.XX = inertia_profiles[inertia_profiles.size()-1][1];
        fInertia.YY = inertia_profiles[inertia_profiles.size()-1][2];
        fInertia.ZZ = inertia_profiles[inertia_profiles.size()-1][3];
        fInertia.XY = 0;
        fInertia.XZ = 0;
        fInertia.YZ = 0;

    for (std::vector<BaseParticle*>::iterator it= particleHandler.begin(); it!=particleHandler.end(); ++it){
        if ((*it)->isClump()) {
            // Get Principal directions/rotation matrix
            Matrix3D Q = static_cast<ClumpParticle*>(*it)->getRotationMatrix();
            Matrix3D Qt = transpose(Q);


            // Prescribe particle TOI (specified in inertia_profiles array) for every moment of simulation
            int ind = (int) floor(inertia_profiles.size()*(getTime()/progDuration));


            if (ind < inertia_profiles.size()-1) {
                cInertia.XX = inertia_profiles[ind][1];
                cInertia.YY = inertia_profiles[ind][2];
                cInertia.ZZ = inertia_profiles[ind][3];
                cInertia.XY = 0;
                cInertia.XZ = 0;
                cInertia.YZ = 0;

                inertiaRate.XX = inertia_profiles[ind][4];
                inertiaRate.YY = inertia_profiles[ind][5];
                inertiaRate.ZZ = inertia_profiles[ind][6];
                inertiaRate.XY = 0;
                inertiaRate.XZ = 0;
                inertiaRate.YZ = 0;

                MatrixSymmetric3D rotatedCInertia = MtoS(Q * (StoM(cInertia) * Qt));
                static_cast<ClumpParticle *>(*it)->setInitInertia(rotatedCInertia);


                // Add extra torques  - I^{dot} * omega
                Vec3D angVel = (*it)->getAngularVelocity();

                Matrix3D S = Matrix3D(0, -angVel.Z, angVel.Y,   angVel.Z, 0, -angVel.X, -angVel.Y, angVel.X, 0);
                Matrix3D St = transpose(S);
                //MatrixSymmetric3D rotatedInertiaRate = MtoS( Q * StoM(inertiaRate) * Qt +
                //                                                S * StoM(rotatedCInertia) * Qt +
                //                                                Q * StoM(rotatedCInertia) * St);
                MatrixSymmetric3D rotatedInertiaRate = MtoS( Q * StoM(inertiaRate) * Qt);

                Vec3D TorqueDueToRate = rotatedInertiaRate * angVel;
                static_cast<ClumpParticle *>(*it)->setTorque(-TorqueDueToRate);

                // Store angular momentum (for validation purposes)
                angularMomentumLog.push_back((rotatedCInertia * angVel).getLength());
            }

            if (ind == inertia_profiles.size()-1){ // The last step of control program

                // Ensure no torques and static TOI at the final part of the simulation
                static_cast<ClumpParticle*>(*it)->setTorque(Vec3D(0,0,0));
                MatrixSymmetric3D rotatedFInertia = MtoS(Q * (StoM(fInertia) * Qt));
                static_cast<ClumpParticle*>(*it)->setInitInertia(rotatedFInertia);

                // Compute functional (has to be done once in the end of simulation)
                Vec3D angVel = (*it)->getAngularVelocity();

                Vec3D e1 = Vec3D(1,0,0); // Cartesian axes
                Vec3D e2 = Vec3D(0,1,0);
                Vec3D e3 = Vec3D(0,0,1);

                Vec3D n1 = Q * Vec3D(1,0,0); // principal axes orientations
                Vec3D n2 = Q * Vec3D(0,1,0);
                Vec3D n3 = Q * Vec3D(0,0,1);

                Vec3D w = angVel / angVel.getLength(); // angular velocity

                c_theta = acos(Vec3D::dot(n3, w));
                Vec3D n_phi = n1 * Vec3D::dot(n1, w) + n2 * Vec3D::dot(n2, w);
                n_phi = n_phi / n_phi.getLength();
                c_phi = atan2(Vec3D::dot(n2, n_phi), Vec3D::dot(n1, n_phi));

                Vec3D f = final_orientation;

                Mdouble f_theta = acos(Vec3D::dot(e3, f));
                Vec3D e_phi = e1 * Vec3D::dot(e1, f) + e2 * Vec3D::dot(e2, f);
                e_phi = e_phi / e_phi.getLength();
                Mdouble f_phi = atan2(Vec3D::dot(e2, e_phi), Vec3D::dot(e1, e_phi));

                Vec3D k1 = Vec3D(sin(c_theta) * cos(c_phi), sin(c_theta) * sin(c_phi), cos(c_theta));
                Vec3D k2 = Vec3D(sin(f_theta) * cos(f_phi), sin(f_theta) * sin(f_phi), cos(f_theta));
                
                functional = acos(Vec3D::dot(k1, k2)); 
                
                

            }
        }
        }

        Mercury3D::actionsAfterTimeStep();
    }


    Vec3D init_orientation = Vec3D(1,0,0);
    Vec3D final_orientation = Vec3D(1,0,0);
    Double2DVector inertia_profiles;
    Mdouble functional = 0; // Public member to store functional
    Mdouble progDuration = 0; // Time of maneuvering
    Mdouble symDuration = 0; // Total simulation time
    Mdouble baseAngVel = 0;  // Initial angular velocity of the body
    Mdouble c_theta;
    Mdouble c_phi;
    DoubleVector angularMomentumLog;

private:
    int clump_index;
    ClumpData data;
    Mdouble clump_mass;
    Mdouble clump_damping = 10;
    
};


Vec3D LoadInitOrient()
{
    std::cout<<"Read init orientation"<<std::endl;
    // Load init orientation
    std::ifstream infile((getMercuryDPMBuildDir() + "/Drivers/Clump/ChangingTOI/opt/init_orientation.txt").c_str(), std::ios::in | std::ios::binary);
    char lin[256];
    infile.getline(lin, 256, '\n');
    std::string line(lin);
    Mdouble theta = std::stof(line);
    infile.getline(lin, 256, '\n');
    std::string line2(lin);
    Mdouble phi = std::stof(line2);
    infile.close();
    return Vec3D(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

Vec3D LoadFinalOrient()
{
    std::cout<<"Read final orientation"<<std::endl;
    // Load final orientation
    std::ifstream infile((getMercuryDPMBuildDir() + "/Drivers/Clump/ChangingTOI/opt/final_orientation.txt").c_str(), std::ios::in | std::ios::binary);
    char lin[256];
    infile.getline(lin, 256, '\n');
    std::string line(lin);
    Mdouble theta = std::stof(line);
    infile.getline(lin, 256, '\n');
    std::string line2(lin);
    Mdouble phi = std::stof(line2);
    infile.close();
    std::cout<<"LOADED FINAL ORIENTATION: theta = "<<theta<<", phi = "<<phi<<std::endl;
    return Vec3D(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

Double2DVector LoadIProfiles()
{
    std::cout<<"read inertia profile"<<std::endl;
    // Load goal orientation
    std::ifstream infile((getMercuryDPMBuildDir() + "/Drivers/Clump/ChangingTOI/opt/inertia_profiles.txt").c_str(), std::ios::in | std::ios::binary);
    char lin[256];
    Double2DVector i_profiles;
    while (infile.getline(lin, 256, '\n')){
    		StringVector val;
    		std::string line(lin);
    		line+=" ";
    		std::string buffer = "";
    		for ( int j = 0; j < line.size(); j++){
    			if (line[j] != ' ') {buffer += line[j];}
        		else {val.push_back(buffer); buffer = "";}
        	}
        	DoubleVector row;
        	for (int i = 0; i < 7; i++){ row.push_back(std::stof(val[i])); }
        	i_profiles.push_back(row);	
    	}
    	infile.close();
    return i_profiles;
}

int main(int argc, char* argv[])
{
    ChangingTOIParticle problem;
    problem.init_orientation = LoadInitOrient();
    problem.final_orientation = LoadFinalOrient();
    problem.inertia_profiles = LoadIProfiles();

    // Get parameters passed through the command line
    std::string a;
    problem.progDuration = stod(helpers::readFromCommandLine(argc, argv, "-p1", a));
    problem.symDuration = stod(helpers::readFromCommandLine(argc, argv, "-p2", a));
    problem.baseAngVel = stod(helpers::readFromCommandLine(argc, argv, "-p3", a));

    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    species->setDensity(1.0); // sets the species type-0 density
    species->setDissipation(0.0);
    species->setStiffness(1e6);
    const Mdouble collisionTime = species->getCollisionTime(problem.getClumpMass());
    problem.setClumpDamping(0);
    problem.setTimeStep(collisionTime/50 );
    std::cout<<"TIMESTEP: "<< collisionTime/50  << std::endl;

    // Quick demonstration
    problem.setSaveCount(500);
    problem.setTimeMax(problem.symDuration);
    problem.removeOldFiles();
    problem.solve();
    // Paraview data
    ExecCommand("rm -rf paraview_ChangingTOI/");
    ExecCommand("mkdir paraview_ChangingTOI/");
    ExecCommand("../../../Tools/data2pvd ChangingTOI.data paraview_ChangingTOI/ChangingTOI");
    std::string command;
    command = "python " + getMercuryDPMSourceDir() + "/Tools/MClump/PlotEnergies.py " + getMercuryDPMBuildDir() + "/Drivers/Clump/ChangingTOI/ " + "ChangingTOI";
    ExecCommand(command.c_str());

    // Return the functional via the text file
    std::ofstream funct;  funct.open ("opt/functional.txt");
    funct << problem.functional;
    funct.close();

    // Return current angles
    std::ofstream angl;  angl.open ("opt/angles.txt");
    angl << problem.c_theta << std::endl;
    angl << problem.c_phi <<std::endl;
    angl.close();

    // Return the log of the angular momentum
    std::ofstream mom;  mom.open ("opt/momentum.txt");
    for (int i = 0; i<problem.angularMomentumLog.size(); i+=500){mom <<problem.angularMomentumLog[i]<<std::endl;}
    mom.close();
    return 0;
}
