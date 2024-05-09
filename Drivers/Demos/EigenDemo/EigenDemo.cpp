//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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
#include <random>
#include "Mercury3D.h"
#include <Boundaries/PeriodicBoundary.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include "Particles/SphericalParticle.h"
#include "DPMBase.h"
#include <vector>
#include <string>
#include <algorithm>
#include <Eigen/Sparse>


class CFDDEMCoupleTest : public Mercury3D {
private:
    double D50_=0.00025, kappa=0.4, rho_f=1.225, nu= 1.46e-5;
    double theta=0.05;//Shields=0.05
    double dz, u_star;
    int N_grid;//number of grid points of flow field
    std::vector<double> z, u_f, l_f, u, l, tau_f_vec, F, phi_p, T;

public:
    typedef Eigen::SparseMatrix<double> mat;//define a sparse matrix type
    typedef Eigen::Triplet<double> Tri;//define a triplet type
    std::vector<Tri> tripletList;
    std::ofstream out;

    void setupInitialConditions() override {

        readRestartFile("GenerateSmallK0Bed");
        setRestarted(false);
        setName("EigenDemo");
        setFileType(FileType::ONE_FILE);
        setZMax(0.2);

        std::string Name=getName();
        logger(INFO, Name);

        auto originalSpecies = dynamic_cast<LinearViscoelasticFrictionSpecies *>(speciesHandler.getObject(
                0));

        double massMin = originalSpecies->getMassFromRadius(particleHandler.getSmallestParticle()->getRadius());
        const Mdouble collisionTime = originalSpecies->getCollisionTime(massMin);
        setTimeMax(0.15);
        setTimeStep(collisionTime/50);
        logger(INFO,"timestep=%",getTimeStep());
        setSaveCount(1000);

        //add impact particles
        SphericalParticle particle;
        particle.setSpecies(speciesHandler.getObject(0));
        particle.setRadius(particleHandler.getLargestParticle()->getRadius());
        particle.setPosition(Vec3D(0.16*getXMax(),0.5*getYMax(), 7*D50_));
        particle.setVelocity(Vec3D(0.2*5,0,-0.2*5));
        particleHandler.copyAndAddObject(particle);

        //fluid field information
        u_star = sqrt(theta * (originalSpecies->getDensity() - rho_f) * 9.81 * D50_ / rho_f);
        N_grid = (getZMax() - getZMin()) / D50_*2 + 1;
        z.resize(N_grid);
        u.resize(N_grid);
        l.resize(N_grid - 1);
        double x;
        u[0] = 0;
        z[0] = 0;
        std::vector<double> data;
        /*input the flow velocity and mixing length from the file of AdaptingTest*/
        std::ifstream filein("FlowShields005.txt");
        if (!filein.is_open())
            logger(ERROR, "Error opening file");
        while (filein >> x) {
            data.push_back(x);
        }

        for (int n = 0; n < data.size(); n++) {
            if (n % 3 == 0) {
                z[n / 3 + 1] = data[n];
            } else if (n % 3 == 1) {
                u[n / 3 + 1] = data[n];
            } else {
                l[n / 3] = data[n];
            }
        }
        dz = z[1] - z[0];
    }


//member functions for updating l
    double itterateL(double l_down, double u) {
        double l_up = l_down;//initial guess
        int itter_l = 0;
        double res_l = 1;
        double dresdl;
        while (fabs(res_l) > 1e-8 && itter_l < 600) {
            itter_l = itter_l + 1;
            res_l = (l_down - l_up) / dz +
                    kappa * (1 - exp(-sqrt(1.0 / 7.0 * (fabs(u) * 0.5 * (l_up + l_down)) / nu)));//(20)rhs - lhs
            dresdl = -1 / dz - kappa * exp(-sqrt(1.0 / 7.0 * (fabs(u) * 0.5 * (l_up + l_down)) / nu)) * 0.5 *
                               pow((1.0 / 7.0 * (fabs(u) * 0.5 * (l_up + l_down)) / nu), -0.5) * 1.0 / 7.0 *
                               fabs(u) * 0.5 / nu;
            l_up = l_up - 0.8 * res_l / dresdl;
        }
        return l_up;
    }

    std::pair<std::vector<double>, std::vector<double>> CalcPhiPAndFz() {
        std::vector<double> phi_p_vec(N_grid, 0);
        std::vector<double> F_z_vec(N_grid, 0);
        for (BaseParticle *p: particleHandler) {
            int z_p = p->getPosition().Z / dz;//find the lower layer id near the particle
            if (p->getPosition().Z+p->getRadius() <= getZMax()-dz*0.5) {
                double u_rel =
                        u[z_p] + (u[z_p + 1] - u[z_p]) * (p->getPosition().Z - z[z_p]) / dz -
                        p->getVelocity().X;
                double f_drag = constants::pi / 8.0 * 1.225 * pow(p->getRadius() * 2.0, 2) *
                                (0.5*fabs(u_rel) + 24.0*nu/(p->getRadius()*2) + 2.0*sqrt(12.0*nu/(p->getRadius()*2))*sqrt(fabs(u_rel))) * u_rel;

                double A = (getXMax() - getXMin()) * (getYMax() - getYMin());//p->getRadius() * 2;//d vary or not?
                //lower part of a particle
                int iz_low = (p->getPosition().Z - p->getRadius()) * 2.0 / dz;
                int iz_l;
                if (iz_low % 2 == 0) {
                    iz_l = iz_low / 2.0;
                } else {
                    iz_l = (iz_low + 1) / 2.0;
                }

                double H_0, H_1, H_2, H_3, V_part0, V_part1, V_part2, V_part3;//the volume of a part of the particle in this layer
                if (p->getPosition().Z + p->getRadius() > z[iz_l] + dz / 2.0) {
                    //V_part = p->getVolume();
                    H_0 = z[iz_l] + dz / 2.0 - (p->getPosition().Z - p->getRadius());
                    V_part0 = constants::pi * pow(H_0, 2) * (p->getRadius() - H_0 / 3.0);
                    phi_p_vec[iz_l] += 1.0 / (A * dz) * V_part0;
                    F_z_vec[iz_l] += 1.0 / (A*dz) * f_drag * V_part0/p->getVolume();

                    if(p->getPosition().Z + p->getRadius() > z[iz_l+2]+dz/2.0 ){
                        H_1 = z[iz_l+1] + dz / 2.0 - p->getPosition().Z + p->getRadius();
                        V_part1 = constants::pi * pow(H_1, 2) * (p->getRadius() - H_1 / 3.0) - V_part0;
                        phi_p_vec[iz_l+1] += 1.0 / (A * dz) * V_part1;
                        F_z_vec[iz_l+1] += 1.0 / (A*dz) * f_drag * V_part1/p->getVolume();

                        H_3 = p->getPosition().Z + p->getRadius() - (z[iz_l+2]+dz/2.0);
                        V_part3 = constants::pi * pow(H_3, 2) * (p->getRadius() - H_3 / 3.0);
                        phi_p_vec[iz_l+3] += 1.0 / (A * dz) * V_part3;
                        F_z_vec[iz_l+3] += 1.0 / (A*dz) * f_drag * V_part3/p->getVolume();

                        H_2 = p->getPosition().Z + p->getRadius() - (z[iz_l+1]+dz/2.0);
                        V_part2 = constants::pi * pow(H_2, 2) * (p->getRadius() - H_2 / 3.0) - V_part3;
                        phi_p_vec[iz_l+2] += 1.0 / (A * dz) * V_part2;
                        F_z_vec[iz_l+2] += 1.0 / (A*dz) * f_drag * V_part2/p->getVolume();
                    }
                    else if(p->getPosition().Z + p->getRadius() > z[iz_l+1]+dz/2.0 ){
                        H_1 = z[iz_l+1] + dz / 2.0 - p->getPosition().Z + p->getRadius();
                        V_part1 = constants::pi * pow(H_1, 2) * (p->getRadius() - H_1 / 3.0) - V_part0;
                        phi_p_vec[iz_l+1] += 1.0 / (A * dz) * V_part1;
                        F_z_vec[iz_l+1] += 1.0 / (A*dz) * f_drag * V_part1/p->getVolume();

                        H_2 = p->getPosition().Z + p->getRadius() - (z[iz_l+1]+dz/2.0);
                        V_part2 = constants::pi * pow(H_2, 2) * (p->getRadius() - H_2 / 3.0);
                        phi_p_vec[iz_l+2] += 1.0 / (A * dz) * V_part2;
                        F_z_vec[iz_l+2] += 1.0 / (A*dz) * f_drag * V_part2/p->getVolume();
                    }
                    else{
                        H_1 = p->getPosition().Z + p->getRadius() - (z[iz_l]+dz/2.0);
                        V_part1 = constants::pi * pow(H_1, 2) * (p->getRadius() - H_1 / 3.0);
                        phi_p_vec[iz_l+1] += 1.0 / (A * dz) * V_part1;
                        F_z_vec[iz_l+1] += 1.0 / (A*dz) * f_drag * V_part1/p->getVolume();
                    }

                }
                else {
                    V_part0 = p->getVolume();
                    phi_p_vec[iz_l] += 1.0 / (A * dz) * V_part0;
                    F_z_vec[iz_l] += 1.0 / (A*dz) * f_drag * V_part0/p->getVolume();
                }
            }
        }
        return std::make_pair(phi_p_vec, F_z_vec);
    }



    void computeExternalForces(BaseParticle *p) override {
        DPMBase::computeExternalForces(p);
        int z_p = p->getPosition().Z / dz;//find the lower layer id near the particle
        double u_rel;
        if (z_p < u.size() - 1) {
            u_rel = u[z_p] + (u[z_p + 1] - u[z_p]) * (p->getPosition().Z - z[z_p]) / dz - p->getVelocity().X;
        } else {
            u_rel = u.back() - p->getVelocity().X;
        }
        double f_drag = constants::pi / 8.0 * 1.225 * pow(p->getRadius() * 2.0, 2) *
                        (0.5 * fabs(u_rel) + 24.0 * nu / (p->getRadius() * 2) +
                         2.0 * sqrt(12.0 * nu / (p->getRadius() * 2)) * sqrt(fabs(u_rel))) * u_rel;
        p->addForce(Vec3D(f_drag, 0, 0));
    }


    void actionsAfterTimeStep() override {
        /**clear in each time step**/
        F.resize(N_grid, 0);
        phi_p.resize(N_grid, 0);
        T.resize(N_grid, 0);
        tripletList.clear();

        phi_p = CalcPhiPAndFz().first;
        F = CalcPhiPAndFz().second;

        /**construct coefficient sparse matrix A for A*u=B**/
        /**Firstly construct T; T_i = dt*(nu+l[i]^2*fabs((u[i+1]-u[i])/dz))/dz^2**/
        for (int i = 0; i < N_grid - 1; i++) {
            T[i] = getTimeStep() * (nu + pow(l[i], 2) * fabs(u[i + 1] - u[i]) / dz) / dz / dz;
        }

        /**diagonal elements**/
        tripletList.push_back(Tri(0, 0, 1));
        for (int i = 1; i < N_grid - 1; i++) {
            tripletList.push_back(Tri(i, i, 1 + T[i] + T[i - 1]));
        }
        tripletList.push_back(Tri(N_grid - 1, N_grid - 1, 1 + T[N_grid - 2]));//last coefficient
        /**right-diagonal elements**/
        for (int j = 1; j < N_grid - 1; j++) {
            tripletList.push_back(Tri(j, j + 1, -T[j]));
        }
        /**left-diagonal elements**/
        for (int k = 0; k < N_grid - 1; k++) {
            tripletList.push_back(Tri(k + 1, k, -T[k]));
        }

        mat A(N_grid, N_grid);
        A.setFromTriplets(tripletList.begin(), tripletList.end());

        /**build vector B**/
        Eigen::VectorXd U(N_grid), B(N_grid);
        B(0) = 0;
        for (int m = 1; m < N_grid - 1; m++) {
            B(m) = u[m] - getTimeStep() * F[m] / (1 - phi_p[m]) / rho_f;
        }
        B(N_grid - 1) =
                u[N_grid - 1] + getTimeStep() * (pow(u_star, 2) / dz - F[N_grid - 1] / (1 - phi_p[N_grid - 1]) / rho_f);

        /**define a solver and solve the matrix equation**/
        Eigen::SimplicialLLT<mat> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        U = solver.solve(B);

        /**transfer the value from U to u**/
        for (int i = 0; i < N_grid; i++) {
            u[i] = U(i);
        }

        l[0] = kappa * dz / 2.0 * (1 - exp(-1 / 26.0 * dz / 2.0 * u_star / nu));//get the first l by (21)???
        for (int i_l = 1; i_l < N_grid - 1; i_l++) {
            l[i_l] = itterateL(l[i_l - 1], u[i_l]);
        }

    }


    void printTime() const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime() << std::endl;
        for (int io = 0; io < N_grid-1; io++) {
            std::cout << z[io+1] << '\t'  << u[io+1] << '\t' << l[io] << '\t' << std::endl;
        }

    }


};

int main(int argc, char *argv[]){

    CFDDEMCoupleTest sp0;
    sp0.solve(argc,argv);

    return 0;
}





