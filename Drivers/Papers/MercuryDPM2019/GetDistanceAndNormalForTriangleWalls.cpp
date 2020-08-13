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

#include <Walls/Coil.h>
#include <Mercury3D.h>
#include <Species/LinearViscoelasticReversibleAdhesiveSpecies.h>
#include <Walls/TriangleWall.h>

/**
 * \brief Tests the contact detection between particles and a set of TriangleWalls.
 * \detail In particular, distinguishing face, edge and vertex contacts is tricky.
 * The most difficult case is when a face is less or equal in size to a particle, so this is tested here.
 **/
class GetDistanceAndNormalForTriangleWalls : public Mercury3D
{
public:
    void setupInitialConditions() override
    {
        setName("GetDistanceAndNormalForTriangleWalls");
        setFileType(FileType::NO_FILE);
        dataFile.setFileType(FileType::ONE_FILE);
        setTimeStep(1e-6);
        setTimeMax(getTimeStep());
        setMin(-1 * Vec3D(1, 1, 1));
        setMax(1.001 * Vec3D(1, 1, 1));
        setXBallsAdditionalArguments("-w 1100 -s 0.9");

        //use an adhesive species to visualise negative overlaps
        LinearViscoelasticReversibleAdhesiveSpecies species;
        species.setDensity(1);
        species.setStiffness(1);
        species.setAdhesionStiffness(.1);
        species.setAdhesionForceMax(1);
        auto s = speciesHandler.copyAndAddObject(species);

        TriangleWall wall;
        wall.setSpecies(s);
        wall.setGroupId(1);
        double x[3]{-.5, 0, .5};
        double y = .5;
        double z[3]{-.5, 0, .5};
        Vec3D P[5]{{x[0], -y, z[0]},
                   {x[0], -y, z[1]},
                   {x[1], -y, z[1]},
                   {x[2], -y, z[1]},
                   {x[2], -y, z[2]}};
        Vec3D Q[5]{{x[0], +y, z[0]},
                   {x[0], +y, z[1]},
                   {x[1], +y, z[1]},
                   {x[2], +y, z[1]},
                   {x[2], +y, z[2]}};

        for (int i = 0; i < 4; i++) {
            wall.setVertices(P[i], P[i+1], Q[i]);
            wallHandler.copyAndAddObject(wall);
            wall.setVertices(Q[i], P[i+1], Q[i+1]);
            wallHandler.copyAndAddObject(wall);
        }
        auto w = wallHandler.getObject(3);
        logger(INFO,"wall %",*w);

        SphericalParticle particle;
        particle.setSpecies(s);
        particle.setRadius(0.5);
        auto p = particleHandler.copyAndAddObject(particle);

        Vec3D pos;
        Vec3D normal;
        Mdouble distance;
        Mdouble h = 0.01;
        //pos.Y=0.5;
        std::ofstream file("GetDistanceAndNormalForTriangleWalls.out");
        file << "x\tz\td\tnx\tnz\n";
        for (pos.X = getXMin(); pos.X <= getXMax(); pos.X += h) {
            for (pos.Z = getZMin(); pos.Z <= getZMax(); pos.Z += h) {
                p->setPosition(pos);
                double overlaps = 0;
                Vec3D force {0,0,0};
                unsigned count = 0;
                double radius = 0.2;
                for (auto wall : wallHandler) {
                    wall->getDistanceAndNormal(*p, distance, normal);
                    double overlap = std::max(radius-distance,0.0);
                    overlaps += overlap;
                    force += -normal*overlap;
                    if (overlap>0) count++;
                }
                w->getDistanceAndNormal(*p, distance, normal);
                file << pos.X
                     << '\t' << pos.Z
                        << '\t' << distance
                        << '\t' << normal.X
                        << '\t' << normal.Z
                        << '\t' << overlaps
                        << '\t' << force.X
                        << '\t' << force.Z
                        << '\t' << count
                     << '\n';
            }
        }
        file.close();

        helpers::writeToFile("GetDistanceAndNormalForTriangleWalls.m",
                "clear variables\n"
                "raw=importdata('GetDistanceAndNormalForTriangleWalls.out','\\t',1)\n"
                "n(1)=length(unique(raw.data(:,1)));\n"
                "n(2)=length(unique(raw.data(:,2)));\n"
                "x=reshape(raw.data(:,1),n);\n"
                "z=reshape(raw.data(:,2),n);\n"
                "d=reshape(raw.data(:,3),n);\n"
                "nx=reshape(raw.data(:,4),n);\n"
                "nz=reshape(raw.data(:,5),n);\n"
                "o=reshape(raw.data(:,6),n);\n"
                "fx=reshape(raw.data(:,7),n);\n"
                "fz=reshape(raw.data(:,8),n);\n"
                "count=reshape(raw.data(:,9),n);\n"
                "\n"
                "% figure(1)\n"
                "% contourf(x,z,d,20,'EdgeColor','none')\n"
                "% hold on\n"
                "% contour(x,z,d,[0 0],'EdgeColor','k')\n"
                "% plotGeometry\n"
                "% hold off\n"
                "% xlabel('x')\n"
                "% ylabel('z')\n"
                "% h=colorbar\n"
                "% %xlabel(h, 'distance')\n"
                "% title('distance')\n"
                "% saveas(gcf,'distanceTriangleWalls.png')\n"
                "% \n"
                "% figure(2)\n"
                "% m=6;\n"
                "% quiver(x(1:m:end,1:m:end),z(1:m:end,1:m:end),nx(1:m:end,1:m:end),nz(1:m:end,1:m:end),0.5)\n"
                "% hold on\n"
                "% contour(x,z,d,[0 0],'EdgeColor','k')\n"
                "% plotGeometry\n"
                "% hold off\n"
                "% xlabel('x')\n"
                "% ylabel('z')\n"
                "% title('normal')\n"
                "% saveas(gcf,'normalTriangleWalls.png')\n"
                "\n"
                "figure(3)\n"
                "contourf(x,z,o,20,'EdgeColor','none')\n"
                "hold on\n"
                "contour(x,z,o,[0 0],'EdgeColor','k')\n"
                "plotGeometry\n"
                "hold off\n"
                "xlabel('x')\n"
                "ylabel('z')\n"
                "h=colorbar\n"
                "%xlabel(h, 'distance')\n"
                "title('distance')\n"
                "saveas(gcf,'distanceMinTriangleWalls.png')\n"
                "\n"
                "figure(4)\n"
                "m=6;\n"
                "quiver(x(1:m:end,1:m:end),z(1:m:end,1:m:end),fx(1:m:end,1:m:end),fz(1:m:end,1:m:end),0.5)\n"
                "hold on\n"
                "contour(x,z,d,[0 0],'EdgeColor','k')\n"
                "plotGeometry\n"
                "hold off\n"
                "xlabel('x')\n"
                "ylabel('z')\n"
                "title('normal')\n"
                "saveas(gcf,'normalMinTriangleWalls.png')\n"
                "%%\n"
                "figure(5)\n"
                "m=6;\n"
                "contourf(x,z,count,20)\n"
                "hold on\n"
                "plotGeometry\n"
                "hold off\n"
                "xlabel('x')\n"
                "ylabel('z')\n"
                "colormap('colorcube')\n"
                "colorbar\n"
                "title('count')\n"
                "saveas(gcf,'countTriangleWalls.png')\n"
                "\n"
                "function plotGeometry\n"
                "plot(.5*[-1 -1 0 1 1],.5*[-1 0 0 0 1],'k')\n"
                "end\n"
                );
        logger(INFO,"Written file");
        setWallsWriteVTK(true);
        writeVTKFiles();
        exit(0);
    }
};

int main(int argc, char* argv[])
{
    GetDistanceAndNormalForTriangleWalls().solve();
    return 0;
}
