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

#include "VTKWriter/WallVTKWriter.h"
#include "DPMBase.h"


void WallVTKWriter::getVTKData(VTKContainer& vtk) const
{
    // set capacity of points and cells based on the previous time step
    // the initial values are based on the minimum, which is one triangle per wall
    static unsigned int capacityPoints = 3*handler_.getSize();
    static unsigned int capacityTriangleStrips = handler_.getSize();
    vtk.triangleStrips.reserve(capacityPoints);
    vtk.points.reserve(capacityTriangleStrips);

    //add all wall data to the point and cell arrays
    for (const auto& w: handler_)
    {
        w->renderWall(vtk);
    }

//    logger(INFO, "size (capacity) of points: % (%), cells: % (%)", vtk.points.size(), capacityPoints, vtk.triangleStrips.size(), capacityTriangleStrips);

    //store from previous time step
    capacityPoints = vtk.points.size();
    capacityTriangleStrips = vtk.triangleStrips.size();
}

void WallVTKWriter::writeVTK() const
{
    if (PROCESSOR_ID!=0) return;
    std::fstream file = makeVTKFileWithHeader();
    
    VTKContainer vtk;
    getVTKData(vtk);
    file << "<Piece NumberOfPoints=\"" << vtk.points.size()
         << "\" NumberOfCells=\"" << vtk.triangleStrips.size()
         << "\">\n"
         << "<Points>\n";
    writeVTKPoints(file, vtk);
    file << "</Points>\n";
    // this cannot be done right now since we can't math Points/Cells and Walls
    //    file << "<PointData  Vectors=\"vector\">\n";
    //    write(file,"indSpecies",[](BaseWall* w){return w->getIndSpecies();});
    //    file << "</PointData>\n";
    file << "<Cells>\n";
    writeVTKCells(file, vtk);
    file << "</Cells>\n";
    file << "<CellData>\n";
    writeVTKCellData(file, vtk);
    file << "</CellData>\n"
         << "</Piece>\n"
         << "</UnstructuredGrid>\n"
         << "</VTKFile>\n";
    file.close();
}

void WallVTKWriter::write(std::fstream& file, std::string name, std::function<double(BaseWall*)> f) const
{
    file << "  <DataArray type=\"Float32\" Name=\""+ name + "\" format=\"ascii\">\n";
    for (const auto& p: handler_)
    {
        file << '\t' << f(p) << '\n';
    }
    file << "  </DataArray>\n";
}

void WallVTKWriter::writeVTKPoints(std::fstream& file, VTKContainer& vtk) const
{
    file << "  <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const Vec3D& p : vtk.points)
    {
        file << '\t' << p << '\n';
    }
    file << "  </DataArray>\n";
}

void WallVTKWriter::writeVTKCells(std::fstream& file, VTKContainer& vtk) const
{
    file << "  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const std::vector<double>& c : vtk.triangleStrips)
    {
        file << '\t';
        for (const double& i : c) file << i << ' ';
        file << '\n';
    }
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    unsigned count = 0;
    for (const std::vector<double>& c : vtk.triangleStrips)
    {
        count += c.size();
        file << '\t' << count << '\n';
    }
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"UInt8\"  Name=\"types\" format=\"ascii\">\n";
    for (const std::vector<double>& c : vtk.triangleStrips)
    {
        if (c.front() == c.back())
        {
            //polygon
            file << "\t7\n";
        }
        else
        {
            //triangle strips
            file << "\t6\n";
        }
    }
    file << "  </DataArray>\n";
}

void WallVTKWriter::writeVTKCellData(std::fstream& file, VTKContainer& vtk) const {
    if(writeWallSurfaceAreaVTK_){
        writeVTKSurfaceArea(file, vtk);
    }
}

/*!
 * \brief Calculates and writes the surface areas of the cells to the vtu file.
 *
 * \details The calculation is based on the rendered VTKContainer and it takes time.
 * This functionality is diasabled by default and can be enabled using WallHandler::setWriteWallSurfaceAreaVTK().
 *
 * \param[in] file Output filestream
 * \param[in] vtk Rendered VTKContainer
 */
void WallVTKWriter::writeVTKSurfaceArea(std::fstream& file, VTKContainer& vtk) const {
    file << "  <DataArray type=\"Float32\"  Name=\"SurfaceArea\" format=\"ascii\">\n";
    /**
     * \todo gmb This only needs to be calculated if the mesh changes.
     */
    for (const std::vector<double>& strip : vtk.triangleStrips){
        double area = 0;
        for(size_t i = 1; i < strip.size() - 1; i++){
            size_t i1 = static_cast<size_t>(strip.at(i-1));
            size_t i2 = static_cast<size_t>(strip.at(i));
            size_t i3 = static_cast<size_t>(strip.at(i+1));

            /**
             * \todo gmb Check if it always true.
             */
            // strip is a polygon
            if(strip.front() == strip.back() && i == strip.size() - 2){
                break;
            }

            Vec3D &p1 = vtk.points.at(i1);
            Vec3D &p2 = vtk.points.at(i2);
            Vec3D &p3 = vtk.points.at(i3);

            Vec3D p12 = p2 - p1;
            Vec3D p13 = p3 - p1;

            double triangleArea = 0.5 * Vec3D::cross(p12, p13).getLength();
            area += triangleArea;
        }
        file << '\t' << area << '\n';
    }

    file << "  </DataArray>\n";
}

void WallVTKWriter::setWriteWallSurfaceAreaVTK(bool writeWallSurfaceAreaVTK)
{
    writeWallSurfaceAreaVTK_ = writeWallSurfaceAreaVTK;
}

bool WallVTKWriter::getWriteWallSurfaceAreaVTK() const 
{
    return writeWallSurfaceAreaVTK_;
}
