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

#include "WallDetailsVTKWriter.h"
#include "DPMBase.h"
#include "Walls/NurbsWall.h"
#include "Walls/WearableNurbsWall.h"

void WallDetailsVTKWriter::writeVTK() const
{
    if (PROCESSOR_ID != 0) return;
    
    // For each type of wall a VTKData object needs to be created, as each type is written to a different vtu file
    VTKData dataNurbsWall;
    VTKData dataWearableNurbsWall;
    //VTKData dataOtherTypeOfWall;
    
    bool nw_sw = shouldWrite(WallHandler::DetailsVTKOptions::NURBSWALL);
    bool wnw_sw = shouldWrite(WallHandler::DetailsVTKOptions::WEARABLENURBSWALL);
    
    // Loop through all wall and check if they belong to a certain type and if the data should be written or not
    for (const auto& w : handler_)
    {
        // WearableNurbsWall inherits from NurbsWall, therefore make sure to differentiate between them.
        // To prevent many of the same castings, first check if any of the two needs to be written.
        if (nw_sw || wnw_sw)
        {
            // Just a single cast to WearableNurbsWall and NurbsWall
            auto wnw_ptr = dynamic_cast<WearableNurbsWall*>(w);
            // NurbsWall: should be written and correct type (a valid cast is possible)
            if (nw_sw && dynamic_cast<NurbsWall*>(w))
            {
                // Make sure the NurbsWall vtk data is written, not the WearableNurbsWall implementation of writeWallDetailsVTK()
                if (wnw_ptr)
                    wnw_ptr->NurbsWall::writeWallDetailsVTK(dataNurbsWall);
                else
                    w->writeWallDetailsVTK(dataNurbsWall);
            }
            // WearableNurbsWall: should be written and of type
            if (wnw_sw && wnw_ptr)
            {
                w->writeWallDetailsVTK(dataWearableNurbsWall);
            }
        }
        // Without inheritance differentiating a simple call could look like:
        //else if (shouldWrite(WallHandler::DetailsVTKOptions::OTHER_TYPE_OF_WALL) && dynamic_cast<OtherTypeOfWall*>(w))
        //{
        //    w->writeWallDetailsVTK(dataOtherTypeOfWall);
        //}
    }
    
    // Again, check if the data should be written and if so write the file(s).
    // When a certain wall type should be written, but none of the walls the handler are of that type, the file(s) will
    // be written anyway. This is closest to what a user would expect, even though the file(s) don't hold any data.
    if (nw_sw)
        dataNurbsWall.writeVTKData(generateFileName("NurbsWall"));
    if (wnw_sw)
        dataWearableNurbsWall.writeVTKData(generateFileName("WearableNurbsWall"));
    //if (shouldWrite(WallHandler::DetailsVTKOptions::OTHER_TYPE_OF_WALL))
    //    dataOtherTypeOfWall.writeVTKData(generateFileName("OtherTypeOfWall"));
    
    // Other stuff not strictly wall specific
    if (shouldWrite(WallHandler::DetailsVTKOptions::BOUNDINGBOX))
    {
        VTKData dataBoundingBox;
        handler_.writeWallDetailsVTKBoundingBox(dataBoundingBox);
        dataBoundingBox.writeVTKData(generateFileName("BoundingBox"));
    }
    
    fileCounter++;
}

bool WallDetailsVTKWriter::shouldWrite(WallHandler::DetailsVTKOptions type) const
{
    FileType fileType = handler_.getWriteDetailsVTK(type);
    return (fileType == FileType::ONE_FILE && fileCounter == 0) || fileType == FileType::MULTIPLE_FILES || fileType == FileType::MULTIPLE_FILES_PADDED;
}

std::string WallDetailsVTKWriter::generateFileName(const std::string& identifier) const
{
    return outputDirectory_ + "/" + handler_.getDPMBase()->getName() + "WallDetails" + identifier + "_" + std::to_string(fileCounter) + ".vtu";
}