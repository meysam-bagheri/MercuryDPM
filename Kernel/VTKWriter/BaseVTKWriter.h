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


#ifndef VTKWRITER_H
#define VTKWRITER_H

#include <string>
#include <fstream>
#include "Logger.h"

/**
 * 
 * @tparam H some handler, e.g. particleHandler
 */
template<typename H>
class BaseVTKWriter
{

public:
    
    BaseVTKWriter(H& handler) : handler_(handler)
    {
        fileCounter = 0;
        outputDirectory_ = ".";
    }
    
    BaseVTKWriter(const BaseVTKWriter& other) : handler_(other.handler_)
    {
        fileCounter = other.fileCounter;
        outputDirectory_ = other.outputDirectory_;
    }
    
    virtual void writeVTK() const = 0;
    
    unsigned getFileCounter() const
    {
        return fileCounter;
    }
    
    void setFileCounter(unsigned fileCounter)
    {
        this->fileCounter = fileCounter;
    }

    void setOutputDirectory(const std::string & dir)
    {
        logger.assert_always(!dir.empty(), "Invalid VTK directory");

        if (dir.at(dir.size() - 1) == '/'){
            outputDirectory_ = dir.substr(0, dir.size() - 1);
        } else {
            outputDirectory_ = dir;
        }
    }

    const std::string & getOutputDirectory() const
    {
        return outputDirectory_;
    } 
    
protected:
    std::fstream makeVTKFileWithHeader() const;
    
    void writeVTKFooterAndClose(std::fstream& file) const;
    
    ///particle handler from which the particles should be written
    H& handler_;
    
    mutable unsigned int fileCounter;

    std::string outputDirectory_;
};

///\todo vtw wall files only need to be written by one processor
template<typename T>
std::fstream BaseVTKWriter<T>::makeVTKFileWithHeader() const
{
    //extract the word "Wall" or "Particle" from the VTK writer name
    std::string name = handler_.getName();
    name = name.substr(0, name.length() - 7);

    //determine name of output file
    std::string fileName;
#ifdef MERCURYDPM_USE_MPI
    if (NUMBER_OF_PROCESSORS > 1 && name != "Wall")
    {
        fileName = handler_.getDPMBase()->getName() + "Processor_" + std::to_string(PROCESSOR_ID) +
        '_' + name + '_' + std::to_string(fileCounter++)  + ".vtu";
    }
    else
    {
        fileName = handler_.getDPMBase()->getName() +
                                 name + '_' +
                                 std::to_string(fileCounter++) + ".vtu";
    }
#else
    fileName = handler_.getDPMBase()->getName() +
               name + '_' +
               std::to_string(fileCounter++) + ".vtu";
#endif

    // open output file
    char filePath[2048];
    const int n = snprintf(filePath, sizeof(filePath), "%s/%s", outputDirectory_.c_str(), fileName.c_str());
    if (n <= 0) {
        logger(ERROR, "unable to create file path");
    }

    std::fstream file;
    file.open(filePath, std::ios_base::out);
    if (file.fail())
    {
        logger(WARN, "File % could not be opened", filePath);
    }

    // write output file header
    file << "<?xml version=\"1.0\"?>\n";
    file << "<!-- time " << handler_.getDPMBase()->getTime() << "-->\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "<UnstructuredGrid>\n";
    return file;
}

template<typename T>
void BaseVTKWriter<T>::writeVTKFooterAndClose(std::fstream& file) const
{
    // write output file footer
    file << "<Cells>\n";
    file << "  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    file << "  </DataArray>\n";
    file << "  <DataArray type=\"UInt8\"  Name=\"types\" format=\"ascii\">\n";
    file << "  </DataArray>\n";
    file << "</Cells>\n";
    file << "</Piece>\n";
    file << "</UnstructuredGrid>\n";
    file << "</VTKFile>\n";
    // close output file
    file.close();
}

#endif
