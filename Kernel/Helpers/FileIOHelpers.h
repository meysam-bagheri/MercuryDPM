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

#ifndef MERCURYDPM_FILEIO_HELPERS_H
#define MERCURYDPM_FILEIO_HELPERS_H

#include "Logger.h"
#include "Math/ExtendedMath.h"

#include <fstream>
#include <string>
#include <vector>

namespace helpers
{
    /*!
     * \brief Writes a string to a file.
     */
    bool writeToFile(const std::string & filename, const std::string & filecontent);

    /*!
     * \brief Writes a string to a file.
     */
    void writeCommandLineToFile(const std::string & filename, const int argc, char * const argv[]);

    /*!
     * \brief Adds a string to an existing file.
     */
    bool addToFile(const std::string & filename, const std::string & filecontent);

    /*!
    * \brief Function to check if a file exists, is used to check if a run has already need done.
    */
    bool fileExists(const std::string & strFilename);

    /*!
     * \brief Provides a simple interface for opening a file.
     */
    bool openFile(std::fstream& file, const std::string & filename, std::fstream::openmode mode);

    std::vector<double> readArrayFromFile(const std::string & filename, int& n, int& m);

    void more(const std::string & filename, unsigned nLines = constants::unsignedMax);

    /**
     * \brief Creates a directory.
     * 
     * \param directory Absoulte/relative path of the directory 
     * \param allowExists Won't fail if directory exists
     * \return Directory created
     */
    bool createDirectory(const std::string & directory, bool allowExists = false);

    std::string getPath();

    /**
     * \brief Reads optional variables in the restart file
     *
     * \details A variable is stored in the restart file by storing the variables name and value, e.g.
     *   " name value"
     * If a variable is always written to the restart file, it can be read-in like this:
     *   is >> dummy >> variable;
     * If a variable is optional, the variable name has to be checked, and should be read in like this:
     *   readOptionalVariable(is,"name",variable);
     */
    template<typename T>
    bool readOptionalVariable(std::istream& is, const std::string& name, T& variable)
    {
        ///\todo readOptionalVariable should check the full variable name, not just the next value.
        /// However, I don't know how to put the location in the ifstream back.
        const auto pos = is.tellg();
        std::string dummy;
        is >> dummy;
        if (dummy == name)
        {
            is >> variable;
            return true;
        } else {
            is.seekg(pos);
            return false;
        }
    }

    /**
     * Allows a quick read-in from a parameter file.
     *
     * For example, the following code reads in the variable time from the file in:
     * \code{.cpp}
     *   double time = readFromFile("in","time",24);
     * \endcode
     *
     * The in file needs to contain the string time, followed by the value, e.g.
     * \code{.sh}
     *   time 20
     * \endcode
     *
     * If the file cannot be opened, or the parameter is not found,
     * the variable is set to the default value specified by the third argument.
     *
     * @param fileName name of input
     * @param varName  variable name as it appears in the input file
     * @param defaultValue    default value (used if the parameter could not be read)
     * @return         value of variable
     */
    template<typename T>
    T readFromFile(const std::string & fileName, const std::string & varName, const T defaultValue)
    {
        //open filestream
        std::ifstream is(fileName.c_str(), std::ios::in);
        if (is.fail())
        {
            logger(INFO, "readFromFile: file % could not be opened, variable % set to default value %",
                fileName, varName, defaultValue);
            return defaultValue;
        }
        
        //read in variables, until the right one is fount
        std::string s;
        while (!is.eof())
        {
            is >> s;
            if (s == varName)
            {
                T value;
                is >> value;
                logger(INFO, "readFromFile: variable % set to % ", varName, value);
                return value;
            }
        }
        
        //if the right variable is never found
        logger(WARN, "readFromFile: variable % not set in file %, using default value % ", varName, fileName, defaultValue);
        return defaultValue;
    }
}

#endif // MERCURYDPM_FILEIO_HELPERS_H
