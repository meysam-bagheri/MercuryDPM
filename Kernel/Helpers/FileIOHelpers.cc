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

#include "Helpers/FileIOHelpers.h"
#include "Logger.h"

#include <errno.h>
#include <fstream>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

/*!
 * \details Provides a simple interface for writing a string to a file.
 * This function is mainly used to create ini or restart file that the code
 * later reads back in.
 *
 * Example of usage:
 * > helpers::writeToFile("RestartUnitTest.ini",
 * > "1 0 0 0 0 1 1 1\n"
 * > "0.5 0.5 0  0 0 0.5  0  0 0 0  0 0 0  0\n");
 *
 * \param[in] filename the name of the file
 * \param[in] filecontent the content
 * \returns true on success.
 * \todo gmb Make this MPI compatible.
 */
bool helpers::writeToFile(const std::string & filename, const std::string & filecontent)
{
    std::fstream file;
    file.open(filename.c_str(), std::ios::out);
    if (file.fail())
    {
        logger(WARN, "Error in writeToFile: file could not be opened");
        return false;
    }
    file << filecontent;
    file.close();
    return true;
}

/*!
 * \todo gmb Make this MPI compatible.
 * \note Not used function.
 */
void helpers::writeCommandLineToFile(const std::string & filename, const int argc, char * const argv[])
{
    std::stringstream ss;
    for (int i=0; i<argc; ++i) {
        ss << argv[i] << ' ';
    }
    writeToFile(filename,ss.str());
}

/*!
 * \todo gmb Make this MPI compatible.
 */
bool helpers::addToFile(const std::string & filename, const std::string & filecontent)
{
    std::fstream file;
    file.open(filename.c_str(), std::ios::app);
    if (file.fail())
    {
        logger(INFO, "Error in writeToFile: file could not be opened");
        return false;
    }
    file << filecontent;
    file.close();
    return true;
}

/*!
 * \details This is a FileExist routine, which is used to test if a run have
 * already need preformed, allows me to plug holes in parm studies.
 */
bool helpers::fileExists(const std::string & strFilename)
{
    struct stat stFileInfo;
    bool blnReturn;
    int intStat;

    // Attempt to get the file attributes

    intStat = stat(strFilename.c_str(), &stFileInfo);
    if (intStat == 0)
    {
        // We were able to get the file attributes
        // so the file obviously exists.
        blnReturn = true;
    }
    else
    {
        // We were not able to get the file attributes.
        // This may mean that we don't have permission to
        // access the folder which contains this file. If you
        // need to do that level of checking, lookup the
        // return values of stat which will give you
        // more details on why stat failed.
        blnReturn = false;
    }

    return blnReturn;
}

/*!
 * \details Provides a simple interface for opening a file, in order to avoid
 * that the user has to learn the syntax for opening a file.
 * \param[out] file The std::fstream object that the user can write to.
 * \param[in] filename The name of the file.
 * \param[in] mode The openmode of the file, typically std::fstream::out or std::fstream::in.
 * \return true is the file was successfully opened, false else.
 * \todo gmb Make this MPI compatible. Is it used?
 */
bool helpers::openFile(std::fstream& file, const std::string & filename, std::fstream::openmode mode)
{
    file.open(filename.c_str(), mode);
    if (file.fail())
        return false;
    else
        return true;
}

/*!
 * \todo gmb Make this MPI compatible.
 * \note Not used function.
 */
std::vector<double> helpers::readArrayFromFile(const std::string & filename, int& n, int& m)
{
    std::fstream file;
    file.open(filename.c_str(), std::ios::in);
    if (file.fail())
    {
        logger(ERROR, "Error in readArrayFromFile: file could not be opened");
    }
    file >> n >> m;
    std::vector<double> v;
    Mdouble val;
    for (int i = 0; i < n * m; i++)
    {
        file >> val;
        v.push_back(val);
    }
    file.close();
    return v;
}

/*!
 * \todo gmb Make this MPI compatible.
 */
void helpers::more(const std::string & filename, unsigned nLines)
{
    if (nLines != constants::unsignedMax)
        logger(INFO, "First % lines of %:\n", Flusher::NO_FLUSH, nLines, filename);
    std::fstream file;
    file.open(filename.c_str(), std::ios::in);
    if (file.fail())
        logger(ERROR, "Error in readArrayFromFile: file could not be opened");
    std::string line;
    for (unsigned i = 0; i < nLines; i++)
    {
        if (file.eof()) break;
        getline(file, line);
        logger(INFO, " %\n", line);
    }
    file.close();
}


/**
 * \brief Creates a directory.
 * 
 * \param directory Absoulte/relative path of the directory
 * \param allowExists Won't fail if directory exists
 * \return true/false
 */
bool helpers::createDirectory(const std::string & directory, bool allowExists) 
{
#ifdef MERCURYDPM_USE_MPI
    if (PROCESSOR_ID == 0)
    {
#endif

    if (directory == ".") {
        return true;
    }

    errno = 0;
    const int err = ::mkdir(directory.c_str(), 0777);
    if (err == 0) {
        return true;
    } else if(errno == EEXIST && allowExists) {
        return true;
    } else {
        logger(ERROR, "Unable to create directory `%`: % (%)", directory, strerror(errno), errno);
        return false;
    }

#ifdef MERCURYDPM_USE_MPI
    } else {
        return true;
    }
#endif
}

std::string helpers::getPath()
{
    char currentPath[FILENAME_MAX];
    if (getcwd(currentPath, sizeof(currentPath)) == NULL) {
        logger(WARN, "Get current dir failed: %", strerror(errno));
    }
    return std::string(currentPath);
}
