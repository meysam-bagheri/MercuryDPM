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

#ifndef MERCURYDPM_COMMANDLINE_HELPERS_H
#define MERCURYDPM_COMMANDLINE_HELPERS_H

#include "Math/ExtendedMath.h"

namespace helpers
{
    /*!
     * \brief Counts the leading dash ('-') characters in a string.
     */
    std::size_t countLeadingDashes(const std::string & s);

    /*!
     * \brief May be used to hide arguments from argc and argv.
     */
    bool removeFromCommandline(int& argc, char* argv[], const std::string & varName, int nArgs);

    /*!
     * \brief Returns true if command line arguments contain varName, false else
     */
    bool readFromCommandLine(int argc, char *argv[], const std::string & varName);

    template<typename T>
    T readFromCommandLine(int argc, char *argv[], const std::string & varName, T value)
    {
        const size_t nDashes = countLeadingDashes(varName);
        for (int i = 0; i < argc-1; ++i) {
            if (varName == argv[i]) {
                value = atof(argv[i+1]);
                logger(INFO, "readFromCommandLine: % set to % ", varName.substr(nDashes), value);
                return value;
            }
        }
        //if the variable is not found
        logger(INFO, "readFromCommandLine: % set to default value % ", varName.substr(nDashes), value);
        return value;
    }

    template<typename T, size_t n>
    std::array<T,n> readArrayFromCommandLine(int argc, char *argv[], const std::string & varName, std::array<T,n> value)
    {
        const size_t nDashes = countLeadingDashes(varName);
        for (int i = 0; i < argc-1; ++i) {
            if (varName == argv[i]) {
                int j = i + 1;
                std::stringstream out;
                for (auto& v : value) {
                    v = atof(argv[j]);
                    out << v << ' ';
                    ++j;
                }
                logger(INFO, "readArrayFromCommandLine: % set to % ", varName.substr(nDashes), out.str());
                return value;
            }
        }
        //if the variable is not found
        std::stringstream out;
        for (auto& v : value) out << v << ' ';
        logger(INFO, "readArrayFromCommandLine: % set to default value % ", varName.substr(nDashes), out.str());
        return value;
    }

    template<typename T>
    std::vector<T> readVectorFromCommandLine(int argc, char *argv[], const std::string & varName, size_t n, std::vector<T> values)
    {
        const size_t nDashes = countLeadingDashes(varName);
        for (int i = 0; i < argc-1; ++i) {
            if (varName == argv[i]) {
                // read until the next argument starts
                values.resize(0);
                std::stringstream out;
                for (int j = i+1; j < argc && argv[j][0] != '-'; ++j) {
                    values.push_back(atof(argv[j]));
                    out << values.back() << ' ';
                }
                logger(INFO, "readVectorFromCommandLine: % set to % ", varName.substr(nDashes), out.str());
                return values;
            }
        }
        //if the variable is not found
        std::stringstream out;
        for (auto& v: values) out << v << ' ';
        logger(INFO, "readVectorFromCommandLine: % set to default value % ", varName.substr(nDashes), out.str());
        return values;
    }

    template<>
    std::string readFromCommandLine<std::string>(int argc, char* argv[], const std::string & varName, std::string value);
}

#endif // MERCURYDPM_COMMANDLINE_HELPERS_H
