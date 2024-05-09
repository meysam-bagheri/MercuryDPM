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

#include "Helpers/MathHelpers.h"
#include "Helpers/StringHelpers.h"

#include <algorithm>
#include <cctype>

std::string helpers::lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

std::string helpers::toString(Mdouble value, unsigned precision)
{
    std::ostringstream stm;
    stm << round(value, precision);
    return stm.str();
}

/*!
 * \details This function is used to avoid errors from reading in old or manually modified restart files.
 * Instead of reading variable by variable directly from the restart stringstream,
 * a full line is read first, from which the variables are read. Thus, if a line
 * has the wrong number of arguments, it might affect the reading of the current
 * line, but correctly reads the next line.
 *
 * Example of usage:
 * > std::stringstream line;
 * > std::stringstream is = restartFile.getFStream();
 * > helpers::getLineFromStringStream(is, line);
 * > std::string dummy;
 * > line >> dummy;
 *
 * \param[in]  in the stringstream from which a line is read out should be initialized as std::stringstream(std::stringstream::out)
 * \param[out] out the stringstream into which the line is read; should be initialized as std::stringstream(std::stringstream::in | std::stringstream::out)
 */
void helpers::getLineFromStringStream(std::istream& in, std::stringstream& out)
{
    std::string line_string;
    getline(in, line_string);
    out.str(std::move(line_string));
    out.clear();
}

/**
 * reads next value in stream as a string and compares it with name.
 * If name is equal to string, the function outputs true.
 * If name is not equal to string, the function undoes the read by setting seekg, and outputs false.
 * @param is
 * @param name
 * @return
 */
bool helpers::isNext(std::istream& is, const std::string name) {
    std::string dummy;
    auto pos = is.tellg();
    is >> dummy;
    if (dummy != name) {
        is.seekg(pos);
        return false;
    } else {
        return true;
    }
}

bool helpers::compare(std::istream& is, std::string s)
{
    // Get current position
    //check if the next line starts with 'interactionFile'; otherwise, skip interaction
    int len = is.tellg();
    std::string dummy;
    is >> dummy;
    if (dummy != s)
    {
        is.seekg(len, std::ios_base::beg);
        logger(VERBOSE, "helpers::compare: Next stream value (%) is not %", dummy, s);
        return false;
    }
    return true;
}
