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

#ifndef MERCURYDPM_STRING_HELPERS_H
#define MERCURYDPM_STRING_HELPERS_H

#include "Math/ExtendedMath.h"

#include <string>
#include <memory>

namespace helpers
{
    /*!
     * \brief returns the input string after converting upper-case characters to lower case
     */
    std::string lower(std::string s);

    /*!
     * \brief Reads a line from one stringstream into another, and prepares the latter for reading in.
     */
    void getLineFromStringStream(std::istream& in, std::stringstream& out);

    bool isNext(std::istream& is, const std::string name);

    /**
     * \brief Checks if the next argument in the input stream is a certain string
     */
    bool compare(std::istream& is, std::string s);

    /*!
     * \brief converts a floating point number into a string with a given precision
     */
    std::string toString(Mdouble value, unsigned precision);

    template<typename T>
    std::string toString(const T& n)
    {
        std::ostringstream stm;
        stm << n;
        return stm.str();
    }

    template<typename T>
    std::string toString(const std::vector<T>& vec)
    {
        std::ostringstream stm;
        for (const auto val : vec) {
            stm << val << ' ';
        }
        return stm.str();
    }

    template <typename... Args>
    std::string stringFormat(const std::string& format, Args... args)
    {
        int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) + 1; // Extra space for '\0'
        if (size_s <= 0) {
            logger(ERROR, "snprintf() failed");
        }
        auto size = static_cast<size_t>(size_s);
        std::unique_ptr<char[]> buf(new char[size]);
        std::snprintf(buf.get(), size, format.c_str(), args...);
        return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
    }
}

#endif // MERCURYDPM_STRING_HELPERS_H
