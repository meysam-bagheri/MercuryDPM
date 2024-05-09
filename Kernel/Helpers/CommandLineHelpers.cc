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

#include "Helpers/CommandLineHelpers.h"

/*!
 * \brief Counts the leading dash ('-') characters in a string.
 * \param[in] s String constaining dash ('-') characters.
 * \returns Number of leading dashes ('-')
 */
std::size_t helpers::countLeadingDashes(const std::string & s){
    std::size_t counter = 0;

    for(char c : s){
        if(c == '-'){
            ++counter;
        } else {
            break;
        }
    }

    return counter;
}

/*!
 * \param[in] argc
 * \param[in] argv
 * \param[in] varName The name of the commandline argument to be removed
 * \param[in] nArgs The number of additional command line arguments to remove
 *
 * \returns boolean True if the argument was found and removed
 * \details This function does hide the desired argument from the supplied argv
 * argc combination. It does it by moving the specified arguent to the end of the
 * supplied argv and reducing argc by the coorect number. Other pieces of code that
 * rely on argc should therefor no longer see the hidden argument.
 *
 * If used in comination with readFromCommandLine(...), this function allows handling
 * of arguments that are not seen by solve(), even if commandline arguments are passed.
 */
bool helpers::removeFromCommandline(int& argc, char* argv[], const std::string & varName, int nArgs)
{
    for (int i = 0; i < argc; ++i) {
        if (varName == argv[i]) 
        {
            char *tmp[nArgs + 1];
            
            for (int j = i; j < argc - 1 - nArgs; j++)
            {
                // Store the pointers of the handled argument
                if (j < i + 1 + nArgs)
                {
                    tmp[j-i] = argv[j];
                }
                
                // Move the pointers after the argument to the front
                argv[j] = argv[j + nArgs + 1];
            }
            
            // Move the stored argument to the end
            for (int j = argc - 1 - nArgs; j < argc; j++)
            {
                argv[j] = tmp[j + 1 + nArgs - argc];
            }
            argc -= nArgs + 1;
            
            return true;
        }
    }
    
    return false;
}

template<>
std::string helpers::readFromCommandLine<std::string>(int argc, char *argv[], const std::string & varName, std::string value)
{
    const size_t nDashes = countLeadingDashes(varName);
    for (int i = 0; i < argc - 1; ++i) {
        if (varName == argv[i]) {
            value = argv[i+1];
            logger(INFO, "readFromCommandLine: % set to % ", varName.substr(nDashes), value);
            return value;
        }
    }
    //if the variable is not found
    logger(INFO, "readFromCommandLine: % set to default value % ", varName.substr(nDashes), value);
    return value;
}

/**
 * Returns true if command line arguments contain varName, false else
 * Usage example:
 *    if (readFromCommandLine(argc, argv, '-verbose')) ...
 * @param argc pass through number of command line arguments
 * @param argv pass through values of command line arguments
 * @param varName name of command line arguments that is required to return true
 * @return true or false
 */
bool helpers::readFromCommandLine(int argc, char *argv[], const std::string & varName)
{
    const size_t nDashes = countLeadingDashes(varName);
    for (int i = 0; i < argc; ++i) {
        if (varName == argv[i]) {
            logger(INFO, "readFromCommandLine: % set to true", varName.substr(nDashes));
            return true;
        }
    }
    //if the variable is not found
    logger(INFO, "readFromCommandLine: % set to default value false", varName.substr(nDashes));
    return false;
}
