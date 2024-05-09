//Copyright (c) 2013-2024, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef MERCURYDPM_PRINTWALLTIMEMIXIN_H
#define MERCURYDPM_PRINTWALLTIMEMIXIN_H

#include <chrono>
#include "DPMBase.h"

/*!
 * \class PrintWallTimeMixin
 * \brief A mixin class for printing wall time during simulations.
 *
 * \details This class provides functionality to print the wall time (real-world time) at different stages of a simulation.
 * It also provides methods to estimate the remaining wall time based on the current progress of the simulation.
 *
 * \note This class is a mixin and should be used as part of multiple inheritance.
 */
class PrintWallTimeMixin : virtual public DPMBase {
public:

    /*!
     * \brief Default constructor.
     */
    PrintWallTimeMixin();

    /*!
     * \brief Prints the current wall time.
     *
     * This function is overridden from DPMBase.
     */
    void printTime() const override;

private:
    /*!
     * \brief Initial simulation time.
     */
    Mdouble initialSimulationTime;

    /*!
     * \brief wall time (imagine a clock on the wall)
     */
    std::chrono::system_clock::time_point initialWallTime;

    /*!
     * \brief Last simulation time. Mutable so it can be modified in the const method printTime.
     */
    mutable Mdouble lastSimulationTime;

    /*!
     * \brief Last wall time. Mutable so it can be modified in the const method printTime.
     */
    mutable std::chrono::system_clock::time_point lastWallTime;

    /*!
     * \brief Formats a time point for display.
     *
     * \param[in] tp The time point to format.
     * \return A string representation of the time point.
     */
    static std::string formatTime(const std::chrono::system_clock::time_point tp);

    /*!
     * \brief Formats a duration for display.
     *
     * \param[in] duration The duration to format.
     * \return A string representation of the duration.
     */
    static std::string formatDuration(const std::chrono::system_clock::duration duration);

    /*!
     * \brief Gets the current wall time.
     *
     * \return The current wall time.
     */
    std::chrono::system_clock::time_point getCurrentWallTime() const;

    /*!
     * \brief Estimates the remaining wall time for the simulation.
     *
     * \return The estimated remaining wall time.
     */
    std::chrono::system_clock::duration estimateRemainingWallTime() const;

    /*!
     * \brief Estimates the remaining wall time for the simulation using dead reckoning.
     *
     * \return The estimated remaining wall time.
     */
    std::chrono::system_clock::duration deadReckoningEstimateRemainingWallTime() const;
};

#endif //MERCURYDPM_PRINTWALLTIMEMIXIN_H