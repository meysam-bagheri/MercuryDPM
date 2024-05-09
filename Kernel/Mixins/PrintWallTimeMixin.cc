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

#include "PrintWallTimeMixin.h"

/*!
 * \brief Default constructor.
 *
 * \details Initializes the initial simulation time and wall time.
 */
PrintWallTimeMixin::PrintWallTimeMixin() {
    initialSimulationTime = getTime();
    initialWallTime = std::chrono::system_clock::now();
    lastWallTime = initialWallTime;
}

/*!
 * \brief Prints the current wall time.
 *
 * \details This function is overridden from DPMBase. It calculates the current simulation time and wall time,
 * calculates the time difference from the last wall time, and logs these details.
 */
void PrintWallTimeMixin::printTime() const {

#ifdef MERCURYDPM_USE_MPI
    MPIContainer& communicator = MPIContainer::Instance();
  if (communicator.getProcessorID() != 0)
    return;
#endif
    auto simulationTime = getTime();
    auto currentWallTime = getCurrentWallTime();
    auto timeDifference = currentWallTime - lastWallTime;
    logger(
            INFO, "t=%3.6 tmax=%3.6 wallTime=% dt=%3.6 timeLeft=%3.6 or %3.6",
            simulationTime, getTimeMax(),
            formatTime(currentWallTime),
            formatDuration(timeDifference),
            formatDuration(estimateRemainingWallTime()),
            formatDuration(deadReckoningEstimateRemainingWallTime())
    );

    // lastWallTime is mutable, so we can change it inside an otherwise
    // const function
    lastWallTime = currentWallTime;
    lastSimulationTime = simulationTime;
}

/*!
 * \brief Formats a time point for display.
 *
 * \details This function takes a time point and formats it into a string for display.
 *
 * \param[in] tp The time point to format.
 * \return A string representation of the time point.
 */
std::string PrintWallTimeMixin::formatTime(const std::chrono::system_clock::time_point tp) {
    auto t = std::chrono::system_clock::to_time_t(tp);
    std::ostringstream os;
    os << std::put_time(localtime(&t), "%Y-%m-%d %H:%M:%S");
    return os.str();
}

/*!
 * \brief Formats a duration for display.
 *
 * \details This function takes a duration and formats it into a string for display.
 *
 * \param[in] duration The duration to format.
 * \return A string representation of the duration.
 */
std::string PrintWallTimeMixin::formatDuration(const std::chrono::system_clock::duration duration) {
    auto dur = duration;  // not const

    std::ostringstream os;
    auto hours = std::chrono::duration_cast<std::chrono::hours>(dur).count();
    dur = dur % std::chrono::hours(1);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(dur).count();
    dur = dur % std::chrono::minutes(1);

    if (hours > 0) {
        os << hours << "h " << minutes << "m";
        return os.str();
    }

    if (minutes > 0) {
        auto seconds = std::chrono::duration_cast<std::chrono::seconds>(dur).count();
        os << minutes << "m " << seconds << "s";
        return os.str();
    }

    // in this case seconds is a double, not a long long int
    double seconds = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count() / 1000.;
    os << seconds << "s";
    return os.str();
}

/*!
 * \brief Gets the current wall time.
 *
 * \details This function returns the current wall time.
 *
 * \return The current wall time.
 */
std::chrono::system_clock::time_point PrintWallTimeMixin::getCurrentWallTime() const {
    return std::chrono::system_clock::now();
}

/*!
 * \brief Estimates the remaining wall time for the simulation.
 *
 * \details This function estimates the remaining wall time based on the current progress of the simulation.
 *
 * \return The estimated remaining wall time.
 */
std::chrono::system_clock::duration PrintWallTimeMixin::estimateRemainingWallTime() const {
    Mdouble time = getTime();
    Mdouble ratio = (getTimeMax() - time) / (time - initialSimulationTime);
    auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(
            getCurrentWallTime() - initialWallTime).count();
    long long int estimatedRemainingSeconds = elapsedSeconds * ratio;
    return std::chrono::seconds(estimatedRemainingSeconds);
}

/*!
 * \brief Estimates the remaining wall time for the simulation
 *
 * \details This function estimates the remaining wall time using a method called dead reckoning.
 *
 * \return The estimated remaining wall time.
 */
std::chrono::system_clock::duration PrintWallTimeMixin::deadReckoningEstimateRemainingWallTime() const {
    Mdouble time = getTime();
    Mdouble ratio = (getTimeMax() - time) / (time - lastSimulationTime);
    auto secondsSinceLastSave = std::chrono::duration_cast<std::chrono::seconds>(
            getCurrentWallTime() - lastWallTime).count();
    long long int estimatedRemainingSeconds = secondsSinceLastSave * ratio;
    return std::chrono::seconds(estimatedRemainingSeconds);
}