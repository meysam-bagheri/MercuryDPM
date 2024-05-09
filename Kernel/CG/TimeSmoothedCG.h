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
#ifndef TimeSmoothedCG_H
#define TimeSmoothedCG_H

#include <GeneralDefine.h>
#include "CG/CG.h"
#include <iostream>
#include <vector>
#include "TimeSmoothedFields.h"

class BaseParticle;

class BaseInteraction;

class DPMBase;

class Function;

/*!
 * \brief Evaluates time-smoothed continuum fields and writes the data into a 
 * stat file.
 * \details
 * Like CG, this class should be used to evaluate time-dependent fields. 
 * However, as the data is averaged over several time steps, the resulting 
 * fields are smoother and more reliable, as they use more information.
 * For time-independent fields, see TimeAveragedCG.
 * 
 * The class is derived from CG, which already contains the basic 
 * functionality to obtain time-dependent statistics. This class then averages 
 * the field values F over time, with a non-uniform weight for each time step,
 * \f[\bar F(t,w_t) = \frac{\sum_i w(t-t_i) F(t_i)}{\sum_i w(t-t_i)},\f]
 * with the weight given by a Gaussian distribution,
 * \f[\bar w(t-t_i) = exp(-(t-t_i)^2/2w_t^2).\f]
 */
template<class Coordinates, template<class> class BaseFunction=CGFunctions::Lucy, class Fields=CGFields::StandardFields>
class TimeSmoothedCG final : public CG<Coordinates, BaseFunction, Fields>
{
public:
    typedef BaseFunction<Coordinates> Function;
    
    /*!
     * \brief Default constructor; does nothing, i.e. no points are created
     * initially. 
     */
    TimeSmoothedCG();
    
    /*!
     * \brief Copy constructor. It copies the TimeSmoothedCGFunction and all objects it contains.
     * \param[in] p the TimeSmoothedCGFunction that has to be copied
     */
    TimeSmoothedCG(const TimeSmoothedCG& p) = default;
    
    /*!
     * \brief Destructor, it simply destructs the TimeSmoothedCGFunction and all the objects it contains.
     */
    virtual ~TimeSmoothedCG() = default;
    
    /*!
     * \brief
     */
    TimeSmoothedCG<Coordinates, BaseFunction, Fields>* copy() const override;
    
    /*!
     * \brief
     */
    void write(std::ostream& os) const override;
    
    void writeAll(std::ostream& os, TimeSmoothedFields<Fields>& average) const;
    
    /*!
     * \brief
     */
    std::string getName() const override;
    
    /*!
     * \brief Called at the beginning of the DPM simulation to initialise the cg 
     * evaluation and to open the statFile.
     */
    void initialise() override;
    
    /*!
     * \brief Called after a given number of time steps (statFile::saveCount_)
     * to evaluate the CG fields.
     */
    void evaluate() override;
    
    /*!
     * \brief Called at the end of the DPM simulation to finish the cg 
     * evaluation and to close the statFile.
     */
    void finish() override;
    
    void setWidthTime(Mdouble widthTime) final;
    
    Mdouble getWidthTime() const final;
    
    void setTimeStep(Mdouble timeStep) final;
    
    Mdouble getTimeStep() const final;

protected:
    
    /*!
     * Width of the Gauss function used to calculate the weights for the 
     * temporal smoothing.
     * Has to be specified by the user.
     */
    Mdouble widthTime_;
    
    /*!
     * Cutoff of the Gauss function used to calculate the weights for the 
     * temporal smoothing.
     * Internal variable, set to 3.0*widthTime_.
     */
    Mdouble cutoffTime_;
    
    /*!
     * Stores the next time at which smoothed fields should be evaluated, but 
     * which is not yet included in averages_. 
     * Internal variable, set DPMBase::time_ in initialize, then increased by 
     * timeStep_ every time a new object is added to averages_.
     */
    Mdouble nextTime_;
    
    /*!
     * Interval between two successive time steps for which time-smoothed fields 
     * are evaluated.
     * Has to be specified by the user.
     */
    Mdouble timeStep_;
    
    /*!
     * Stores the smoothed fields that are currently being evaluated.
     * Internal variable, see TimeSmoothedFields for details.
     */
    std::vector<TimeSmoothedFields<Fields> > averages_;
};

#include "TimeSmoothedCG.hcc"

#endif
